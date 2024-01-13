import base64
import gzip
from io import StringIO
import os
import traceback
from typing import Optional

from Bio.SeqRecord import SeqRecord
from Bio.PDB import PDBParser
from Bio.PDB.SASA import ShrakeRupley
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqUtils import seq1 as amino3to1

from df.bio_helper import column_to_sequences, string_to_sequence, column_to_structures, generate_color_gradient
from df.data_transfer import ColumnData, TableData, DataFunctionRequest, DataFunctionResponse, DataFunction, DataType, \
                             string_input_field, boolean_input_field, integer_input_field,  input_field_to_column, \
                             Notification, NotificationLevel

from ruse.bio.antibody import align_antibody_sequences, NumberingScheme, CDRDefinitionScheme, ANTIBODY_NUMBERING_COLUMN_PROPERTY
from ruse.bio.bio_data_table_helper import sequence_to_genbank_base64_str

class ProteinRASA(DataFunction):
    """
    Computes relative solvent accessible surface area of a given protein
    """

    def execute(self, request: DataFunctionRequest) -> DataFunctionResponse:

        # get settings from UI
        structure_column = input_field_to_column(request, 'uiStructureColumn')
        structures, notifications = column_to_structures(structure_column)

        sequence_column = input_field_to_column(request, 'uiSequenceColumn')
        sequences = column_to_sequences(sequence_column)

        rasa_cutoff = integer_input_field(request, 'uiRASACutoff')
        exposedColor = string_input_field(request, 'uiExposedColor')

        labelAll = boolean_input_field(request, 'uiLabelAll')
        startColor = string_input_field(request, 'uiStartColor')
        endColor = string_input_field(request, 'uiEndColor')

        if labelAll:
            RASA_color_gradient = generate_color_gradient(startColor, endColor, 21)

        # setup BioPython tools
        sr = ShrakeRupley()
        backbone_atoms = ['N', 'H', 'CA', 'HA', 'C', 'O', 'H2', 'H3']

        # output storage
        output_sequences = []

        # process each structure
        for index, (structure, sequence) in enumerate(zip(structures, sequences)):
            if not (structure and sequence):
                notifications.append(Notification(level = NotificationLevel.ERROR,
                                                  title = 'Relative Accessible Surface Area',
                                                  summary = f'Error for structure or sequence in Row {index}/nNull values found.',
                                                  details = f'Either the structure or sequence were null.'))
                continue

            # create a mapping between the structure residues and potentially gapped sequence residues
            # would this be a useful biotools function rather than buried here?
            sequence_number_mapping = {}
            structure_residue_number = 0
            structure_residues = [amino3to1(res.get_resname()) for res in structure.get_residues()]
            for sequence_residue_number, sequence_residue in enumerate(sequence.seq):
                if sequence_residue != '-':
                    if sequence_residue != structure_residues[structure_residue_number]:
                        # send a notification and then fail
                        pass
                    sequence_number_mapping[structure_residue_number] = sequence_residue_number
                    structure_residue_number += 1

            # compute RASA for each residue/atom in context of the intact protein
            try:
                sr.compute(structure, level = 'A')
            except Exception as ex:
                notifications.append(Notification(level = NotificationLevel.ERROR,
                                                  title = 'Relative Accessible Surface Area',
                                                  summary = f'Error for structure in Row {index}.\n' +
                                                            'Shrake-Rupley calculation failed for complete structure.',
                                                  details = f'{ex.__class__} - {ex}\n' +
                                                            f'{traceback.format_exc()}'))
                output_sequences.append(sequence)
                continue

            success = True  # for scoping outside loop
            for residue_index, residue in enumerate(structure.get_residues()):
                # compute RASA for each residue isolated from the protein structure
                residue_copy = residue.copy()

                success = True  # optimistic for success
                ex_storage = None  # use for possible storage of exception for use outside loop
                try:
                    sr.compute(residue_copy, level = 'A')
                except Exception as ex:
                    ex_storage = ex.copy()
                    success = False
                    continue

                in_context_sa = sum([atom.sasa for atom in residue.get_atoms() if atom.get_id() not in backbone_atoms])
                reference_sa = sum([atom.sasa for atom in residue_copy.get_atoms() if atom.get_id() not in backbone_atoms])
                rasa = round(in_context_sa / reference_sa * 100, 1)

                # create annotations
                if labelAll:
                    feature = SeqFeature(FeatureLocation(sequence_number_mapping[residue_index],
                                                         sequence_number_mapping[residue_index] + 1),
                                         type = 'misc_feature',
                                         qualifiers = {'feature_name': 'Relative Accessible Surface Area',
                                                       'note': ['RASA',
                                                                'glysade_annotation_type: RASA',
                                                                f'RASA: {rasa}%',
                                                                f'Residue ID:  {residue.resname}',
                                                                f'ld_style:{{"color": "{RASA_color_gradient[int(rasa // 5)]}", "shape": "rounded-rectangle"}}',
                                                                'ld_track: RASA']})
                    sequence.features.append(feature)

                if rasa >= rasa_cutoff:
                    feature = SeqFeature(FeatureLocation(sequence_number_mapping[residue_index],
                                                         sequence_number_mapping[residue_index] + 1),
                                         type = 'misc_feature',
                                         qualifiers = {'feature_name': 'Solvent Exposed Residue',
                                                       'note': [f'Exposed >= {rasa_cutoff}',
                                                                'glysade_annotation_type: RASA_exposed',
                                                                f'RASA: {rasa}%',
                                                                f'Residue ID:  {residue.resname}',
                                                                f'ld_style:{{"color": "{exposedColor}", "shape": "rounded-rectangle"}}',
                                                                'ld_track: RASA_exposed']})
                    sequence.features.append(feature)

            if not success:
                notifications.append(Notification(level = NotificationLevel.ERROR,
                                                  title = 'Relative Accessible Surface Area',
                                                  summary = f'Error for structure in Row {index}.\n' +
                                                            f'Shrake-Rupley calculation failed for for some residues.',
                                                  details = f'Example error: {ex.__class__} - {ex}\n' +
                                                            f'{traceback.format_exc()}'))
            output_sequences.append(sequence)

        # construct output column
        rows = [sequence_to_genbank_base64_str(s) for s in output_sequences]
        output_column = ColumnData(name = f'Relative Accessible Surface Area Annotations',
                                   dataType = DataType.BINARY,
                                   contentType = 'chemical/x-genbank', values = rows)

        return DataFunctionResponse(outputColumns = [output_column])
