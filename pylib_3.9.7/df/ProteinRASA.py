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

from df.bio_helper import column_to_sequences, string_to_sequence, column_to_structures
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

        id_column = input_field_to_column(request, 'uiIDColumn')

        structure_column = input_field_to_column(request, 'uiStructureColumn')
        structures = column_to_structures(structure_column, id_column)

        sequence_column = input_field_to_column(request, 'uiSequenceColumn')
        sequences = column_to_sequences(sequence_column, id_column)

        rasa_cutoff = integer_input_field(request, 'uiRASACutoff')

        # setup BioPython tools
        sr = ShrakeRupley()
        parser = PDBParser(QUIET = True)
        backbone_atoms = ['N', 'H', 'CA', 'HA', 'C', 'O', 'H2', 'H3']

        # output storage
        notifications = []
        output_sequences = []

        # process each structure
        for index, (structure, sequence) in enumerate(zip(structures, sequences)):
            if id_column:
                identifier = id_column.values[index]
            else:
                identifier = index

            if not (structure and sequence):
                notifications.append(Notification(level = NotificationLevel.ERROR,
                                                  title = 'Relative Accessible Surface Area',
                                                  summary = f'Error for structure {identifier}/nNull values found.',
                                                  details = f'Either the structure or sequence were null.'))

            with StringIO(structures[index]) as structure_fh:
                structure = parser.get_structure(identifier, structure_fh)

            # compute RASA for each residue/atom in context of the intact protein
            sr.compute(structure, level = 'A')

            for residue_index, residue in enumerate(structure.get_residues()):
                # compute RASA for each residue isolated from the protein structure
                residue_copy = residue.copy()
                sr.compute(residue_copy, level = 'A')

                in_context_sa = sum([atom.sasa for atom in residue.get_atoms() if atom.get_id() not in backbone_atoms])
                reference_sa = sum([atom.sasa for atom in residue_copy.get_atoms() if atom.get_id() not in backbone_atoms])

                # create annotations
                feature = SeqFeature(FeatureLocation(residue_index + 1, residue_index + 1),
                                     type = 'misc_feature',
                                     qualifiers = {'note': ['glysade_annotation_type: RASA',
                                                            f'RASA: {round(in_context_sa / reference_sa, 2) / 100}%',
                                                            f'Residue ID:  {residue.resname}',
                                                            f'"ld_style:{{"color": "{color}", "shape": "rounded-rectangle"}}"']})
                sequence.features.append(feature)

            output_sequences.append(sequence)

        # contruct output column
        rows = [sequence_to_genbank_base64_str(s) for s in output_sequences]
        output_column = ColumnData(name = f'Relative Accessible Surface Area Annotations',
                                   dataType = DataType.BINARY,
                                   contentType = 'chemical/x-genbank', values = rows)
        return DataFunctionResponse(outputColumns=[output_column])

        #     except Exception as ex:
        #         notifications.append(Notification(level = NotificationLevel.ERROR,
        #                                           title = 'Antibody Structure Prediction',
        #                                           summary = f'Error for ID {ab_id}/n{ex.__class__} - {ex}',
        #                                           details = f'{traceback.format_exc()}'))
        #