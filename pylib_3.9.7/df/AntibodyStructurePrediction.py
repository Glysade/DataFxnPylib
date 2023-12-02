import base64
import os
import gzip
from typing import Optional

from Bio.SeqRecord import SeqRecord

from df.bio_helper import column_to_sequences, string_to_sequence
from df.data_transfer import ColumnData, TableData, DataFunctionRequest, DataFunctionResponse, DataFunction, DataType, \
                             string_input_field, boolean_input_field, input_field_to_column, \
                             Notification, NotificationLevel

from ruse.bio.antibody import align_antibody_sequences, NumberingScheme, CDRDefinitionScheme, ANTIBODY_NUMBERING_COLUMN_PROPERTY
from ruse.bio.bio_data_table_helper import sequence_to_genbank_base64_str

from ImmuneBuilder import ABodyBuilder2
from ImmuneBuilder.ABodyBuilder2 import header as ABB2_HEADER
from ImmuneBuilder.refine import refine
from ImmuneBuilder.util import add_errors_as_bfactors

class AntibodyStructurePrediction(DataFunction):
    """
    Predicts antibody structure from heavy and light chain sequences
    """

    def _antibody_numbering(self, ab_sequences: list[Optional[SeqRecord]], num_scheme: NumberingScheme,
                            cdr_def: CDRDefinitionScheme) -> (list[Optional[str]], dict):

        if not cdr_def:
            if num_scheme == NumberingScheme.IMGT:
                cdr_def = CDRDefinitionScheme.IMGT
            elif num_scheme == NumberingScheme.KABAT:
                cdr_def = CDRDefinitionScheme.KABAT
            else:
                cdr_def = CDRDefinitionScheme.CHOTHIA

        align_information = align_antibody_sequences(ab_sequences, num_scheme, cdr_def)
        output_sequences = align_information.aligned_sequences
        numbering_json = align_information.to_column_json()

        rows = [sequence_to_genbank_base64_str(s) for s in output_sequences]
        properties = {ANTIBODY_NUMBERING_COLUMN_PROPERTY: numbering_json}

        return (rows, properties)

    def execute(self, request: DataFunctionRequest) -> DataFunctionResponse:

        sequence_column = input_field_to_column(request, 'uiAbSeqCol')
        sequence_column.remove_nulls()
        ab_sequences = column_to_sequences(sequence_column)

        ab_ids = input_field_to_column(request, 'uiIDCol').values

        save_all = boolean_input_field(request, 'uiSaveAll')
        refine_all = boolean_input_field(request, 'uiRefineAll')
        if save_all:
            row_multiplier = 4
        else:
            row_multiplier = 1

        do_numbering = boolean_input_field(request, 'uiDoNumbering')
        num_scheme = string_input_field(request, 'uiNumberingScheme')
        cdr_def = string_input_field(request, 'uiCDRdef')

        # create an ABodyBuilder2 structure predictor with the selected numbering scheme
        predictor = ABodyBuilder2(numbering_scheme = num_scheme)

        # new table columns
        ids = []
        model_number = []
        model_rank = []
        orig_seq = []
        HL_concat_seq = []
        heavy_chain_seq = []
        light_chain_seq = []
        embedded_structures = []

        notifications = []

        for ab_seq, ab_id in zip(ab_sequences, ab_ids):
            ids.extend([ab_id] * row_multiplier)
            orig_seq.extend([str(ab_seq.seq)] * row_multiplier)

            # concatenated sequence is provided for heavy and light chains
            # ABB algorithm identifies individual chains
            sequences = {'H': str(ab_seq.seq).upper(), 'L': str(ab_seq.seq).upper()}
            try:
                antibody = predictor.predict(sequences)
            except Exception as ex:
                notifications.append(Notification(level = NotificationLevel.ERROR,
                                                  title = 'Antibody Structure Prediction',
                                                  summary = f'Error for ID {ab_id}/n{ex.__class__} - {ex}',
                                                  details = f'{ex.__traceback__}'))

                # remove list elements for this broken item
                ids = ids[:-row_multiplier]
                orig_seq = orig_seq[:-row_multiplier]

                continue

            heavy_chain_seq.extend([''.join(residue[1] for residue in antibody.numbered_sequences['H'])] * row_multiplier)
            light_chain_seq.extend([''.join(residue[1] for residue in antibody.numbered_sequences['L'])] * row_multiplier)
            HL_concat_seq.extend([''.join([heavy_chain_seq[-1], light_chain_seq[-1]])] * row_multiplier)

            if save_all:
                model_number.extend([idx + 1 for idx in range(row_multiplier)])
                model_rank.extend([idx + 1 for idx in antibody.ranking])

                for model_idx in range(len(antibody.atoms)):
                    pdb_filename = ''.join([ab_id, f'_Model_{model_idx}_Rank_{antibody.ranking.index(model_idx)}.pdb'])
                    error_caught = False
                    try:
                        antibody.save_single_unrefined(pdb_filename, index = model_idx)
                        if refine_all or antibody.ranking[model_idx] == 0:
                            refine(pdb_filename, pdb_filename, check_for_strained_bonds = True, n_threads = -1)
                        add_errors_as_bfactors(pdb_filename, antibody.error_estimates.mean(0).sqrt().cpu().numpy(),
                                               header = [ABB2_HEADER])

                        # open the file, compress, and encode it
                        with open(pdb_filename, 'r', encoding='utf-8') as pdb_file:
                            pdb_data = pdb_file.read()
                            pdb_zip = gzip.compress(pdb_data.encode())
                            pdb_enc = base64.b64encode(pdb_zip).decode('utf8')
                            embedded_structures.append(pdb_enc)
                    except Exception as ex:
                        notifications.append(Notification(level=NotificationLevel.ERROR,
                                                          title='Antibody Structure Prediction',
                                                          summary=f'Error saving ID {ab_id}, model #{model_idx}/n{ex.__class__} - {ex}',
                                                          details=f'{ex.__traceback__}'))

                        # remove list elements for this broken item
                        ids = ids[:-row_multiplier]
                        orig_seq = orig_seq[:-row_multiplier]
                        heavy_chain_seq = heavy_chain_seq[:-row_multiplier]
                        light_chain_seq = light_chain_seq[:-row_multiplier]
                        HL_concat_seq = HL_concat_seq[:-row_multiplier]
                        model_rank = model_rank[:-row_multiplier]
                        model_number = model_number[:-row_multiplier]

                        error_caught = True
                    finally:
                        if os.path.exists(pdb_filename):
                            os.remove(pdb_filename)
                        if error_caught:
                            continue
            else:
                pdb_filename = '_'.join([ab_id, 'predicted.pdb'])
                error_caught = False
                try:
                    antibody.save(pdb_filename)

                    # open the file, compress, and encode it
                    with open(pdb_filename, 'r', encoding='utf-8') as pdb_file:
                        pdb_data = pdb_file.read()
                        pdb_zip = gzip.compress(pdb_data.encode())
                        pdb_enc = base64.b64encode(pdb_zip).decode('utf8')
                        embedded_structures.append(pdb_enc)
                except Exception as ex:
                    notifications.append(Notification(level=NotificationLevel.ERROR,
                                                      title='Antibody Structure Prediction',
                                                      summary=f'Error saving ID {ab_id}/n{ex.__class__} - {ex}',
                                                      details=f'{ex.__traceback__}'))

                    # remove list elements for this broken item
                    ids = ids[:-row_multiplier]
                    orig_seq = orig_seq[:-row_multiplier]
                    heavy_chain_seq = heavy_chain_seq[:-row_multiplier]
                    light_chain_seq = light_chain_seq[:-row_multiplier]
                    HL_concat_seq = HL_concat_seq[:-row_multiplier]
                    model_rank = model_rank[:-row_multiplier]
                    model_number = model_number[:-row_multiplier]

                    error_caught = True
                finally:
                    if os.path.exists(pdb_filename):
                        os.remove(pdb_filename)
                    if error_caught:
                        continue

        # antibody numbering or not
        # assume not, and set standard values for sequence column
        HL_chain_values = HL_concat_seq
        HL_chain_props = {}
        HL_chain_content_type = 'chemical/x-sequence'
        HL_chain_data_type = DataType.STRING

        if do_numbering:
            try:
                HL_chain_values, HL_chain_props = self._antibody_numbering([string_to_sequence(s, index)
                                                                            for index, s in enumerate(HL_concat_seq)],
                                                                           NumberingScheme.from_str(num_scheme),
                                                                           CDRDefinitionScheme.from_str(cdr_def))
                HL_chain_content_type = 'chemical/x-genbank'
                HL_chain_data_type = DataType.BINARY
            except Exception as ex:
                notifications.append(Notification(level=NotificationLevel.ERROR,
                                                  title='Antibody Structure Prediction - Antibody Numbering',
                                                  summary=f'An unexpected error occurred\n{ex.__class__} - {ex}',
                                                  details=f'{ex.__traceback__}'))

        columns = [ColumnData(name = 'ID', dataType = DataType.STRING,
                              values = ids),
                   ColumnData(name = 'Compressed Structures', dataType = DataType.BINARY,
                              contentType = 'chemical/x-pdb', values = embedded_structures,
                              properties={'Dimension': '3'}),
                   ColumnData(name = 'Concatenated Chains (Heavy + Light)', dataType = HL_chain_data_type,
                              contentType = HL_chain_content_type, values = HL_chain_values, properties = HL_chain_props),
                   ColumnData(name = 'Original Sequence', dataType = DataType.STRING,
                              contentType='chemical/x-sequence', values = orig_seq),
                   ColumnData(name = 'Heavy Chain', dataType = DataType.STRING,
                              contentType = 'chemical/x-sequence', values = heavy_chain_seq),
                   ColumnData(name = 'Light Chain', dataType = DataType.STRING,
                              contentType = 'chemical/x-sequence', values = light_chain_seq)]

        if save_all:
            columns[1:1] = [ColumnData(name = 'Model Number', dataType = DataType.INTEGER,
                                       values = model_number),
                            ColumnData(name = 'Model Rank', dataType = DataType.INTEGER,
                                       values = model_rank)]

        output_table = TableData(tableName = 'Antibody Structure Predictions',
                                 columns = columns)

        return DataFunctionResponse(outputTables = [output_table], notifications = notifications)
