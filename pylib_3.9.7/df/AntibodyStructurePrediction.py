import os

from df.bio_helper import column_to_sequences
from df.data_transfer import ColumnData, TableData, DataFunctionRequest, DataFunctionResponse, DataFunction, DataType, \
                             string_input_field, input_field_to_column

from ImmuneBuilder import ABodyBuilder2

class AntibodyStructurePrediction(DataFunction):
    """
    Predicts antibody structure from heavy and light chanin sequences
    """

    def execute(self, request: DataFunctionRequest) -> DataFunctionResponse:

        ab_sequence = column_to_sequences(input_field_to_column(request, 'uiAbSeqCol'))
        ab_id = input_field_to_column(request, 'uiIDCol').values
        num_scheme = string_input_field(request, 'uiNumberingScheme')
        output_dir = string_input_field(request, 'uiOutputDirectory')

        predictor = ABodyBuilder2(numbering_scheme = num_scheme)

        # new table columns
        ids = []
        filenames = []
        orig_seq = []
        HL_concat_seq = []
        heavy_chain_seq = []
        light_chain_seq = []

        for ab_seq, ab_id in zip(ab_sequence, ab_id):
            ids.append(ab_id)
            orig_seq.append(str(ab_seq.seq))

            if ab_seq is None:
                filenames.append(None)
                HL_concat_seq.append(None)
                heavy_chain_seq.append(None)
                light_chain_seq.append(None)
            else:
                # concatenated sequence is provided for heavy and light chains
                # ABB algorithm identifies individual chains
                sequences = {'H': str(ab_seq.seq).upper(), 'L': str(ab_seq.seq).upper()}
                antibody = predictor.predict(sequences)

                filename = os.path.join(output_dir, '_'.join([ab_id, 'predicted.pdb']))
                filenames.append(filename)

                heavy_chain_seq.append(''.join(residue[1] for residue in antibody.numbered_sequences['H']))
                light_chain_seq.append(''.join(residue[1] for residue in antibody.numbered_sequences['L']))
                HL_concat_seq.append(''.join([heavy_chain_seq[-1], light_chain_seq[-1]]))

                antibody.save(filename)

        columns = [ColumnData(name = 'ID', dataType = DataType.STRING,
                              values = ids),
                   ColumnData(name = 'Structure Files', dataType = DataType.STRING,
                              contentType = 'chemical/x-uri', values = filenames,
                              properties = {'Dimension': '3'}),
                   ColumnData(name = 'Concatenated Chains (Heavy + Light)', dataType = DataType.STRING,
                              contentType = 'chemical/x-sequence', values = HL_concat_seq),
                   ColumnData(name = 'Heavy Chain', dataType = DataType.STRING,
                              contentType = 'chemical/x-sequence', values = heavy_chain_seq),
                   ColumnData(name = 'Light Chain', dataType = DataType.STRING,
                              contentType = 'chemical/x-sequence', values = HL_concat_seq),
                   ColumnData(name = 'Original Sequence', dataType = DataType.STRING,
                              contentType='chemical/x-sequence', values = orig_seq)]

        output_table = TableData(tableName = 'Antibody Structure Predictions',
                                 columns = columns)

        return DataFunctionResponse(outputTables = [output_table])