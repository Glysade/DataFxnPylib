import json

from df.bio_helper import column_to_sequences
from df.data_transfer import ColumnData, DataFunctionRequest, DataFunctionResponse, DataFunction, string_input_field, \
    DataType, input_field_to_column
from ruse.bio.antibody import align_antibody_sequences, ANTIBODY_NUMBERING_COLUMN_PROPERTY
from ruse.bio.bio_data_table_helper import sequence_to_genbank_base64_str

from ImmuneBuilder import ABodyBuilder2

class AntibodyStructurePrediction(DataFunction):
    """
    Predicts antibody structure from heavy and light chanin sequences
    """

    def execute(self, request: DataFunctionRequest) -> DataFunctionResponse:

        predictor = ABodyBuilder2()

        output_file = "c:\Temp\my_antibody.pdb"
        sequences = {
            'H': 'EVQLVESGGGVVQPGGSLRLSCAASGFTFNSYGMHWVRQAPGKGLEWVAFIRYDGGNKYYADSVKGRFTISRDNSKNTLYLQMKSLRAEDTAVYYCANLKDSRYSGSYYDYWGQGTLVTVS',
            'L': 'VIWMTQSPSSLSASVGDRVTITCQASQDIRFYLNWYQQKPGKAPKLLISDASNMETGVPSRFSGSGSGTDFTFTISSLQPEDIATYYCQQYDNLPFTFGPGTKVDFK'}

        antibody = predictor.predict(sequences)
        antibody.save(output_file)

        return DataFunctionResponse(outputColumns=[])

        # input_column = input_field_to_column(request, 'sequenceColumn')
        # input_column.remove_nulls()
        # input_sequences = column_to_sequences(input_column)
        # numbering_scheme = string_input_field(request, 'numberingScheme', 'chothia')
        # cdr_definition = string_input_field(request, 'cdrDefinition', 'chothia')
        # if not cdr_definition:
        #     if numbering_scheme == 'imgt':
        #         cdr_definition = 'imgt'
        #     elif numbering_scheme == 'kabat':
        #         cdr_definition = 'kabat'
        #     else:
        #         cdr_definition = 'chothia'

        # align_information = align_antibody_sequences(input_sequences, numbering_scheme, cdr_definition)
        # output_sequences = align_information.aligned_sequences
        # numbering_json = align_information.to_column_json()

        # rows = [sequence_to_genbank_base64_str(s) for s in output_sequences]
        # properties = {ANTIBODY_NUMBERING_COLUMN_PROPERTY: numbering_json}

        # output_column = ColumnData(name='Aligned Sequence', dataType=DataType.BINARY, properties=properties,
        #                            contentType='chemical/x-genbank', values=rows)
        # output_column.insert_nulls(input_column.missing_null_positions)
        # response = DataFunctionResponse(outputColumns=[output_column])
        # return response
