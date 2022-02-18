from df.bio_helper import values_to_sequences, query_from_request, sequences_to_column
from df.data_transfer import DataFunction, DataFunctionRequest, DataFunctionResponse, ColumnData, \
    DataType, boolean_input_field
from ruse.bio.blast_parse import BlastResults, build_common_alignments
from ruse.bio.blast_search import BlastCreateAndSearch


class BlastTableSearch(DataFunction):
    def execute(self, request: DataFunctionRequest) -> DataFunctionResponse:
        options = {}
        args = {}

        query = query_from_request(request)
        show_multiple_alignments = boolean_input_field(request, 'showMultipleAlignments', False)
        # if method in ['TBLASTX', 'BLASTX']:
        #     show_multiple_alignments = False

        input_column = next(iter(request.inputColumns.values()))
        input_column.remove_nulls()
        input_sequences = values_to_sequences(input_column)

        blast = BlastCreateAndSearch()
        blast.search_blast_sequences(query, input_sequences, **args, options=options)
        if blast.search.error is not None:
            raise ValueError(f'Blast search exception: {blast.search.error}')

        results = BlastResults()
        results.parse(blast.search.output_file())

        name_to_index = {s.id: i for i, s in enumerate(input_sequences)}
        size = len(input_column.values)
        aligned_sequences = [None] * size
        e_values = [None] * size
        scores = [None] * size
        bits = [None] * size
        hits_in_order = [None] * size

        for hit in results.hits:
            name = hit.data_table_name()
            if name not in name_to_index:
                raise ValueError('Cannot find row for sequence name {}'.format(name))
            row_number = name_to_index[name]  # int
            aligned_sequences[row_number] = f'{hit.query}|{hit.target}'
            e_values[row_number] = hit.evalue
            scores[row_number] = hit.score
            bits[row_number] = hit.bits
            hits_in_order[row_number] = hit

        aligned_sequences_column = ColumnData(name='Sequence Pairs', dataType=DataType.STRING,
                                              contentType='chemical/x-sequence-pair', values=aligned_sequences)
        e_values_column = ColumnData(name='Evalue', dataType=DataType.FLOAT, values=e_values)
        score_column = ColumnData(name='Score', dataType=DataType.INTEGER, values=scores)
        bits_column = ColumnData(name='Bits', dataType=DataType.FLOAT, values=bits)

        null_positions = input_column.missing_null_positions
        aligned_sequences_column.insert_nulls(null_positions)
        e_values_column.insert_nulls(null_positions)
        score_column.insert_nulls(null_positions)
        bits_column.insert_nulls(null_positions)
        columns = [aligned_sequences_column, e_values_column, score_column, bits_column]
        if show_multiple_alignments:
            multiple_alignments = build_common_alignments(query, hits_in_order)
            multiple_alignment_column = sequences_to_column(multiple_alignments[1:], 'Aligned Sequence', True)
            multiple_alignment_column.insert_nulls(null_positions)
            columns.insert(0, multiple_alignment_column)
        response = DataFunctionResponse(outputColumns=columns)
        return response
