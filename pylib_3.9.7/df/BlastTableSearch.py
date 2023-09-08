from df.bio_helper import column_to_sequences, query_from_request, sequences_to_column
from df.data_transfer import DataFunction, DataFunctionRequest, DataFunctionResponse, ColumnData, \
    DataType, boolean_input_field, Notification, NotificationLevel
from ruse.bio.blast_parse import BlastResults, build_common_alignments
from ruse.bio.blast_search import BlastCreateAndSearch


class BlastTableSearch(DataFunction):
    """
    Performs a blast search of the sequences in the provided table using the provided query sequence
    """

    def execute(self, request: DataFunctionRequest) -> DataFunctionResponse:
        options = {}
        args = {}
        notifications = []  # store notification for user in this list using .append(Notification(...))

        query = query_from_request(request)
        show_multiple_alignments = boolean_input_field(request, 'showMultipleAlignments', False)
        # if method in ['TBLASTX', 'BLASTX']:
        #     show_multiple_alignments = False

        input_column = next(iter(request.inputColumns.values()))
        input_column.remove_nulls()
        input_sequences = column_to_sequences(input_column)

        if input_sequences:
            blast = BlastCreateAndSearch()
            blast.search_blast_sequences(query, input_sequences, **args, options=options)
            if blast.search.error is not None:
                raise ValueError(f'Blast search exception: {blast.search.error}')

            results = BlastResults()
            results.parse(blast.search.output_file())
            hits = results.hits
        else:
            hits = []

        name_to_index = {s.id: i for i, s in enumerate(input_sequences)}
        size = len(input_column.values)
        aligned_sequences = [None] * size
        e_values = [None] * size
        scores = [None] * size
        bits = [None] * size
        hits_in_order = [None] * size

        for hit in hits:
            name = hit.data_table_name()
            if name not in name_to_index:
                raise ValueError('Cannot find row for sequence name {}'.format(name))
            row_number: int = name_to_index[name]
            aligned_sequences[row_number] = f'{hit.query}|{hit.target}'
            e_values[row_number] = hit.evalue
            scores[row_number] = hit.score
            bits[row_number] = hit.bits
            hits_in_order[row_number] = hit

        aligned_sequences_column = ColumnData(name='Sequence Pairs', dataType=DataType.STRING,
                                              contentType='chemical/x-sequence-pair', values=aligned_sequences)
        e_values_column = ColumnData(name='Evalue', dataType=DataType.DOUBLE, values=e_values)
        score_column = ColumnData(name='Score', dataType=DataType.LONG, values=scores)
        bits_column = ColumnData(name='Bits', dataType=DataType.DOUBLE, values=bits)

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
            columns.append(multiple_alignment_column)

        if e_values.count(None) == len(e_values) and \
           scores.count(None) == len(scores) and \
           bits.count(None) == len(bits):
            notifications.append(Notification(level = NotificationLevel.INFORMATION,
                                              title = 'BLAST Table Search',
                                              summary = 'No hits identified in search. Data columns will be empty.',
                                              details = 'The search did not return results for any query sequence. ' +
                                                        'The resulting data columns will not contain any values. ' +
                                                        'This is expected behavior when no hits are found for any query.'))

        response = DataFunctionResponse(outputColumns = columns,
                                        notifications = notifications)

        return response
