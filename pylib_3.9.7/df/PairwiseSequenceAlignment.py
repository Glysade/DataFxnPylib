from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from df.bio_helper import values_to_sequences, sequences_to_column, query_from_request
from df.data_transfer import DataFunction, DataFunctionRequest, DataFunctionResponse, string_input_field, DataType, \
    ColumnData
from ruse.bio.pairwise_alignment import PairwiseAlignmentMethod, needle_pairwise_alignment, water_pairwise_alignment


class PairwiseSequenceAlignment(DataFunction):

    def execute(self, request: DataFunctionRequest) -> DataFunctionResponse:
        input_column = next(iter(request.inputColumns.values()))
        input_column.remove_nulls()
        genbank = input_column.contentType == 'chemical/x-genbank'
        input_sequences = values_to_sequences(input_column)
        query = query_from_request(request, 'query')

        method_name = string_input_field(request, 'method', 'needle')
        method = PairwiseAlignmentMethod[method_name.upper()]

        if method == PairwiseAlignmentMethod.NEEDLE:
            alignments = needle_pairwise_alignment(query, input_sequences)
        elif method == PairwiseAlignmentMethod.WATER:
            alignments = water_pairwise_alignment(query, input_sequences)
        else:
            raise ValueError()

        index_to_name = {i: s.id for i, s in enumerate(input_sequences)}

        alignment_map = {alignment.target.id: alignment for alignment in alignments}
        alignment_rows = [alignment_map.get(index_to_name[row_num]) for row_num in range(len(input_sequences))]
        aligned_sequences = ['{}|{}'.format(alignment.query.seq, alignment.target.seq) if alignment else None
                             for alignment in alignment_rows]
        scores = [alignment.score if alignment else None for alignment in alignment_rows]
        queries = [alignment.query if alignment else None for alignment in alignment_rows]
        targets = [alignment.target if alignment else None for alignment in alignment_rows]
        for target, query in zip(targets, queries):
            query.id = target.id

        aligned_sequences_column = ColumnData(name='Pairwise Alignment', dataType=DataType.STRING,
                                              contentType='chemical/x-sequence-pair', values=aligned_sequences)
        scores_column = ColumnData(name='Pairwise Alignment Score', dataType=DataType.DOUBLE, values=scores)

        target_column = sequences_to_column(targets, 'Pairwise Alignment Target', genbank_output=genbank)
        query_column = sequences_to_column(queries, 'Pairwise Alignment Query', genbank_output=False)

        aligned_sequences_column.insert_nulls(input_column.missing_null_positions)
        scores_column.insert_nulls(input_column.missing_null_positions)
        target_column.insert_nulls(input_column.missing_null_positions)
        query_column.insert_nulls(input_column.missing_null_positions)

        response = DataFunctionResponse(
            outputColumns=[aligned_sequences_column, scores_column, target_column, query_column])
        return response
