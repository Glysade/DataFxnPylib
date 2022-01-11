import uuid

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from df.bio_helper import values_to_sequences
from df.data_transfer import DataFunction, DataFunctionRequest, DataFunctionResponse, string_input_field, ColumnData, \
    DataType
from ruse.bio.blast_parse import BlastResults
from ruse.bio.blast_search import BlastCreateAndSearch


class BlastTableSearch(DataFunction):
    def execute(self, request: DataFunctionRequest) -> DataFunctionResponse:
        options = {}
        args = {}
        query_input_field = string_input_field(request, 'query')
        if not query_input_field:
            raise ValueError()
        query = SeqRecord(Seq(query_input_field), id=str(uuid.uuid4()), name='', description='')

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

        for hit in results.hits:
            name = hit.data_table_name()
            if name not in name_to_index:
                raise ValueError('Cannot find row for sequence name {}'.format(name))
            row_number = name_to_index[name]
            aligned_sequences[row_number] = f'{hit.query}|{hit.target}'
            e_values[row_number] = hit.evalue
            scores[row_number] = hit.score
            bits[row_number] = hit.bits

        aligned_sequences_column = ColumnData(name='Aligned Sequence', dataType=DataType.STRING,
                                              contentType='chemical/x-sequence-pair', values=aligned_sequences)
        e_values_column = ColumnData(name='Evalue', dataType=DataType.FLOAT, values=e_values)
        score_column = ColumnData(name='Score', dataType=DataType.INTEGER, values=scores)
        bits_column = ColumnData(name='Bits', dataType=DataType.FLOAT, values=bits)

        aligned_sequences_column.insert_nulls(input_column.missing_null_positions)
        e_values_column.insert_nulls(input_column.missing_null_positions)
        score_column.insert_nulls(input_column.missing_null_positions)
        bits_column.insert_nulls(input_column.missing_null_positions)

        response = DataFunctionResponse(
            outputColumns=[aligned_sequences_column, e_values_column, score_column, bits_column])
        return response
