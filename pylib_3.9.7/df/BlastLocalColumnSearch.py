import os

from df.BlastLocalTextSearch import run_blast_search
from df.bio_helper import column_to_sequences
from df.data_transfer import DataFunction, DataFunctionRequest, DataFunctionResponse, string_input_field


class BlastLocalColumnSearch(DataFunction):
    """
    Performs a BLAST search of user selected databases using queries from a sequence column
    """

    def execute(self, request: DataFunctionRequest) -> DataFunctionResponse:
        sequence_column_id = string_input_field(request, 'sequenceColumn')
        id_column_id = string_input_field(request, 'idColumn')
        sequence_column = request.inputColumns[sequence_column_id]
        id_column = None if id_column_id is None else request.inputColumns[id_column_id]
        input_sequences = column_to_sequences(sequence_column, id_column)
        input_sequences = [s for s in input_sequences if s]
        blastdb = string_input_field(request, 'blastDbPath')
        os.environ['BLASTDB'] = blastdb
        return run_blast_search(input_sequences, request, 'Blast column search results', sequence_column.contentType,
                                True)
