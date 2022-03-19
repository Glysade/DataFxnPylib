import os

from df.BlastLocalTextSearch import run_blast_search
from df.bio_helper import values_to_sequences
from df.data_transfer import DataFunction, DataFunctionRequest, DataFunctionResponse, string_input_field


class BlastLocalColumnSearch(DataFunction):

    def execute(self, request: DataFunctionRequest) -> DataFunctionResponse:
        input_column = next(iter(request.inputColumns.values()))
        input_column.remove_nulls()
        input_sequences = values_to_sequences(input_column)
        blastdb = string_input_field(request, 'blastDbPath')
        os.environ['BLASTDB'] = blastdb
        return run_blast_search(input_sequences, request, 'Blast column search results', input_column.contentType)
