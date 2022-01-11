from df.BlastLocalTextSearch import run_blast_search
from df.bio_helper import values_to_sequences
from df.data_transfer import DataFunction, DataFunctionRequest, DataFunctionResponse


class BlastLocalColumnSearch(DataFunction):

    def execute(self, request: DataFunctionRequest) -> DataFunctionResponse:
        input_column = next(iter(request.inputColumns.values()))
        input_column.remove_nulls()
        input_sequences = values_to_sequences(input_column)
        return run_blast_search(input_sequences, request, input_column.contentType)
