import os

from df.BlastLocalTextSearch import run_blast_search
from df.bio_helper import column_to_sequences
from df.data_transfer import DataFunction, DataFunctionRequest, DataFunctionResponse, string_input_field, \
                             input_field_to_column


class BlastLocalColumnSearch(DataFunction):
    """
    Performs a BLAST search of user selected databases using queries from a sequence column
    """

    def execute(self, request: DataFunctionRequest) -> DataFunctionResponse:
        # when marking/filtering is enabled, the name column ID
        # is suffixed with _both, _filteredRows, _markedRow_<Marking ID>, or _both_<Marking ID>
        # use helper function to access column to avoid ambiguity
        sequence_column = input_field_to_column(request, 'sequenceColumn')
        id_column = input_field_to_column(request, 'idColumn')

        input_sequences = column_to_sequences(sequence_column, id_column)
        input_sequences = [s for s in input_sequences if s]

        blastdb = string_input_field(request, 'blastDbPath')
        os.environ['BLASTDB'] = blastdb

        response = run_blast_search(input_sequences, request, 'Blast column search results',
                                sequence_column.contentType,
                                True)

        return response
