from Bio.Data.CodonTable import standard_dna_table
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from df.bio_helper import values_to_sequences, sequences_to_column
from df.data_transfer import DataFunction, DataFunctionRequest, DataFunctionResponse
from ruse.bio.bio_util import is_dna_record


def translate(rec: SeqRecord) -> SeqRecord:
    if is_dna_record(rec):
        return rec.translate()
    # create single naive DNA sequence- maybe should error here instead
    s = rec.seq
    mapping = standard_dna_table.back_table
    codons = [mapping[r] if r in mapping else '-' for r in s]
    s = Seq(''.join(codons))
    rec = SeqRecord(s, rec.id, rec.name)
    return rec


class TranslateSequences(DataFunction):

    def execute(self, request: DataFunctionRequest) -> DataFunctionResponse:
        input_column = next(iter(request.inputColumns.values()))
        input_sequences = values_to_sequences(input_column)
        output_sequences = [None if s is None else translate(s) for s in input_sequences]
        output_column = sequences_to_column(output_sequences, f'Translated {input_column.name}', genbank_output=False)
        response = DataFunctionResponse(
            outputColumns=[output_column])
        return response
