from Bio.Data import CodonTable
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from df.bio_helper import column_to_sequences, sequences_to_column
from df.data_transfer import DataFunction, DataFunctionRequest, DataFunctionResponse, string_input_field
from ruse.bio.bio_util import is_dna_record


def translate(rec: SeqRecord, codon_table_name: str) -> SeqRecord:
    if is_dna_record(rec):
        return rec.translate(codon_table_name)
    # create single naive DNA sequence- maybe should error here instead
    s = rec.seq
    codon_table = CodonTable.unambiguous_dna_by_name[codon_table_name]
    mapping = codon_table.back_table
    codons = [mapping[r] if r in mapping else '-' for r in s]
    s = Seq(''.join(codons))
    rec = SeqRecord(s, rec.id, rec.name)
    return rec


class TranslateSequences(DataFunction):

    def execute(self, request: DataFunctionRequest) -> DataFunctionResponse:
        input_column = next(iter(request.inputColumns.values()))
        input_sequences = column_to_sequences(input_column)
        codon_table_name = string_input_field(request, 'codonTableName', 'Standard')
        output_sequences = [None if s is None else translate(s, codon_table_name) for s in input_sequences]
        output_column = sequences_to_column(output_sequences, f'Translated {input_column.name}', genbank_output=False)
        response = DataFunctionResponse(
            outputColumns=[output_column])
        return response
