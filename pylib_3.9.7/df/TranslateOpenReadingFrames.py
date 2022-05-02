from typing import Tuple, List

from Bio.Data import CodonTable
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from df.bio_helper import column_to_sequences, sequences_to_column
from df.data_transfer import DataFunction, DataFunctionRequest, DataFunctionResponse, ColumnData, DataType, TableData, \
    integer_input_field, string_input_field


# adapted from BioPython tutorial
def find_orfs_with_trans(seq: Seq, min_protein_length: int = 50,
                         trans_table: str = 'Standard') -> List[Tuple[int, int, int, Seq, int]]:
    answer = []
    seq_len = len(seq)
    for strand, nuc in [(+1, seq), (-1, seq.reverse_complement())]:
        for frame in range(3):
            trans = nuc[frame:].translate(trans_table)
            trans_len = len(trans)
            aa_start = 0
            while aa_start < trans_len:
                aa_end = trans.find("*", aa_start)
                if aa_end == -1:
                    aa_end = trans_len
                if aa_end - aa_start >= min_protein_length:
                    if strand == 1:
                        start = frame + aa_start * 3
                        end = frame + aa_end * 3 + 3
                        if end >= seq_len:
                            end = seq_len - 3
                    else:
                        start = seq_len - frame - aa_end * 3 - 3
                        if start < 0:
                            start += 3
                        end = seq_len - frame - aa_start * 3
                    answer.append((start, end, strand, trans[aa_start:aa_end], frame))
                aa_start = aa_end + 1
    answer.sort()
    return answer


class TranslateOpenReadingFrames(DataFunction):
    """
    Translates DNA from an input column using all open reading frames, and outputs the protein found into a new table
    """

    def execute(self, request: DataFunctionRequest) -> DataFunctionResponse:
        sequence_column_id = string_input_field(request, 'sequenceColumn')
        id_column_id = string_input_field(request, 'idColumn')
        sequence_column = request.inputColumns[sequence_column_id]
        id_column = None if id_column_id is None else request.inputColumns[id_column_id]
        input_sequences = column_to_sequences(sequence_column, id_column)
        input_sequences = [s for s in input_sequences if s]
        min_protein_length = integer_input_field(request, 'minimumProteinLength', 100)
        codon_table_name = string_input_field(request, 'codonTableName', 'Standard')

        sequence_number = []
        sequence_id = []
        sequence = []
        length = []
        strands = []
        starts = []
        ends = []
        frames = []

        for index, seq in enumerate(input_sequences):
            proteins = find_orfs_with_trans(seq.seq, min_protein_length, codon_table_name)
            for start, end, strand, pro, frame in proteins:
                sequence_number.append(index + 1)
                sequence_id.append(seq.id)
                sequence.append(SeqRecord(pro))
                strands.append('+' if strand == 1 else '-')
                length.append(len(pro))
                starts.append(start + 1)
                ends.append(end + 1)
                frames.append(frame + 1)

        protein_column = sequences_to_column(sequence, 'Protein', genbank_output=False)
        output_columns = [
            ColumnData(name='Sequence Number', dataType=DataType.LONG, values=sequence_number),
            ColumnData(name='Sequence Id', dataType=DataType.STRING, values=sequence_id),
            protein_column,
            ColumnData(name='Length', dataType=DataType.LONG, values=length),
            ColumnData(name='Strand', dataType=DataType.STRING, values=strands),
            ColumnData(name='Frame', dataType=DataType.LONG, values=frames),
            ColumnData(name='Start', dataType=DataType.LONG, values=starts),
            ColumnData(name='End', dataType=DataType.LONG, values=ends)
        ]
        output_table = TableData(tableName='Open Reading Frames', columns=output_columns)
        return DataFunctionResponse(outputTables=[output_table])
