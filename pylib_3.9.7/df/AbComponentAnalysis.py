import json
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import substitution_matrices

from df.bio_helper import column_to_sequences, sequences_to_column
from df.data_transfer import DataFunction, DataFunctionRequest, DataFunctionResponse, string_input_field, DataType, \
    ColumnData, TableData


class AbComponentAnalysis(DataFunction):
    """
    A data function that creates an analysis of property changes in one part of an antibody while other parts are held
    constant
    """

    def score_match(self, pair, matrix, gap_mismatch):
        if pair in matrix:
            return matrix[pair]
        elif (tuple(reversed(pair))) in matrix:
            return matrix[(tuple(reversed(pair)))]
        elif pair == ('-', '-'):
            return 1
        else:
            return gap_mismatch

    def score_pairwise(self, seq1, seq2, matrix, gap_mismatch):
        score = 0
        gap = False
        seq1 = seq1.upper()
        seq2 = seq2.upper()

        # make sequences same length
        if len(seq1) != len(seq2):
            if len(seq1) > len(seq2):
                seq2 += "-" * (len(seq1) - len(seq2))
            else:
                seq1 += "-" * (len(seq2) - len(seq1))

        for i in range(len(seq1)):
            pair = (seq1[i], seq2[i])
            if not gap:
                score += self.score_match(pair, matrix, gap_mismatch)

        return score / len(seq1)

    def execute(self, request: DataFunctionRequest) -> DataFunctionResponse:
        inColumns = request.inputColumns
        seqColumnIds = request.inputFields['sequenceColumns'].data
        dataColumnIds = request.inputFields['dataColumns'].data
        gapMismatch = request.inputFields['gapMismatch'].data
        idColumn = request.inputFields['idColumn'].data
        blosum = substitution_matrices.load("BLOSUM62")

        src = []
        lhsId = []
        rhsId = []
        lhsSeq = []
        rhsSeq = []
        seqPair = []
        blosumScore = []
        dataColumnVals = {dataColumnId: [] for dataColumnId in dataColumnIds}
        print(f'dfvs = {dataColumnVals}')

        uniqueSeq = {}
        for seqColumnId in seqColumnIds:
            unique = {}
            uniqueSeq[seqColumnId] = unique
            values = inColumns[seqColumnId].values
            for idx, val in enumerate(values):
                if val:
                    if val in unique:
                        unique[val].append(idx)
                    else:
                        unique[val] = [idx]

        for seqColumnId in seqColumnIds:
            # create list of sequences per row, excluding current column
            values = []
            for id in seqColumnIds:
                if id != seqColumnId:
                    values.append(inColumns[id].values)
            rowsSeq = [list(x) for x in zip(*values)]

            uniqueColumnSeq = list(uniqueSeq[seqColumnId].keys())
            for idx1, seq1 in enumerate(uniqueColumnSeq):
                for idx2, seq2 in enumerate(uniqueColumnSeq):
                    if idx1 < idx2:
                        seq1Rows = uniqueSeq[seqColumnId][seq1]
                        seq2Rows = uniqueSeq[seqColumnId][seq2]
                        for ridx1 in seq1Rows:
                            for ridx2 in seq2Rows:
                                if rowsSeq[ridx1] == rowsSeq[ridx2]:
                                    src.append(inColumns[seqColumnId].name)
                                    lhsId.append(inColumns[idColumn].values[ridx1])
                                    rhsId.append(inColumns[idColumn].values[ridx2])
                                    lhsSeq.append(seq1.upper())
                                    rhsSeq.append(seq2.upper())
                                    seqPair.append(f'{seq1.upper()}|{seq2.upper()}')
                                    blosumScore.append(self.score_pairwise(seq1, seq2, blosum, gapMismatch))

                                    for dcId in dataColumnIds:
                                        data1 = inColumns[dcId].values[ridx1]
                                        data2 = inColumns[dcId].values[ridx2]
                                        if data1 is not None and data2 is not None:
                                            dataColumnVals[dcId].append(
                                                data1 - data2)
                                        else:
                                            dataColumnVals[dcId].append(None)

        columns = []
        columns.append(ColumnData(name='Source Column', dataType=DataType.STRING, values=src))
        columns.append(ColumnData(name='LHS ID', dataType=inColumns[idColumn].dataType, values=lhsId))
        columns.append(ColumnData(name='RHS ID', dataType=inColumns[idColumn].dataType, values=rhsId))
        columns.append(
            ColumnData(name='LHS Sequence', dataType=DataType.STRING, contentType='chemical/x-sequence', values=lhsSeq))
        columns.append(
            ColumnData(name='RHS Sequence', dataType=DataType.STRING, contentType='chemical/x-sequence', values=rhsSeq))
        columns.append(
            ColumnData(name='Sequence Pair', dataType=DataType.STRING, contentType='chemical/x-sequence-pair',
                       values=seqPair))
        columns.append(ColumnData(name='BLOSUM Score', dataType=DataType.DOUBLE, values=blosumScore))
        for id in dataColumnVals:
            columns.append(
                ColumnData(name=inColumns[id].name, dataType=inColumns[id].dataType, values=dataColumnVals[id]))

        response = DataFunctionResponse(outputTables=[TableData(tableName='Ab Analysis Result', columns=columns)])
        return response
