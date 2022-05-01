import json
import ruse.bio
from df.bio_helper import column_to_sequences, sequences_to_column

from df.data_transfer import ColumnData, DataFunctionRequest, DataFunctionResponse, DataFunction, string_input_field, \
    DataType
from ruse.bio.phylo_tree import PhyloTree
from ruse.bio.sequence_align import SequenceAlignmentMethod


class MultipleSequenceAlignment(DataFunction):
    def execute(self, request: DataFunctionRequest) -> DataFunctionResponse:
        input_column = next(iter(request.inputColumns.values()))
        input_column.remove_nulls()
        input_sequences = column_to_sequences(input_column)
        alignment_method_str = string_input_field(request, 'alignmentMethod', 'clustalo')

        alignment_method = SequenceAlignmentMethod[alignment_method_str.upper()]

        msa = ruse.bio.sequence_align.MultipleSequenceAlignment()
        msa.align_sequences({}, input_sequences, alignment_method=alignment_method)

        output_sequences = msa.aligned_sequences
        genbank = input_column.contentType == 'chemical/x-genbank'
        if genbank:
            msa.copy_features_to_aligned_sequences()
        output_column = sequences_to_column(output_sequences, 'Aligned Sequence', genbank_output=genbank)

        id_values = [s.id if s else None for s in output_sequences]
        tree = PhyloTree(msa.tree)
        id_column = ColumnData(name='Dendrogram', dataType=DataType.STRING, properties={'tree': json.dumps(tree.data_tree)},
                               values=id_values)

        size = len(input_column.values)
        distances = [None] * size
        msa.add_distances()
        for index, id in zip(range(size), id_values):
            distances[index] = msa.distance_lookup.get(id)
        distance_column = ColumnData(name='Distance', dataType=DataType.DOUBLE, values=distances)
        missing = input_column.missing_null_positions
        output_column.insert_nulls(missing)
        id_column.insert_nulls(missing)
        distance_column.insert_nulls(missing)
        response = DataFunctionResponse(outputColumns=[output_column, id_column, distance_column])
        return response
