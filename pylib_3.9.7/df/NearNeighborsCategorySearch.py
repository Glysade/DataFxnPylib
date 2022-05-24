import json
from typing import List, Any, Optional

from df.data_transfer import DataFunction, DataFunctionRequest, DataFunctionResponse, string_input_field, \
    double_input_field, DataType, ColumnData


class NearNeighborsCategorySearch(DataFunction):
    """
    Perform a remote Blast sequence search against an Entrez database
    """

    def execute(self, request: DataFunctionRequest) -> DataFunctionResponse:
        nn_column_id = string_input_field(request, 'nearNeighborsColumn')
        nn_values = request.inputColumns[nn_column_id].values
        nn_lists: List[Any] = [json.loads(v) if v else None for v in nn_values]
        categories: List[str] = request.inputColumns[string_input_field(request, 'categoryColumn')].values
        experimental: List[Optional[float]] = request.inputColumns[
            string_input_field(request, 'experimentalActivityColumn')].values
        predicted: List[float] = request.inputColumns[string_input_field(request, 'predictedActivityColumn')].values
        min_activity: Optional[float] = double_input_field(request, 'minActivity')
        max_activity: Optional[float] = double_input_field(request, 'maxActivity')

        virtual_neighbors: List[bool] = []
        for idx, nn in enumerate(nn_lists):
            virtual_neighbors.append(False)
            category = categories[idx]
            if category.lower() == 'virtual':
                continue
            exp = experimental[idx]
            if min_activity is not None and exp < min_activity:
                continue
            if max_activity is not None and exp > max_activity:
                continue
            found_virtual = False
            for nbr in nn['nbrs']:
                if found_virtual:
                    break
                other_idx = nbr['idx']
                other_category = categories[other_idx]
                if other_category.lower() == 'real':
                    continue
                pre = predicted[other_idx]
                if min_activity is not None and pre < min_activity:
                    continue
                if max_activity is not None and pre > max_activity:
                    continue
                found_virtual = True
            virtual_neighbors[idx] = found_virtual

        virtual_neighbors_column = ColumnData(name='Virtual neighbors', dataType=DataType.BOOLEAN,
                                              values=virtual_neighbors)
        response = DataFunctionResponse(outputColumns=[virtual_neighbors_column])
        return response
