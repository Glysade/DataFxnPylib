from abc import ABC, abstractmethod
from enum import Enum
from typing import Optional, Any, Dict, List

from pydantic import BaseModel


# I wanted to use this type instead of any for column/inputField values, but it coerces values to strings with used
# with pydantic
# DataValue = Union[str, int, float, bool]


class DataType(Enum):
    BOOLEAN = 'boolean'
    STRING = 'string'
    BINARY = 'binary'
    INTEGER = 'integer'
    FLOAT = 'float'
    STRING_LIST = 'list(string)'
    INTEGER_LIST = 'list(integer)'
    FLOAT_LIST = 'list(float)'


class InputFieldSelectorType(Enum):
    COLUMN = 'column'
    TABLE = 'table'


class InputField(BaseModel):
    id: str
    contentType: Optional[str]
    dataType: DataType
    data: Any
    selectorType: Optional[InputFieldSelectorType]


class ColumnData(BaseModel):
    name: str
    dataType: DataType
    properties: Dict[str, str] = {}
    contentType: Optional[str]
    values: List[Any]
    missing_null_positions: List[int] = []

    # TODO: maybe just use column data as transfer object move this functionality to another class
    def remove_nulls(self) -> None:
        positions = [i for i, v in enumerate(self.values) if v is None]
        positions.reverse()
        for p in positions:
            self.values.pop(p)
        positions.reverse()
        self.missing_null_positions = positions

    def insert_nulls(self, positions: Optional[List[int]] = None, offset: int = 0) -> None:
        if not positions:
            positions = self.missing_null_positions
        for p in positions:
            self.values.insert(p + offset, None)
        self.missing_null_positions = []


class TableData(BaseModel):
    tableName: str
    columns: List[ColumnData]


class DataFunctionRequest(BaseModel):
    inputColumns: Dict[str, ColumnData] = []
    inputFields: Dict[str, InputField] = {}
    serviceName: str
    id: str
    script: Optional[str]


class DataFunctionResponse(BaseModel):
    outputColumns: List[ColumnData] = []
    outputTables: List[TableData] = []


class DataFunction(ABC):
    @abstractmethod
    def execute(self, request: DataFunctionRequest) -> DataFunctionResponse:
        pass


def _input_field_value(request: DataFunctionRequest, field: str) -> Optional[Any]:
    input_fields = request.inputFields
    if field not in input_fields:
        return None
    return input_fields[field].data


def _check_data_type(request: DataFunctionRequest, field: str, expected_type: DataType) -> None:
    data_type = request.inputFields[field].dataType
    if data_type != expected_type:
        raise ValueError(f'Input field {field} has data type {data_type} expected {expected_type}')


def string_input_field(request: DataFunctionRequest, field: str, default_value: Optional[str] = None) -> Optional[str]:
    data = _input_field_value(request, field)
    if not data:
        return default_value
    _check_data_type(request, field, DataType.STRING)
    return data


def integer_input_field(request: DataFunctionRequest, field: str, default_value: Optional[int] = None) -> Optional[int]:
    data = _input_field_value(request, field)
    if not data:
        return default_value
    _check_data_type(request, field, DataType.INTEGER)
    return data


def boolean_input_field(request: DataFunctionRequest, field: str, default_value: Optional[bool] = None) -> Optional[
    int]:
    data = _input_field_value(request, field)
    if not data:
        return default_value
    _check_data_type(request, field, DataType.BOOLEAN)
    return data


def float_input_field(request: DataFunctionRequest, field: str, default_value: Optional[float] = None) -> Optional[
    float]:
    data = _input_field_value(request, field)
    if not data:
        return default_value
    _check_data_type(request, field, DataType.FLOAT)
    return data


def binary_input_field(request: DataFunctionRequest, field: str, default_value: Optional[str] = None) -> Optional[str]:
    data = _input_field_value(request, field)
    if not data:
        return default_value
    _check_data_type(request, field, DataType.BINARY)
    return data


def string_list_input_field(request: DataFunctionRequest, field: str, default_value: Optional[List[str]] = None) -> \
        Optional[List[str]]:
    data = _input_field_value(request, field)
    if not data:
        return default_value
    _check_data_type(request, field, DataType.STRING_LIST)
    return data
