import base64
import gzip
from typing import Optional, List

from rdkit.Chem.rdchem import Mol

from df.data_transfer import ColumnData, DataType
from ruse.rdkit.rdkit_utils import type_to_format, string_to_mol, mol_to_string


def column_to_molecules(column: ColumnData) -> List[Optional[Mol]]:
    mols = [None if v is None else value_to_molecule(v, column.dataType, column.contentType) for v in column.values]
    return mols


def _decode_binary(data: str) -> str:
    binary_data = base64.b64decode(data)
    value = gzip.decompress(binary_data)
    for charset in ['ascii', 'utf-8', 'latin-1']:
        try:
            value = value.decode(charset)
            return value
        except UnicodeDecodeError:
            pass
    raise UnicodeDecodeError('Unable to decode binary data')


def _default_content_type(data_type: DataType) -> str:
    if data_type == DataType.STRING:
        content_type = 'chemical/x-smiles'
    elif data_type == DataType.BINARY:
        content_type = 'chemical/x-mdl-molfile'
    return content_type


def value_to_molecule(v: str, data_type: DataType, content_type: Optional[str] = None) -> Optional[Mol]:
    if data_type != DataType.BINARY and data_type != data_type.STRING:
        raise ValueError(f'Can\'t convert DataType {data_type} to molecule')
    if content_type is None:
        content_type = _default_content_type(data_type)

    data = data_type == _decode_binary(v) if data_type == DataType.BINARY else v
    fmt = type_to_format(content_type)
    return string_to_mol(fmt, data)


def _encode_binary(data: str) -> str:
    data_bytes = bytes(data, 'utf-8')
    zipped_data = gzip.compress(data_bytes)
    encoded_data = base64.b64encode(zipped_data).decode('utf-8')
    return encoded_data


def molecules_to_column(mols: List[Optional[Mol]], column_name: str, data_type: DataType,
                        content_type: Optional[str] = None) -> ColumnData:
    if content_type is None:
        content_type = _default_content_type(data_type)
    values = [None if m is None else molecule_to_value(m, data_type, content_type) for m in mols]
    return ColumnData(name=column_name, dataType=data_type, content_type=content_type, values=values)


def molecule_to_value(mol: Mol, data_type: DataType, content_type: Optional[str]) -> str:
    if data_type != DataType.BINARY and data_type != data_type.STRING:
        raise ValueError(f'Can\'t convert molecule to DataType {data_type}')
    if content_type is None:
        content_type = _default_content_type(data_type)
    fmt = type_to_format(content_type)
    value = mol_to_string(fmt, mol)
    data = data_type == _decode_binary(value) if data_type == DataType.BINARY else value
    return data
