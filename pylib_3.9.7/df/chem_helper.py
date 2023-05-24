"""
==============
chem_helper.py
==============

Helper functions for data functions that manipulate chemical structures.

The `RDKit <https://rdkit.org>`_ is used to handle chemistry

For those data functions that encode or decode molecules using a :class:`~df.data_transfer.DataType` argument
with an optional content type argument, if the content_type is not set, chemical/x-smiles will be used if the
data_type is :attr:`~df.data_transfer.DataType.STRING`, otherwise chemical/x-mdl-molfile if the data_type is
:attr:`~df.data_transfer.DataType.BINARY`.

"""

import base64
import gzip
from typing import Optional, List

from rdkit.Chem.rdchem import Mol

from df.data_transfer import ColumnData, DataType, DataFunctionRequest
from ruse.rdkit.rdkit_utils import type_to_format, string_to_mol, mol_to_string


def column_to_molecules(column: ColumnData,
                        do_sanitize_mols: Optional[bool]=True) -> List[Optional[Mol]]:
    """
    Converts a structure column to a list of molecules

    :param column: the input column
    :param do_sanitize_mols: whether to run sanitize_mol on the final mols
    :return: a list of molecules
    """

    mols = [None if v is None else value_to_molecule(v, column.dataType, column.contentType, do_sanitize_mols) for v in column.values]
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
    raise ValueError('Unable to decode binary data')


def _default_content_type(data_type: DataType) -> str:
    content_type = None
    if data_type == DataType.STRING:
        content_type = 'chemical/x-smiles'
    elif data_type == DataType.BINARY:
        content_type = 'chemical/x-mdl-molfile'
    if not content_type:
        raise ValueError(f'Unable to determine contentType for dataType {data_type}')
    return content_type


def value_to_molecule(v: str, data_type: DataType, content_type: Optional[str] = None,
                        do_sanitize_mols: Optional[bool]=True) -> Optional[Mol]:
    """
    Converts a string value to a molecule.

    :param v: the string value
    :param data_type: the data type for the string should be :attr:`~df.data_transfer.DataType.STRING` or :attr:`~df.data_transfer.DataType.BINARY`
    :param content_type: molecular content type
    :param do_sanitize_mols: whether to run sanitize_mol on the final mols
    :return: a molecule
    """
    if data_type != DataType.BINARY and data_type != data_type.STRING:
        raise ValueError(f'Can\'t convert DataType {data_type} to molecule')
    if content_type is None:
        content_type = _default_content_type(data_type)

    data = _decode_binary(v) if data_type == DataType.BINARY else v
    fmt = type_to_format(content_type)
    return string_to_mol(fmt, data, do_sanitize_mols)


def input_field_to_molecule(request: DataFunctionRequest, query_data_field: str = 'queryData') -> Optional[Mol]:
    """
    Converts an input field request to a molecule.

    :param request: the request
    :param query_data_field: the input field name
    :return: a molecule
    """
    input_field = request.inputFields[query_data_field]
    assert input_field.dataType == DataType.STRING
    mol_string = str(input_field.data)
    return value_to_molecule(mol_string, DataType.STRING, input_field.contentType)


def _encode_binary(data: str) -> str:
    zipped_data = gzip.compress(data.encode('utf-8'))
    encoded_data = base64.b64encode(zipped_data).decode('utf-8')
    return encoded_data


def molecules_to_column(mols: List[Optional[Mol]], column_name: str, data_type: DataType,
                        content_type: Optional[str] = None) -> ColumnData:
    """
    Converts a list of molecules to a Spotfire column.

    :param mols: the molecules
    :param column_name: the name for the column
    :param data_type: the data type for the column should be :attr:`~df.data_transfer.DataType.STRING` or :attr:`~df.data_transfer.DataType.BINARY`
    :param content_type: an optional content type for encoding the molecules
    :return: a new column
    """
    if content_type is None:
        content_type = _default_content_type(data_type)
    values = [None if m is None else molecule_to_value(m, data_type, content_type) for m in mols]
    return ColumnData(name=column_name, dataType=data_type, contentType=content_type, values=values)


def molecule_to_value(mol: Mol, data_type: DataType, content_type: Optional[str]) -> str:
    """
    Converts a molecule to an encoded string

    :param mol: the molecule
    :param data_type: the data type for the output string should be :attr:`~df.data_transfer.DataType.STRING` or :attr:`~df.data_transfer.DataType.BINARY`
    :param content_type: an optional content type for encoding the molecule
    :return: encoded molecule
    """
    if data_type != DataType.BINARY and data_type != data_type.STRING:
        raise ValueError(f'Can\'t convert molecule to DataType {data_type}')
    if content_type is None:
        content_type = _default_content_type(data_type)
    fmt = type_to_format(content_type)
    value = mol_to_string(fmt, mol)
    data = _encode_binary(value) if data_type == DataType.BINARY else value
    return data
