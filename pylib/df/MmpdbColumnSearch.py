import glob
import hashlib
import json
import multiprocessing
import os.path
import re
import tempfile
import uuid
from datetime import datetime
from typing import List, Final, Optional

from df.MmpdbDatabaseSearch import mmpdb_query, search_mmpdb

try:
    from mmpdblib import commandline as mmp
except ImportError:
    from mmpdblib import cli as mmp
from pydantic import BaseModel
from rdkit import Chem

from df.chem_helper import column_to_molecules
from df.data_transfer import DataFunction, DataFunctionRequest, DataFunctionResponse, string_input_field, \
    string_list_input_field

VERSION: Final[str] = '0.0.2'


class MmpdbProperties(BaseModel):
    date_created: datetime
    checksum: str
    property_names: List[str]
    database_file: str
    version: str = VERSION


def clean_old_databases(mmpdb_dir: str, database_file: Optional[str] = None) -> None:
    database_files = glob.glob(os.path.join(mmpdb_dir, 'mmpdb*.mmpdb'))
    for f in database_files:
        if not database_file or f != database_file:
            # the internals of MMPDB do not close database connections
            # so if we try to delete the database after accessing it we'll get an error
            # this is an issue for Pyro on Windows.
            try:
                os.remove(f)
            except PermissionError:
                pass


def build_mmp_database(mmpdb_dir: str, request: DataFunctionRequest) -> MmpdbProperties:
    fragments_file = os.path.join(mmpdb_dir, 'mmpdb.fragments.gz')
    property_file = os.path.join(mmpdb_dir, 'mmpdb.props')
    smiles_file = os.path.join(mmpdb_dir, 'mmpdb.smiles')
    settings_file = os.path.join(mmpdb_dir, 'mmpdb.settings')

    structure_column_id = string_input_field(request, 'structureColumn')
    structures = column_to_molecules(request.inputColumns[structure_column_id])
    smiles = [Chem.MolToSmiles(m) if m else None for m in structures]
    property_column_ids = string_list_input_field(request, 'propertyFields')
    properties = [request.inputColumns[column_id].values for column_id in property_column_ids]
    property_names = ['_'.join(request.inputColumns[column_id].name.split()) for column_id in property_column_ids]

    input_data = {'smiles': smiles, 'properties': properties}
    input_data_str = json.dumps(input_data)
    input_data_checksum = hashlib.md5(input_data_str.encode()).hexdigest()

    if not os.path.exists(mmpdb_dir):
        os.makedirs(mmpdb_dir)

    if os.path.exists(settings_file):
        with open(settings_file, 'r') as fh:
            settings_json = fh.read()

        settings = MmpdbProperties.parse_raw(settings_json)
        if settings.checksum == input_data_checksum and settings.version == VERSION:
            clean_old_databases(mmpdb_dir, settings.database_file)
            return settings
        else:
            os.remove(settings_file)

    clean_old_databases(mmpdb_dir)
    uid = str(uuid.uuid4())
    database_file = os.path.join(mmpdb_dir, f'mmpdb_{uid}.mmpdb')

    with open(property_file, 'w') as prop_fh, open(smiles_file, 'w') as smi_fh:
        property_name_str = ' '.join(property_names)
        prop_fh.write(f'ID {property_name_str}\n')
        for compound_num in range(0, len(smiles)):
            smi = smiles[compound_num]
            if not smi:
                continue
            props = [p[compound_num] for p in properties]
            if all(p is None for p in props):
                continue
            prop_str = ' '.join('*' if p is None else str(p) for p in props)
            smi_fh.write(f'{smi} {compound_num}\n')
            prop_fh.write(f'{compound_num} {prop_str}\n')

    # build fragments file
    n_proc = multiprocessing.cpu_count()
    # args = ['fragment', '--in', 'smi', '--out', 'fragments.gz', '--num-jobs', str(n_proc), '--output',
    args = ['fragment', '--in', 'smi', '--num-jobs', str(n_proc), '--output',
            fragments_file, '--delimiter',
            'space', smiles_file]
    try:
        mmp.main(args=args, standalone_mode=False)
    except TypeError:
        mmp.main(args)

    # index fragments
    args = ['index', fragments_file, '-o', database_file]
    try:
        mmp.main(args=args, standalone_mode=False)
    except TypeError:
        mmp.main(args)

    # add properties to database
    args = ['loadprops', '-p', property_file, database_file]
    try:
        mmp.main(args=args, standalone_mode=False)
    except TypeError:
        mmp.main(args)

    os.remove(smiles_file)
    os.remove(property_file)
    os.remove(fragments_file)

    settings = MmpdbProperties(checksum=input_data_checksum, date_created=datetime.now(), property_names=property_names,
                               database_file=database_file)
    with open(settings_file, 'w') as fh:
        fh.write(settings.json())

    return settings


class MmpdbColumnSearch(DataFunction):
    def execute(self, request: DataFunctionRequest) -> DataFunctionResponse:
        structure_column_id = string_input_field(request, 'structureColumn')
        mmpdb_dir = re.sub('[^A-Za-z0-9 ]+', '', structure_column_id)
        mmpdb_dir = re.sub(' ', '_', mmpdb_dir)
        mmpdb_dir = os.path.join(tempfile.gettempdir(), 'ChemChartsData', f'mmpdb_{mmpdb_dir}')

        settings = build_mmp_database(mmpdb_dir, request)

        query_smiles, query_mol = mmpdb_query(request)
        return search_mmpdb(query_smiles, query_mol, settings.property_names, settings.database_file,
                            'MMPDB table search')
