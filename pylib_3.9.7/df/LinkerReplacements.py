import concurrent.futures as cf
from os import cpu_count
from pathlib import Path

from df.chem_helper import column_to_molecules, molecules_to_column
from df.data_transfer import (DataFunction, DataFunctionRequest,
                              DataFunctionResponse, ColumnData, DataType,
                              TableData, boolean_input_field,
                              integer_input_field, string_input_field)
from rdkit import Chem

from df.replace_bioisostere_linkers import replace_linkers

LINKER_DB = Path(__file__).parent.parent.parent / 'Data' / 'chembl_32_bioisostere_linkers.db'
# These 2 values should match those used to make LINKER_DB.  10 and 5
# are the numbers suggested by Ertl et al.
MAX_HEAVIES = 10
MAX_BONDS = 5


def build_output_table(new_mols: list[Chem.Mol], parent_mols: list[Chem.Mol],
                       parent_ids: list[str],
                       parent_id_type: DataType,
                       all_linker_smis: list[list[str]]) -> TableData:
    parent_col = molecules_to_column(parent_mols, "Parent Mols",
                                     DataType.BINARY)
    new_col = molecules_to_column(new_mols, "New Linker Mols",
                                  DataType.BINARY)
    id_col = ColumnData(name='Parent ID', dataType=parent_id_type,
                        values=parent_ids)
    num_linker_cols = 0
    for als in all_linker_smis:
        if len(als) > num_linker_cols:
            num_linker_cols = len(als)
    columns = [parent_col, id_col, new_col]
    for i in range(num_linker_cols):
        linker_smis = []
        for als in all_linker_smis:
            if i < len(als):
                linker_smis.append(als[i])
            else:
                linker_smis.append('')
        columns.append(ColumnData(name=f'Linker {i + 1}',
                                  dataType=DataType.STRING,
                                  values=linker_smis))
    table_data = TableData(tableName='New Linker Mols',
                           columns=columns)
    return table_data


class LinkerReplacements(DataFunction):

    def __init__(self) -> None:
        self._plus_delta_bonds = -1
        self._minus_delta_bonds = -1
        self._match_hbonds = False
        self._structure_column_name = None
        self._parent_mols = None
        self._parent_ids = None
        self._ids_type = None
        self._max_mols_per_input = 1000
        super().__init__()

    def extract_input_data(self, request: DataFunctionRequest) -> None:
        column_id = string_input_field(request, 'structureColumn')
        structure_column = request.inputColumns[column_id]
        self._parent_mols = column_to_molecules(structure_column)
        self._plus_delta_bonds = integer_input_field(request,
                                                     'plusDeltaBonds')
        self._minus_delta_bonds = integer_input_field(request,
                                                      'minusDeltaBonds')
        self._match_hbonds = boolean_input_field(request, 'matchHbonds')
        self._max_mols_per_input = integer_input_field(request,
                                                       'maxMolsPerInput')
        column_id = string_input_field(request, 'idColumn')
        id_column = request.inputColumns[column_id]
        self._parent_ids = id_column.values
        self._ids_type = id_column.dataType

    def do_replacements(self) -> tuple[list[Chem.Mol], list[Chem.Mol], list[str], list[list[str]]]:
        """
        Replace all the linkers in the input set, returning a list of
        new molecules and a corresponding list of the input molecules.
        Note that the input molecules aren't copied, just the reference
        to them, so multiple entries in parent_mols may be pointing to
        the same object.
        """
        all_new_mols = {}
        parent_mols = {}
        all_linker_smis = {}
        # When the dataset is big enough for parallelisation to make a
        # difference, the risk of a combinatorial explosion and hence
        # excessive memory use is very large, so don't go wild on the
        # number of CPUs.
        num_cpus = cpu_count() // 2
        if not num_cpus:
            num_cpus = 1
        with cf.ProcessPoolExecutor(max_workers=num_cpus) as pool:
            futures_to_mol_id = {}
            for mol, mol_id in zip(self._parent_mols, self._parent_ids):
                if mol is None or not mol:
                    continue
                fut = pool.submit(replace_linkers, mol, LINKER_DB,
                                  max_heavies=MAX_HEAVIES,
                                  max_bonds=MAX_BONDS,
                                  plus_length=self._plus_delta_bonds,
                                  minus_length=self._minus_delta_bonds,
                                  match_donors=self._match_hbonds,
                                  match_acceptors=self._match_hbonds,
                                  max_mols_per_input=self._max_mols_per_input)
                futures_to_mol_id[fut] = mol_id
            for fut in cf.as_completed(futures_to_mol_id):
                new_mols, query_cp, linker_smis = fut.result()
                mol_id = futures_to_mol_id[fut]
                if new_mols:
                    all_new_mols[mol_id] = new_mols
                    parent_mols[mol_id] = query_cp
                    all_linker_smis[mol_id] = linker_smis

        # put the output in input order
        new_new_mols = []
        new_parent_mols = []
        new_parent_ids = []
        new_linker_smis = []
        for pid in self._parent_ids:
            try:
                new_new_mols.extend(all_new_mols[pid])
                new_parent_mols.extend([parent_mols[pid]] * len(all_new_mols[pid]))
                new_parent_ids.extend([pid] * len(all_new_mols[pid]))
                new_linker_smis.extend(all_linker_smis[pid])
            except KeyError:
                # For some reason, such as there were no linkers in
                # the molecule, it produced no output.
                pass

        return new_new_mols, new_parent_mols, new_parent_ids, new_linker_smis

    def execute(self, request: DataFunctionRequest) -> DataFunctionResponse:

        self.extract_input_data(request)
        new_mols, parent_mols, parent_ids, new_linker_smis = self.do_replacements()
        output_table = build_output_table(new_mols, parent_mols, parent_ids,
                                          self._ids_type, new_linker_smis)

        return DataFunctionResponse(outputTables=[output_table])
