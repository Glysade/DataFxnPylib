from pathlib import Path

from df.chem_helper import column_to_molecules, molecules_to_column
from df.data_transfer import (DataFunction, DataFunctionRequest,
                              DataFunctionResponse, DataType, TableData,
                              boolean_input_field,
                              integer_input_field, string_input_field)
from rdkit import Chem

from df.replace_bioisostere_linkers import replace_linkers

LINKER_DB = Path(__file__).parent.parent.parent / 'Data' / 'chembl_32_bioisostere_linkers.db'
# These 2 values should match those used to make LINKER_DB.  10 and 5
# are the numbers suggested by Ertl et al.
MAX_HEAVIES=10
MAX_BONDS=5


def build_output_table(new_mols: list[Chem.Mol], parent_mols: list[Chem.Mol]) -> TableData:
    parent_col = molecules_to_column(parent_mols, "Parent Mols",
                                      DataType.BINARY)
    new_col = molecules_to_column(new_mols, "New Linker Mols",
                                  DataType.BINARY)
    table_data = TableData(tableName='New Linker Mols',
                           columns=[parent_col, new_col])
    return table_data

class LinkerReplacements(DataFunction):

    def __init__(self) -> None:
        self._plus_delta_bonds = -1
        self._minus_delta_bonds = -1
        self._match_hbonds = False
        self._structure_column_name = None
        self._parent_mols = None
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
        # column_id = string_input_field(request, 'idColumn')
        # id_column = request.inputColumns[column_id]
        # print(id_column.values)

    def do_replacements(self) -> tuple[list[Chem.Mol], list[Chem.Mol]]:
        """
        Replace all the linkers in the input set, returning a list of
        new molecules and a corresponding list of the input molecules.
        Note that the input molecules aren't copied, just the reference
        to them, so multiple entries in parent_mols may be pointing to
        the same object.
        """
        all_new_mols = []
        parent_mols = []
        for mol in self._parent_mols:
            if mol is None or not mol:
                continue
            print(f'doing {Chem.MolToSmiles(mol)}')
            new_mols, query_cp = \
                replace_linkers(mol, LINKER_DB, max_heavies=MAX_HEAVIES,
                                max_bonds=MAX_BONDS,
                                plus_length=self._plus_delta_bonds,
                                minus_length=self._minus_delta_bonds,
                                match_donors=self._match_hbonds,
                                match_acceptors=self._match_hbonds,
                                max_mols_per_input=self._max_mols_per_input)
            print(f'number of replacements : {len(new_mols)}')
            all_new_mols.extend(new_mols)
            parent_mols.extend([query_cp] * len(new_mols))

        return all_new_mols, parent_mols

    def execute(self, request: DataFunctionRequest) -> DataFunctionResponse:

        self.extract_input_data(request)
        new_mols, parent_mols = self.do_replacements()
        output_table = build_output_table(new_mols, parent_mols)

        return DataFunctionResponse(outputTables=[output_table])
