import re
from collections import defaultdict
from itertools import product
from json import load
from pathlib import Path
from typing import Any, Optional

from rdkit import Chem, rdBase
from rdkit.Chem import rdDepictor, rdmolops

from df.chem_helper import column_to_molecules, input_field_to_molecule, \
    molecules_to_column
from df.data_transfer import DataFunction, DataFunctionRequest, \
    DataFunctionResponse, ColumnData, DataType, \
    TableData, boolean_input_field, string_input_field, string_list_input_field


def load_replacements_data() -> dict[str: tuple[list[str], list[str]]]:
    """
    Load the replacements data from the JSON file.  Returns a dict
    keyed on the SMILES of the thing to be replaced, containing a tuple
    with the lists of SMILES of first and second layer replacements.
    :return:
    """
    im_here = Path(__file__)
    substs_file = im_here.parent.parent.parent / 'Data' / 'top500_R_replacements.json'
    with open(substs_file, 'r') as f:
        raw_substs_data = load(f)

    # it's a lot more convenient in use if the outer list is turned
    # into a dict.
    substs_data = {}
    for raw_subst in raw_substs_data:
        smi = raw_subst['@SMILES']
        layer1 = []
        layer2 = []
        # a vagary of the conversion from the original XML
        # to the JSON used here is that if the layer
        # only had 1 entry, it's not put in a list, just
        # into a dict
        first_layers = raw_subst['first_layer']
        if isinstance(first_layers, dict):
            first_layers = [first_layers]
        for first_layer in first_layers:
            layer1.append(first_layer['@SMILES'])
            try:
                second_layers = first_layer['second_layer']
                if isinstance(second_layers, dict):
                    second_layers = [second_layers]
                for second_layer in second_layers:
                    layer2.append(second_layer['@SMILES'])
            except KeyError:
                # not all substituents have second layers
                pass
        substs_data[smi] = (layer1, layer2)

    return substs_data


def make_rgroup_lookup_smi(mol: Chem.Mol) -> tuple[str, int, bool]:
    """
    Take the RGroup, which is expected to have atom mappings on the
    dummy atoms, and return a new SMILES with the first atom mapping
    removed and its atom map number.  There will sometimes
    be polydentate RGroups which won't have a replacement in the table
    so don't need to be handled - it's fine if they are left with
    dangling atom maps as they will be placed straight back onto the
    core.
    Removes the atom map num from molecule, so the molecule is
    altered by the process.
    """
    atom_map_num = 0
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 0:
            atom_map_num = atom.GetAtomMapNum()
            atom.SetAtomMapNum(0)
            break
    smi = Chem.MolToSmiles(mol)
    num_stars = sum(1 for i in range(len(smi)) if smi[i] == '*')
    return smi, atom_map_num, bool(num_stars > 1)


def make_mapped_rgroups(smis: list[str], atom_map_num: int, rgroup_smi: str,
                        layer_number: str) -> list[Chem.Mol]:
    """
    Take the SMIlES strings, create a Mol for each one and add the atom
    map num to the dummy atom.  If the SMILES for the R Group isn't the
    same as the original, flog the atoms with the property GLYS_R_GROUP_
    plus layer_number for subsequent coloured rendering.
    :param smis:
    :param atom_map_num:
    :param rgroup_smi:
    :param layer_number:
    :return:
    """
    ret_mols = []
    prop_name = 'GLYS_R_GROUP_' + layer_number
    for smi in smis:
        # The R SMILES '*[H]' gives an irritating warning
        # 'not removing hydrogen atom with dummy atom neighbors'
        # and it crops up a lot.  Turn warnings off temporarily.
        if smi == '*[H]':
            rdBase.DisableLog('rdApp.warning')
        frag_mol = Chem.MolFromSmiles(smi)
        if smi == '*[H]':
            rdBase.EnableLog('rdApp.warning')
        if smi != rgroup_smi:
            for at in frag_mol.GetAtoms():
                if at.GetAtomicNum() == 0 and not at.GetAtomMapNum():
                    at.SetAtomMapNum(atom_map_num)
                at.SetProp(prop_name, f'{at.GetIdx()}')
        ret_mols.append(frag_mol)

    return ret_mols


def make_rgroups_for_substs(rgroup_smi: str, atom_map_num: int,
                            substs_data: dict[str: tuple[list[str], list[str]]],
                            use_layer1: bool, use_layer2: bool) -> tuple[list[Chem.Mol], list[Chem.Mol]]:
    """
    Lookup the replacements for rgroup_smi in substs_data and build
    molecules from them with the atom_map_num added to the dummy atom
    so they are ready to go straight into molzip.  Returns lists of
    the layer1 and layer2 replacements as requested, which might
    both be empty if rgroup_smi isn't in substs_data.

    :param rgroup_smi:
    :param atom_map_num:
    :param substs_data:
    :param use_layer1:
    :param use_layer2:
    :return:
    """
    layer1_mols = []
    layer2_mols = []
    if use_layer1:
        try:
            layer1_smis = substs_data[rgroup_smi][0]
            layer1_mols = make_mapped_rgroups(layer1_smis, atom_map_num,
                                              rgroup_smi, '1')
        except KeyError:
            pass
    if use_layer2:
        try:
            layer2_smis = substs_data[rgroup_smi][1]
            layer2_mols = make_mapped_rgroups(layer2_smis, atom_map_num,
                                              rgroup_smi, '2')
        except KeyError:
            pass

    return layer1_mols, layer2_mols


def bonds_between_atoms(mol: Chem.Mol, atoms: list[int]) -> list[int]:
    """
    Return the bonds in the molecule between the atoms whose indices
    are supplied.
    """
    bonds = []
    for a1 in atoms:
        for a2 in atoms:
            if a1 > a2:
                bond = mol.GetBondBetweenAtoms(a1, a2)
                if bond is not None:
                    bonds.append(bond.GetIdx())
    return bonds


def align_analogue_to_core(analogue: Chem.Mol, core_query: Chem.Mol) -> None:
    """
    Use the information in atom props GLYS_CORE to align the analogue
    so that corresponding atoms are on top of the parent, and also
    highlight the analogue atoms and bonds of the core.  The parent
    is altered to match query_mol, so everything lines up from parent
    to parent and analogue to analogue.
    Colour core blue, layer 1 r groups red, layer 2 r groups orange.
    My father, who is very red/green colourblind, says he finds these
    colours easy to distinguish.
    :param analogue:
    :param parent:
    :param query_mol
    :return:
    """
    core_map = []
    layer_1_ats = []
    layer_2_ats = []
    for at in analogue.GetAtoms():
        try:
            core_map.append((int(at.GetProp('GLYS_CORE')), at.GetIdx()))
        except KeyError:
            pass
        if at.HasProp('GLYS_R_GROUP_1'):
            layer_1_ats.append(at.GetIdx())
        if at.HasProp('GLYS_R_GROUP_2'):
            layer_2_ats.append(at.GetIdx())
    rdDepictor.GenerateDepictionMatching2DStructure(analogue, core_query,
                                                    atomMap=core_map)
    core_map.sort(key=lambda p: p[0])
    core_ats = [a[1] for a in core_map]
    core_bonds = bonds_between_atoms(analogue, core_ats)
    layer_1_bonds = bonds_between_atoms(analogue, layer_1_ats)
    layer_2_bonds = bonds_between_atoms(analogue, layer_2_ats)

    # add bonds between core and r group atoms
    for ca in core_ats:
        for rga in layer_1_ats:
            bond = analogue.GetBondBetweenAtoms(ca, rga)
            if bond is not None:
                layer_1_bonds.append(bond.GetIdx())
        for rga in layer_2_ats:
            bond = analogue.GetBondBetweenAtoms(ca, rga)
            if bond is not None:
                layer_2_bonds.append(bond.GetIdx())

    core_bonds_str = ' '.join([str(b + 1) for b in core_bonds])
    prop_text = f'COLOR #00bfff\nATOMS\nBONDS {core_bonds_str}'
    if layer_1_ats:
        layer_1_bonds_str = ' '.join([str(b + 1) for b in layer_1_bonds])
        prop_text += f'\nCOLOR #dc143c\nATOMS\nBONDS {layer_1_bonds_str}'
    if layer_2_ats:
        layer_2_bonds_str = ' '.join([str(b + 1) for b in layer_2_bonds])
        prop_text += f'\nCOLOR #ffbf00\nATOMS\nBONDS {layer_2_bonds_str}'
    analogue.SetProp('Renderer_Highlight', prop_text)


def isotopes_to_atommaps(mol: Chem.Mol) -> Optional[Chem.Mol]:
    """
    Take the molecule and return a copy with all the isotope-labelled
    dummy atoms replaced by atom-mapped dummy atoms of 0 isotope.
    """
    if mol is None:
        return None

    mol_cp = Chem.Mol(mol)
    for at in mol_cp.GetAtoms():
        if at.GetAtomicNum() == 0:
            at.SetAtomMapNum(at.GetIsotope())
            at.SetIsotope(0)
    return mol_cp


def assemble_molecule(core: Chem.Mol, rgroups: tuple[Chem.Mol],
                      polydentate: bool, rgroup_col_names: list[str]) -> Chem.Mol:
    """
    Take the core and r groups and make them into a molecule, labelling
    the core with GLYS_CORE props on the atoms.  Add a prop
    GLYS_CHANGED_R_GROUPS listing the r groups that were changed from
    the original value for that positition.
    """
    mol = Chem.Mol(core)
    for i, anat in enumerate(mol.GetAtoms()):
        anat.SetProp('GLYS_CORE', str(i))

    rgroups_smis = []
    changed_r_group_list = []
    for rg, rg_col_name in zip(rgroups, rgroup_col_names):
        try:
            _ = rg.GetProp('GLYS_ORIG_R_GROUP')
        except KeyError:
            changed_r_group_list.append(rg_col_name)

        if polydentate:
            rg_smi = Chem.MolToSmiles(rg)
        else:
            rg_smi = 'dummy'
        if rg_smi not in rgroups_smis:
            Chem.AssignStereochemistry(rg)
            mol = Chem.CombineMols(mol, rg)
            if polydentate:
                rgroups_smis.append(rg_smi)

    mol = Chem.molzip(mol)
    mol.SetProp('GLYS_CHANGED_R_GROUPS', ':'.join(changed_r_group_list))
    mol = rdmolops.RemoveHs(mol)
    return mol


def build_analogues(core: Chem.Mol, rgroups_line: list[Chem.Mol],
                    rgroup_col_names: list[str],
                    substs_data: dict[str: tuple[list[str], list[str]]],
                    use_layer1: bool, use_layer2: bool) -> list[Chem.Mol]:
    """
    Build all the analogues for a single core and set of R Groups.
    Returns a list of unique molecules.
    """
    # The core and rgroups come in with the attachment points as dummy
    # atoms with isotope labels.  For historical reasons it's easier if
    # these are replaced by atom maps
    new_core = isotopes_to_atommaps(core)
    rgroups = [isotopes_to_atommaps(rg) for rg in rgroups_line]

    analogues = []
    rgroup_repls = []
    polydentate = False
    for rgroup, rgroup_col_name in zip(rgroups, rgroup_col_names):
        if rgroup is None:
            continue
        rgroup_lookup, atom_map_num, pd = make_rgroup_lookup_smi(rgroup)
        if pd:
            polydentate = True
        layer1_mols, layer2_mols = \
            make_rgroups_for_substs(rgroup_lookup, atom_map_num,
                                    substs_data, use_layer1, use_layer2)
        # Add in the original R Group so we also get products with no
        # change at each position. Use a layer number of 0 so it isn't
        # highlighted. rgroup is altered by make_group_lookup_smi, so
        # we need a new copy.  Flag the copy as being the original
        # r group for later use.
        orig_rgroup = make_mapped_rgroups([Chem.MolToSmiles(rgroup)],
                                          atom_map_num, '', '0')
        orig_rgroup[0].SetProp(f'GLYS_ORIG_R_GROUP', rgroup_col_name)
        layer1_mols.append(orig_rgroup[0])
        rgroup_repls.append(layer1_mols + layer2_mols)

    for substs in product(*rgroup_repls):
        analogues.append(assemble_molecule(new_core, substs, polydentate,
                                           rgroup_col_names))

    return analogues


def rebuild_parents(cores: list[Chem.Mol], rgroups: list[list[Chem.Mol]],
                    rgroup_col_names: list[str]) -> list[Optional[Chem.Mol]]:
    """
    Take the cores and rgroups and rebuild the parents, aligned on the
    core and coloured.  Because the RGD loses chirality information at
    the centres the rgroups are attached to, the resultant molecules
    may not be identical to the original parent. For example,
    CONC(=O)c1cc(N2CC[C@H](NC(=O)c3[nH]c(C)c(Cl)c3Cl)C2)nc(Cl)n1
    becomes
    CONC(=O)c1cc(N2CCC(NC(=O)c3[nH]c(C)c(Cl)c3Cl)C2)nc(Cl)n1 with core
    [*:3]N1CCC(C2)NC(=O)[*:4] because the RGD loses the chirality info
    where the amide attaches to the tetrahydropyrrolidine.
    If the core is None, the rebuilt parent will also be None.
    """
    # Since we're not doing R Group substitution, we don't need the
    # real substitutions database.
    dummy_substs_data = {}
    rebuilt_parents = []
    for core, rgroup_line in zip(cores, rgroups):
        if core is None:
            rebuilt_parents.append(None)
        else:
            new_parent = build_analogues(core, rgroup_line, rgroup_col_names,
                                         dummy_substs_data, False, False)
            align_analogue_to_core(new_parent[0], core)
            rebuilt_parents.append(new_parent[0])
    return rebuilt_parents


def build_all_analogues(parent_mols: list[Chem.Mol], parent_ids: list[Any],
                        cores: list[Chem.Mol], rgroups: list[list[Chem.Mol]],
                        rgroup_col_names: list[str],
                        use_layer1: bool, use_layer2: bool) \
        -> tuple[dict[str: tuple[str, Chem.Mol, Chem.Mol]],
                 list[int],
                 dict[str: Chem.Mol]]:
    """
    Make all the analogues from the cores and rgroups.
    Returns a tuple containing:
        a dict keyed on the analogue SMILES that holds the analogue's
        parent id, the analogue molecule itself, and the core used;
        a list of the number of analogues each parent made;
        a dict, keyed on parent id containing a copy of the parent
        that has been constructed from the appropriate core and
        input r groups, aligned on the core and with details to colour
        the core.
    """
    substs_data = load_replacements_data()

    rdDepictor.SetPreferCoordGen(True)

    all_analogues_by_smi = {}
    analogue_counts = defaultdict(int)
    rebuilt_parents = rebuild_parents(cores, rgroups, rgroup_col_names)

    parent_smis_list = []
    for pm in parent_mols:
        Chem.AssignStereochemistry(pm)
        parent_smis_list.append(Chem.MolToSmiles(pm))
    parent_smis = set(parent_smis_list)
    rebuilt_parent_smis = set([Chem.MolToSmiles(pm) for pm in rebuilt_parents if pm is not None])
    parents_dict = {}
    # Give each core a unique number
    cores_dict = {}

    for rb_parent, parent_id, core, rgroups_line in \
            zip(rebuilt_parents, parent_ids, cores, rgroups):
        if core is None:
            continue
        core_smi = Chem.MolToSmiles(core)
        try:
            core_num = cores_dict[core_smi]
        except KeyError:
            core_num = len(cores_dict) + 1
            cores_dict[core_smi] = core_num
        core.SetProp('GLYS_CORE_NUM', f'{core_num}')

        parents_dict[parent_id] = rb_parent
        analogues = build_analogues(core, rgroups_line, rgroup_col_names,
                                    substs_data, use_layer1, use_layer2)
        for an in analogues:
            an_smi = Chem.MolToSmiles(an)
            if (an_smi not in parent_smis and an_smi not in rebuilt_parent_smis
                    and an_smi not in all_analogues_by_smi):
                align_analogue_to_core(an, core)
                all_analogues_by_smi[an_smi] = (parent_id, an, core)
                analogue_counts[parent_id] += 1

    analogue_count_vals = [analogue_counts[id] for id in parent_ids]
    return all_analogues_by_smi, analogue_count_vals, parents_dict


def build_output_objects(all_analogues_by_smi: dict[str: tuple[str, Chem.Mol, Chem.Mol]],
                         analogue_count_vals: list[int],
                         parents_dict: dict[str, Chem.Mol],
                         input_column_name: str,
                         id_type: DataType) -> tuple[TableData, ColumnData]:
    """
    Takes the output from build_all_analogues and makes the output objects
    for return.  The analogues go into a table, and an extra column is made
    of the analogues each parent produced that is to go into the calling
    table.
    """
    analogue_count_col = ColumnData(name=f'Num R Group Subs {input_column_name}',
                                    dataType=DataType.INTEGER,
                                    values=analogue_count_vals)

    analogue_parents = []
    analogue_parent_ids = []
    analogue_col_vals = []
    analogue_cores = []
    analogue_changed_r_groups = []
    core_numbers = []
    rgc_prop_name = 'GLYS_CHANGED_R_GROUPS'
    cn_prop_name = 'GLYS_CORE_NUM'
    for an_smi, (parent_id, analogue, core) in all_analogues_by_smi.items():
        analogue_parents.append(parents_dict[parent_id])
        analogue_parent_ids.append(parent_id)
        analogue_col_vals.append(analogue)
        analogue_cores.append(core)
        analogue_changed_r_groups.append(analogue.GetProp(rgc_prop_name))
        core_numbers.append(analogue.GetProp(cn_prop_name))

    parent_col = molecules_to_column(analogue_parents,
                                     f'Parent {input_column_name}',
                                     DataType.BINARY)
    parent_ids_col = ColumnData(name='Parent ID', dataType=id_type,
                                values=analogue_parent_ids)
    cores_col = molecules_to_column(analogue_cores, 'Core', DataType.BINARY)
    # Leave the core numbers as strings as that gives a more convenient
    # filter in Spotfire.
    core_nums_col = ColumnData(name='Core Number', dataType=DataType.STRING,
                               values=core_numbers)
    analogue_col = molecules_to_column(analogue_col_vals, 'Analog',
                                       DataType.BINARY)
    analogue_changes_col = ColumnData(name='Changed R Groups',
                                      dataType=DataType.STRING,
                                      values=analogue_changed_r_groups)
                                      
    #clean_col_name = input_column_name.replace('..010','>').replace('..00Y','<').replace('..00w',' ')
    #table_name = f'{clean_col_name} Analogs'
    table_name = 'R-Group Replacements'
    table_data = TableData(tableName=table_name,
                           columns=[parent_col, parent_ids_col, analogue_col, analogue_changes_col, core_nums_col, cores_col])
                                    
                                    

    return table_data, analogue_count_col


def r_group_replacement(parent_mols: list[Chem.Mol], parent_ids: list[Any],
                        cores: list[Chem.Mol], rgroups: list[list[Chem.Mol]],
                        rgroup_col_names: list[str],
                        id_type: DataType, use_layer1: bool, use_layer2: bool,
                        input_column_name: str) -> tuple[TableData, ColumnData]:
    """
    Takes a set of cores and RGroups and makes all possible
    substitutions for each R Group using the substitution table and
    the appropriate layers.  Each list (outer list for rgroups) should
    be the same length and the same index in each corresponds to the
    same parent mol.
    Returns a new table with the analogues, and a column to be added to
    the original table giving the number of analogues each parent
    produced. Analogues that are the same as any parent aren't returned
    and if a parent produces no analogues, nothing is added to the
    output table for it.
    The analogues are aligned with the core and colour coded to show
    the core and substititions.
    :param parent_mols: The original parents that the cores and r
                        groups were derived from.
    :param parent_ids: the ids corresponding to the parent_mols
    :param cores: A core for each parent, that the rgroups attach to
    :param rgroups: All the rgroups in the original parent
    :param use_layer1: Use the layer 1 replacements in the table (bool)
    :param use_layer2: Use the layer 1 replacements in the table (bool)
    :param input_column_name: (str)
    :return: The table of analogues (TableData) and a column
    (ColumnData) to be added to the original column.
    """
    all_analogues_by_smi, analogue_count_vals, parents_dict = \
        build_all_analogues(parent_mols, parent_ids, cores, rgroups,
                            rgroup_col_names, use_layer1, use_layer2)

    return build_output_objects(all_analogues_by_smi, analogue_count_vals,
                                parents_dict, input_column_name,
                                id_type)


class RGroupReplacement(DataFunction):

    def __init__(self) -> None:
        self._parent_mols = None
        self._parent_ids = None
        self._structure_column_name = None
        self._ids_type = None
        self._cores = None
        self._rgroup_lines = None
        self._rgroup_col_names = []
        self._use_layer1 = False
        self._use_layer2 = False
        super().__init__()

    def extract_input_data(self, request: DataFunctionRequest) -> None:
        column_id = string_input_field(request, 'structureColumn')
        structure_column = request.inputColumns[column_id]
        self._parent_mols = column_to_molecules(structure_column)
        self._structure_column_name = structure_column.name
        column_id = string_input_field(request, 'idColumn')
        id_column = request.inputColumns[column_id]
        self._ids_type = id_column.dataType
        self._parent_ids = id_column.values
        cores_column_id = string_input_field(request, 'coresColumn')
        cores_column = request.inputColumns[cores_column_id]
        self._cores = column_to_molecules(cores_column)

        rgroup_column_ids = string_list_input_field(request,
                                                    'rGroupColumns')
        rgroups = []
        for rgroup_col_id in rgroup_column_ids:
            rgroup_col = request.inputColumns[rgroup_col_id]
            rgroups.append(column_to_molecules(rgroup_col))
            self._rgroup_col_names.append(rgroup_col.name)

        # pivot the R Group columns so that each row is the R Groups for an
        # input molecule, core
        self._rgroup_lines = list(zip(*rgroups))
        self._use_layer1 = boolean_input_field(request, 'useLayer1')
        self._use_layer2 = boolean_input_field(request, 'useLayer2')

    def execute(self, request: DataFunctionRequest) -> DataFunctionResponse:

        self.extract_input_data(request)

        analogues_table, analogue_count_col = \
            r_group_replacement(self._parent_mols, self._parent_ids,
                                self._cores, self._rgroup_lines,
                                self._rgroup_col_names, self._ids_type,
                                self._use_layer1, self._use_layer2,
                                self._structure_column_name)
        return DataFunctionResponse(outputTables=[analogues_table],
                                    outputColumns=[analogue_count_col])
