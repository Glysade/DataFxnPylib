import re
from collections import defaultdict
from itertools import product
from json import load
from pathlib import Path
from typing import Any, Optional

from rdkit import Chem
from rdkit.Chem import rdDepictor, rdmolops
from rdkit.Chem.rdchem import MolSanitizeException

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


def make_core_and_rgroups(mol: Chem.Mol, core_query: Chem.Mol) -> list[tuple[Chem.Mol, Chem.Mol, list[str]]]:
    """
    Use the core_query to do an R Group strip on mol.  That involves
    matching the core_query onto the molecule, and removing any
    substituents off it, using Chem.FragmentOnBonds.  If the
    core_query matches more than once, all possibilities are returned.
    For these purposes, we can discount symmetrical matches of the same
    atoms.
    Returns a list of the cores as a molecule, a copy of the parent
    molecule and the R Groups as SMILES strings, with atom mappings set
    up so molzip can be used conveniently.  The core and parent
    molecules will have the atoms that matched core_query flagged with
    property _GL_CORE_ and the sequence number of the core atom it
    matched.  These are used in alignment and colouring later.

    :param mol:
    :param core_query:
    :return:
    """
    ret_mols = []
    matches = mol.GetSubstructMatches(core_query)
    for match in matches:
        core_map = []
        # take a copy so markings of core are preserved
        mol_cp = Chem.Mol(mol)
        bonds_to_go = []
        dummy_labels = []
        for i, match_at in enumerate(match):
            at = mol_cp.GetAtomWithIdx(match_at)
            at.SetProp('_GL_CORE_', f'{i}')
            core_map.append((i, at.GetIdx()))
            for bond in mol_cp.GetAtomWithIdx(match_at).GetBonds():
                other_atom = bond.GetOtherAtomIdx(match_at)
                if other_atom not in match:
                    btg_num = len(bonds_to_go)
                    dummy_labels.append((btg_num + 1001, btg_num + 1001))
                    bonds_to_go.append(bond.GetIdx())
        rdDepictor.GenerateDepictionMatching2DStructure(mol_cp, core_query,
                                                        core_map)
        frag_mol = Chem.FragmentOnBonds(mol_cp, bonds_to_go,
                                        dummyLabels=dummy_labels)
        # the fragmentation labels the dummy atoms at the break points
        # with isotope numbers.  Convert these to atom map numbers so
        # that molzip can be used
        for atom in frag_mol.GetAtoms():
            if atom.GetIsotope() > 1000:
                atom.SetAtomMapNum(atom.GetIsotope())
                atom.SetIsotope(0)
        # Don't sanitize the fragments as they may contain partial
        # aromatic rings that will throw an exception.  They will
        # be bidentate (by definition), so added back on below and
        # the problem avoided.
        frags = Chem.GetMolFrags(frag_mol, asMols=True,
                                 sanitizeFrags=False)

        frag_core = None
        r_group_smis = []
        bidentates = []
        stars = re.compile(r'(\[\*:\d+\])')
        for frag in frags:
            for atom in frag.GetAtoms():
                if atom.HasProp('_GL_CORE_'):
                    frag_core = frag
                    break
            if frag != frag_core:
                smi = Chem.MolToSmiles(frag)
                num_stars = len(stars.findall(smi))
                if num_stars > 1:
                    bidentates.append(frag)
                else:
                    # Non-bidentate fragments can now be sanitized
                    # which will make the lookup later more reliable.
                    # If they won't sanitize, keep them for now but it
                    # is more than likely that they will not have any
                    # replacements so will be put straight back on
                    # again.
                    try:
                        Chem.SanitizeMol(frag)
                        smi = Chem.MolToSmiles(frag)
                    except MolSanitizeException:
                        pass
                    r_group_smis.append(smi)

        # If we have bidentate substituents, add them back to the core.
        # This may generate an error "Incomplete atom labelling,
        # cannot make bond" which is (I think) molzip complaining that
        # not all atom mappings in frag_core have a match in bid, which
        # is expected because we are only repairing the bidentate
        # removal and leaving the rest of the fragmentation for the r
        # group replacement.
        if bidentates:
            for bid in bidentates:
                frag_core = Chem.CombineMols(frag_core, bid)
            frag_core = Chem.molzip(frag_core)
        ret_mols.append((frag_core, mol_cp, r_group_smis))
    return ret_mols


def make_rgroup_lookup_smi(mol: Chem.Mol) -> tuple[str, int]:
    """
    Take the RGroup, which is expected to have atom mappings on the
    dummy atoms, and return a new SMILES with the first atom mapping
    removed and its atom map number.  There will sometimes
    be polydentate RGroups which won't have a replacement in the table
    so don't need to be handled - it's fine if they are left with
    dangling atom maps as they will be placed straight back onto the
    core.

    """
    atom_map_num = 0
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 0:
            atom_map_num = atom.GetAtomMapNum()
            atom.SetAtomMapNum(0)
            break
    return Chem.MolToSmiles(mol), atom_map_num


def make_mapped_rgroups(smis: list[str], atom_map_num: int, rgroup_smi: str,
                        layer_number: str) -> list[Chem.Mol]:
    """
    Take the SMIlES strings, create a Mol for each one and add the atom
    map num to the dummy atom.  If the SMILES for the R Group isn't the
    same as the original, flog the atoms with the property _GL_R_GROUP_
    plus layer_number for subsequent coloured rendering.
    :param smis:
    :param atom_map_num:
    :param rgroup_smi:
    :param layer_number:
    :return:
    """
    ret_mols = []
    prop_name = '_GL_R_GROUP_' + layer_number + '_'
    for smi in smis:
        map_smi = smi.replace('*', f'[*:{atom_map_num}]')
        frag_mol = Chem.MolFromSmiles(map_smi)
        if smi != rgroup_smi:
            for at in frag_mol.GetAtoms():
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
    the layer1 and layer2 replacements as requested.  If rgroup_smi
    isn't in substs_data, just returns rgroup_smi suitably mackled so
    that no replacement is performed in the end.  That's a bit
    inefficient, but probably good enough.

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
        except KeyError:
            layer1_smis = [rgroup_smi]
        layer1_mols = make_mapped_rgroups(layer1_smis, atom_map_num,
                                          rgroup_smi, '1')
    if use_layer2:
        try:
            layer2_smis = substs_data[rgroup_smi][1]
        except KeyError:
            layer2_smis = [rgroup_smi]
        layer2_mols = make_mapped_rgroups(layer2_smis, atom_map_num,
                                          rgroup_smi, '2')

    return layer1_mols, layer2_mols


def make_analogues(core_and_rgroups: list[tuple[Chem.Mol, list[str]]],
                   substs_data: dict[str: tuple[list[str], list[str]]],
                   use_layer1: bool, use_layer2: bool,
                   include_orig_rgroup: bool) -> list[Chem.Mol, Chem.Mol]:
    """
    Take the list of core and r groups for the molecule, and make
    analogues of each using the relevant replacement in substs_data.
    If include_orig_rgroup is True, the substitutions at each point
    will include the original R Group, so that combinations are
    produced that include no change at a position.  The case where
    none of the R Groups changed so the original molecule is
    returned will be weeded out later.
    Returns a list of tuples, with the analogue and its parent.  Two
    analogues can have the same parent structure, but be based on
    different matches of the core onto the parent, and hence different
    parent structures - the difference being what's in the _GL_CORE_
    properties on the matching atoms.
    :param core_and_rgroups:
    :param substs_data:
    :param use_layer1:
    :param usr_layer2:
    :param include_orig_rgroup:
    :return:
    """
    analogues = []
    for core, parent, rgroups in core_and_rgroups:
        rgroup_repls = []
        for rgroup in rgroups:
            rgroup_lookup, atom_map_num = make_rgroup_lookup_smi(rgroup)
            layer1_mols, layer2_mols = \
                make_rgroups_for_substs(rgroup_lookup, atom_map_num,
                                        substs_data, use_layer1, use_layer2)
            if include_orig_rgroup:
                layer1_mols.append(Chem.MolFromSmiles(rgroup))
            rgroup_repls.append(layer1_mols + layer2_mols)
        for substs in product(*rgroup_repls):
            analogue = Chem.Mol(core)
            for s in substs:
                analogue = Chem.CombineMols(analogue, s)
            analogue = Chem.molzip(analogue)
            analogues.append((analogue, parent))

    return analogues


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
    Use the information in atom props _GL_CORE_ to align the analogue
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
            core_map.append((int(at.GetProp('_GL_CORE_')), at.GetIdx()))
        except KeyError:
            pass
        if at.HasProp('_GL_R_GROUP_1_'):
            layer_1_ats.append(at.GetIdx())
        if at.HasProp('_GL_R_GROUP_2_'):
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


def replace_rgroups(mols: list[Chem.Mol], ids: list[Any],
                    id_type: DataType, core_query: Chem.Mol,
                    use_layer1: bool, use_layer2: bool,
                    include_orig_rgroup: bool,
                    input_column_name: str) -> tuple[TableData, list[int]]:
    """
    Take the list of molecules, use the core_query to define R Groups
    of it, and then create a table of analogues where the R Groups are
    replaced using Bajorath's table.  Returns a table where each row is
    a new analogue, with the first column the parent molecule.
    :param mols: list of molecules from which analogues are to be
                 derived ([Chem.mol, ]).
    :param ids: list of names of molecules from input table column
    :param id_type: type of the ID data
    :param core_query: query molecule defining the core.  Not all
                       molecules must match, but clearly only those
                       that do will produce analogues (Chem.Mol).
    :param use_layer1: Use the layer 1 replacements in the table (bool)
    :param use_layer2: Use the layer 1 replacements in the table (bool)
    :param include_orig_rgroup: include the original R Group in the
                                substitution list i.e. include products
                                with no change at each position.
    :param input_column_name: (str)
    :return: The table of analogues (TableData)
    """
    substs_data = load_replacements_data()

    rdDepictor.SetPreferCoordGen(True)
    rdDepictor.Compute2DCoords(core_query)

    all_analogues = []
    analogue_parents = []
    analogue_parent_ids = []
    analogue_count = {}
    input_smiles = []
    for mol, id in zip(mols, ids):
        if id not in analogue_count:
            analogue_count[id] = 0
        if mol is not None and mol and mol.HasSubstructMatch(core_query):
            input_smiles.append(Chem.MolToSmiles(mol))
            rdDepictor.Compute2DCoords(mol)
            core_and_rgroups = make_core_and_rgroups(mol, core_query)
            analogues = make_analogues(core_and_rgroups, substs_data,
                                       use_layer1, use_layer2,
                                       include_orig_rgroup)
            parent_smi = Chem.MolToSmiles(mol)
            for analogue in analogues:
                if Chem.MolToSmiles(analogue[0]) != parent_smi:
                    all_analogues.append(analogue[0])
                    analogue_parents.append(analogue[1])
                    analogue_parent_ids.append(id)

    analogue_smiles = defaultdict(int)
    final_analogues = []
    final_parents = []
    final_parent_ids = []
    input_smiles_set = set(input_smiles)
    for analogue, parent, parent_id in zip(all_analogues, analogue_parents, analogue_parent_ids):
        smi = Chem.MolToSmiles(analogue)
        if not analogue_smiles[smi] and smi not in input_smiles_set:
            align_analogue_to_core(analogue, core_query)
            final_analogues.append(analogue)
            final_parents.append(parent)
            final_parent_ids.append(parent_id)
            analogue_count[parent_id] += 1

        analogue_smiles[smi] += 1

    table_name = f'Analogues of {input_column_name}'
    parent_col = molecules_to_column(final_parents, f'Parent {input_column_name}',
                                     DataType.BINARY)
    parent_ids_col = ColumnData(name='Parent IDs', dataType=id_type,
                                values=final_parent_ids)
    analogue_col = molecules_to_column(final_analogues, 'Analogues',
                                       DataType.BINARY)
    analogue_count_col = [analogue_count[id] for id in ids]
    return TableData(tableName=table_name, columns=[parent_col, parent_ids_col,
                                                    analogue_col]), analogue_count_col


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


def build_analogues(core: Chem.Mol, rgroup_line: list[Chem.Mol],
                    substs_data: dict[str: tuple[list[str], list[str]]],
                    use_layer1: bool, use_layer2: bool,
                    include_orig_rgroup:bool) -> list[Chem.Mol]:
    """
    Build all the analogues for a single core and set of R Groups.
    Returns a list of unique molecules.
    """
    # The core and rgroups come in with the attachment points as dummy
    # atoms with isotope labels.  For historical reasons it's easier if
    # these are replaced by atom maps
    new_core = isotopes_to_atommaps(core)
    new_rgroups = [isotopes_to_atommaps(rg) for rg in rgroup_line]
    analogues = []
    rgroup_repls = []
    for rgroup in new_rgroups:
        if rgroup is None:
            continue
        rgroup_lookup, atom_map_num = make_rgroup_lookup_smi(rgroup)
        layer1_mols, layer2_mols = \
            make_rgroups_for_substs(rgroup_lookup, atom_map_num,
                                    substs_data, use_layer1, use_layer2)
        # If there are no replacements for the R Group, set it up so
        # it will go straight back on.  Use a layer number of 0 so it
        # isn't highlighted.
        if include_orig_rgroup or (not layer1_mols and not layer2_mols):
            orig_rgroup = make_mapped_rgroups([Chem.MolToSmiles(rgroup)],
                                              atom_map_num, '', '0')
            layer1_mols.append(orig_rgroup[0])
        rgroup_repls.append(layer1_mols + layer2_mols)

    for substs in product(*rgroup_repls):
        analogue = Chem.Mol(new_core)
        for i, anat in enumerate(analogue.GetAtoms()):
            anat.SetProp('_GL_CORE_', str(i))
        for s in substs:
            analogue = Chem.CombineMols(analogue, s)
        # Note that if an RGroup was an H (rgroup_lookup = *[H]),
        # molzip will complain but will do things correctly. Look out
        # for
        # WARNING: not removing hydrogen atom with dummy atom neighbors
        analogue = rdmolops.RemoveHs(Chem.molzip(analogue))
        analogues.append(analogue)

    return analogues


def build_all_analogues(parent_mols: list[Chem.Mol], parent_ids: list[Any],
                        cores: list[Chem.Mol], rgroups: list[list[Chem.Mol]],
                        id_type: DataType, use_layer1: bool, use_layer2: bool,
                        include_orig_rgroup: bool,
                        input_column_name: str) -> tuple[TableData, ColumnData]:
    """
    Takes a set of cores and RGroups and makes all possible
    substitutions for each R Group using the substitution table and
    the appropriate layers.  Each list (outer list for rgroups) should
    be the same length and the same index in each corresponds to the
    same input mol.
    Returns a new table with the analogues, and a column to be added to
    the original table giving the number of analogues each parent
    produced. Analogues that are the same as any parent aren't returned
    and if a parent produces no analogues, nothing is added to the
    output table for it.
    The analogues are aligned with the core and colour coded to show
    the core and substititions.
    :param parent_mols:
    :param parent_ids:
    :param cores:
    :param rgroups:
    :param use_layer1: Use the layer 1 replacements in the table (bool)
    :param use_layer2: Use the layer 1 replacements in the table (bool)
    :param include_orig_rgroup: include the original R Group in the
                                substitution list i.e. include products
                                with no change at each position.
    :param input_column_name: (str)
    :return: The table of analogues (TableData) and a column
    (ColumnData) to be added to the original column.
    """
    substs_data = load_replacements_data()

    rdDepictor.SetPreferCoordGen(True)

    all_analogues_by_smi = {}
    all_analogues = []
    analogue_counts = defaultdict(int)

    parent_smis = set([Chem.MolToSmiles(pm) for pm in parent_mols])
    parents_dict = {}
    for parent, parent_id, core, rgroup_line in \
            zip(parent_mols, parent_ids, cores, rgroups):
        # Build the parent back up from the core and original R Groups,
        # so it can be aligned and the core coloured for easy
        # comparison to the analogues.
        new_parent = build_analogues(core, rgroup_line, substs_data, False, False,
                                     True)
        align_analogue_to_core(new_parent[0], core)
        parents_dict[parent_id] = new_parent[0]
        analogues = build_analogues(core, rgroup_line, substs_data, use_layer1,
                                    use_layer2, include_orig_rgroup)
        for an in analogues:
            an_smi = Chem.MolToSmiles(an)
            if an_smi not in parent_smis and an_smi not in all_analogues_by_smi:
                align_analogue_to_core(an, core)
                all_analogues_by_smi[an_smi] = (parent_id, an)
                analogue_counts[parent_id] += 1
                all_analogues.append((an, parent_id))

    analogue_count_vals = [analogue_counts[id] for id in parent_ids]
    analogue_count_col = ColumnData(name=f'Num R Group Subs {input_column_name}',
                                    dataType=DataType.INTEGER,
                                    values=analogue_count_vals)

    analogue_parents = []
    analogue_parent_ids = []
    analogue_col_vals = []
    for an_smi, (parent_id, analogue) in all_analogues_by_smi.items():
        analogue_parents.append(parents_dict[parent_id])
        analogue_parent_ids.append(parent_id)
        analogue_col_vals.append(analogue)

    parent_col = molecules_to_column(analogue_parents,
                                     f'Parent {input_column_name}',
                                     DataType.BINARY)
    parent_ids_col = ColumnData(name='Parent IDs', dataType=id_type,
                                values=analogue_parent_ids)
    analogue_col = molecules_to_column(analogue_col_vals, 'Analogues',
                                       DataType.BINARY)

    table_name = f'Analogues of {input_column_name}'
    table_data = TableData(tableName=table_name,
                           columns=[parent_col, parent_ids_col, analogue_col])

    return table_data, analogue_count_col


class RGroupReplacement(DataFunction):

    def execute(self, request: DataFunctionRequest) -> DataFunctionResponse:
        column_id = string_input_field(request, 'structureColumn')
        input_column = request.inputColumns[column_id]
        parent_mols = column_to_molecules(input_column)
        column_id = string_input_field(request, 'idColumn')
        id_column = request.inputColumns[column_id]
        parent_ids = id_column.values
        ids_type = id_column.dataType

        cores_column_id = string_input_field(request, 'coresColumn')
        cores_column = request.inputColumns[cores_column_id]
        cores = column_to_molecules(cores_column)

        rgroup_column_ids = string_list_input_field(request, 'rGroupColumns')
        rgroups = []
        for rgroup_col_id in rgroup_column_ids:
            rgroup_col = request.inputColumns[rgroup_col_id]
            rgroups.append(column_to_molecules(rgroup_col))

        # pivot the R Group columns so that each row is the R Groups for an
        # input molecule, core
        rgroup_lines = list(zip(*rgroups))

        use_layer1 = boolean_input_field(request, 'useLayer1')
        use_layer2 = boolean_input_field(request, 'useLayer2')
        include_orig_rgroup = boolean_input_field(request, 'incOrigRGroups')
        analogues_table, analogue_count_col =\
            build_all_analogues(parent_mols, parent_ids, cores, rgroup_lines,
                                ids_type, use_layer1, use_layer2,
                                include_orig_rgroup, input_column.name)
        response = DataFunctionResponse(outputTables=[analogues_table],
                                        outputColumns=[analogue_count_col])
        return response
