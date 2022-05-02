""""
======
mmp.py
======

Copyright (C) 2017-2022 Glysade, LLC

Utility functions and classes for performing MMP analysis

Uses a modified version of Andrew Dalke's mmpdb rdkit software (https://github.com/rdkit/mmpdb).
Clone the modified version from https://github.com/jones-gareth/mmpdb, and install to python using setup.py.

"""
import base64
import gzip
import json
import multiprocessing
import sys
from array import array
from typing import List, NamedTuple, Union, Optional, Tuple

try:
    from mmpdblib import commandline as mmp, command_support, dbutils, analysis_algorithms
except ImportError:
    from mmpdblib import cli as mmp, dbutils, analysis_algorithms
    from mmpdblib.analysis_algorithms import get_transform_tool
from mmpdblib.analysis_algorithms import TransformResult
from mmpdblib.dbutils import open_database, DBFile
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.rdchem import Mol, RWMol, BondType, Atom, BondStereo, Bond, Conformer
from rdkit.Geometry.rdGeometry import Point3D

from ruse.rdkit.rdkit_utils import RDKitFormat
from ruse.rdkit.mcs import mcs_with_constant_smiles, mcs_match_fragments
from ruse.rdkit.rdkit_utils import encode_mol, rescale_bond_lengths
from ruse.util.data_table import DataTable


class MmpTransformProduct(NamedTuple):
    """
    Base class for transform identified by mmpdb

    Attributes:

        - from_smiles: smiles with mapped wildcard atoms encoding variable part of query

        - to_smiles: smiles with mapped wildcard atoms encoding variable part of product

        - query:  :class:`rdkit.Chem.rdchem.Mol` molecule containing query molecule.  Annotated with bonds broken and variable core

        - product: :class:`rdkit.Chem.rdchem.Mol` molecule containing product molecule.  Annotated with bonds broken and variable core

        - reversed: set true if the transformation is reversed compared to the database rule

        - value: mean numeric difference observed for the transform

        - n_pairs: number of pair values that contribute to the mean difference

        - radius: radius of match between query and pairs in rule

        - rule_environment_id: mmpdb database id for rule

        - constant_smiles: smiles for conserved portion of rule

    """
    from_smiles: str
    to_smiles: str
    query: Mol
    product: Mol
    reversed: bool
    value: float
    n_pairs: int
    radius: int
    rule_environment_id: int
    constant_smiles: str
    query_reference_missing: List[int]


class RadiusAndRuleGroup(NamedTuple):
    rule_env_id: int
    radius: int
    query: Mol
    product: Mol
    from_smiles: str
    to_smiles: str
    constant_smiles: str
    properties: List[float]
    query_reference_missing: List[int]

    def to_grid_row(self, molecules: bool = False, include_constant_smiles: bool = False,
                    build_fingerprints: bool = False, reference_query: Mol = None) -> List:
        if molecules:
            values = [self.query, self.product, Chem.MolFromSmiles(self.from_smiles),
                      Chem.MolFromSmiles(self.to_smiles), self.radius, self.rule_env_id]
        else:
            values = [encode_mol(RDKitFormat.sdf, self.query), encode_mol(RDKitFormat.sdf, self.product),
                      self.from_smiles, self.to_smiles, self.radius, self.rule_env_id]
        values = values + self.properties
        if build_fingerprints and reference_query:
            values.append(build_fingerprint(reference_query, self.query, 'variable'))
        if include_constant_smiles:
            values.append(self.constant_smiles)
        if build_fingerprints and reference_query:
            values.append(build_fingerprint(reference_query, self.query))
        return values


class MmpProduct(NamedTuple):
    """
    Class to store a query and product and associated transforms.  There will be one transform for each property
    requested

    Attributes:

        - query:  query molecule of :class:`rdkit.Chem.rdchem.Mol`

        - product: product molecule :class:`rdkit.Chem.rdchem.Mol`

        - transform_products: A list of :class:`MmpTransformProduct` containing query to product transforms

    """
    query: Mol
    product: Mol
    transform_products: List[MmpTransformProduct]

    def group_by_rule_and_radius(self) -> List[RadiusAndRuleGroup]:
        pairs = list({(-tp.rule_environment_id if tp.reversed else tp.rule_environment_id, tp.radius)
                      for tp in self.transform_products if tp})
        groups = []
        for rule_environment_id, radius in pairs:
            transform_product = None
            properties = []
            for tp in self.transform_products:
                if not tp:
                    properties.append(None)
                    continue
                rule_env_id = -tp.rule_environment_id if tp.reversed else tp.rule_environment_id
                if rule_env_id == rule_environment_id and tp.radius == radius:
                    if not transform_product:
                        transform_product = tp
                    else:
                        assert tp.from_smiles == transform_product.from_smiles
                        assert tp.to_smiles == transform_product.to_smiles
                    properties.append(tp.value)
                else:
                    properties.append(None)
            group = RadiusAndRuleGroup(rule_environment_id, radius, transform_product.query, transform_product.product,
                                       transform_product.from_smiles, transform_product.to_smiles,
                                       transform_product.constant_smiles, properties,
                                       transform_product.query_reference_missing)
            groups.append(group)
        return groups


class MmpTransform(NamedTuple):
    """
    Class to store all transforms found for a query

    Attributes:

        - query: query molecule of :class:`rdkit.Chem.rdchem.Mol`

        - property_names: A list of the property names fior which MMP delta are requested

        - property_ids: MMPDB ids for the property names

        - products: List of products for the query as a list of :class:`MmmpProduct` objects

    """

    query: Mol
    property_names: List[str]
    property_ids: List[int]
    products: List[MmpProduct]

    def to_data_table(self, input_format: RDKitFormat, input_query: str,
                      groups: Optional[List[RadiusAndRuleGroup]] = None) -> DataTable:
        """
        Creates a ruse data table (:class:` ruse.util.data_table.DataTable`), with columns for query, product,
        property deltas, rule id and environment radius

        :return:
        """

        if input_format == RDKitFormat.smi:
            query_property = input_query
        elif input_format == RDKitFormat.sdf:
            mol_gzip = gzip.compress(bytes(input_query, 'utf-8'))
            query_property = base64.b64encode(mol_gzip).decode('utf-8')
        else:
            raise ValueError('Unknown  mmp query type {}'.format(input_format.name))

        columns = [DataTable.column_definition('Query', 'binary', 'chemical/x-mdl-molfile',
                                               properties={'mmpInputQuery': query_property}),
                   DataTable.column_definition('Product', 'binary', 'chemical/x-mdl-molfile'),
                   DataTable.column_definition('FROM', 'string', 'chemical/x-smiles'),
                   DataTable.column_definition('TO', 'string', 'chemical/x-smiles'),
                   DataTable.column_definition('RADIUS', 'int'),
                   DataTable.column_definition('ENV_RULE', 'int')]

        for p in self.property_names:
            property_column = DataTable.column_definition(p, 'float', properties={'mmpPropertyColumn': 'true'})
            columns.append(property_column)

        columns.append(
            DataTable.column_definition('TransformFingerprint', 'binary', 'application/octet-stream'))
        columns.append(DataTable.column_definition('Fingerprint', 'binary', 'application/octet-stream'))

        rows = self.to_grid(build_fingerprints=True, groups=groups)
        data_table = DataTable(columns=columns, data=rows)
        return data_table

    def to_split_files(self, input_format: RDKitFormat, input_query: str):

        file_no = 0
        atom_ids_to_file = {}
        file_groups = {}
        rows = []

        for p in self.products:
            for group in p.group_by_rule_and_radius():
                query_ids = tuple(group.query_reference_missing)
                if query_ids in atom_ids_to_file:
                    file = atom_ids_to_file[query_ids]
                else:
                    file = "mmp_data_table_{}.json".format(file_no)
                    file_no += 1
                    atom_ids_to_file[query_ids] = file
                fingerprint = build_fingerprint_from_atoms(self.query, group.query_reference_missing)
                group_list = file_groups.setdefault(file, [])
                group_no = len(group_list)
                group_list.append(group)
                rows.append([file, group_no, fingerprint])

        file_data_tables = {file: self.to_data_table(input_format, input_query, groups=groups) for file, groups in
                            file_groups.items()}
        for file, data_table in file_data_tables.items():
            with open(file, 'w', encoding='utf8') as fh:
                json.dump(data_table.to_data(), fh)

        columns = [DataTable.column_definition('File', 'string'),
                   DataTable.column_definition('Row', 'int'),
                   DataTable.column_definition('Fingerprint', 'binary', 'application/octet-stream')]
        data_table = DataTable(columns=columns, data=rows)
        return data_table

    def to_grid(self, molecules: bool = False, column_major: bool = False, include_constant_smiles: bool = False,
                build_fingerprints: bool = False, groups: Optional[List[RadiusAndRuleGroup]] = None) \
            -> List[List[Union[str, int, float, Mol]]]:

        """
        Converts the query and products to a matrix data structure, with columns for query, product,
        property deltas, rule id and environment radius.  This can be used to build ruse or pandas
        data tables

        :param molecules: Set true to use RDKit molecules for cells containing structures
        :param column_major: Set true for column major indexing (default is row major)
        :return: matrix as list of lists
        """

        if not groups:
            groups = [g for p in self.products for g in p.group_by_rule_and_radius()]
        rows = [g.to_grid_row(molecules, include_constant_smiles, build_fingerprints, self.query) for g in groups]

        return rows if not column_major else [list(c) for c in zip(*rows)]


def map_labelled_atoms_to_reference(reference: Mol, transform: Mol, atom_property: str = 'missing') -> \
        List[int]:
    mapping = transform.GetSubstructMatch(reference)
    assert mapping
    transform_to_query = {i: m for m, i in enumerate(mapping)}
    labelled_transform_atoms = [atom.GetIdx() for atom in transform.GetAtoms()
                                if atom.HasProp(atom_property) and atom.GetBoolProp(atom_property)]
    labelled_reference_atoms = [transform_to_query[t] for t in labelled_transform_atoms]
    for input_idx, transform_idx in zip(labelled_reference_atoms, labelled_transform_atoms):
        assert reference.GetAtomWithIdx(input_idx).GetAtomicNum() == \
               transform.GetAtomWithIdx(transform_idx).GetAtomicNum()
    labelled_reference_atoms.sort()
    return labelled_reference_atoms


def build_fingerprint_from_atoms(reference_query: Mol, labelled_reference_atoms: List[int]) -> str:
    changed_input_bonds = [bond.GetIdx() for bond in reference_query.GetBonds()
                           if bond.GetBeginAtomIdx() in labelled_reference_atoms
                           and bond.GetEndAtomIdx() in labelled_reference_atoms]
    data = array('b')
    for i in range(reference_query.GetNumAtoms()):
        val = 0x1 if i in labelled_reference_atoms else 0
        data.append(val)
    for i in range(reference_query.GetNumBonds()):
        val = 0x1 if i in changed_input_bonds else 0
        data.append(val)
    compress_data = gzip.compress(data)
    encoded_data = base64.b64encode(compress_data).decode('utf-8')
    return encoded_data


def atoms_and_bonds_from_fingerprint(reference_query: Mol, fingerprint: str) -> Tuple[List[int], List[int]]:
    compress_data = base64.b64decode(fingerprint.encode('utf-8'))
    data = gzip.decompress(compress_data)
    n_atoms = reference_query.GetNumAtoms()
    assert (len(data) == n_atoms + reference_query.GetNumBonds())
    atoms = [i for i, v in enumerate(data) if i < n_atoms and v == 0x1]
    bonds = [i - n_atoms for i, v in enumerate(data) if i >= n_atoms and v == 0x1]
    return atoms, bonds


def build_fingerprint(reference_query: Mol, transform_query: Mol, atom_property: str = 'missing',
                      ) -> str:
    labelled_reference_atoms = map_labelled_atoms_to_reference(reference_query, transform_query, atom_property)
    return build_fingerprint_from_atoms(reference_query, labelled_reference_atoms)


class MmpPair(NamedTuple):
    """
    A pair from the MMPDB database

    Attributes:

        - from_id: database id of LHS of pair

        - from_smiles: structure of LHS

        - from_property: LHS property value

        - to_id: database id of RHS of pair

        - to_smiles: structure of RHS

        - to_property: RHS property value
    """
    from_id: str
    from_smiles: str
    from_property: float
    to_id: str
    to_smiles: str
    to_property: float


class MmpEnvironment(NamedTuple):
    """
    Base class to define the pairs for a property transform.  This class is only used to validate Ruse webservice

    Attributes:

        - pairs: a list of pairs used to define a transform
    """
    pairs: List[MmpPair]

    def delta(self):
        """
        Determines overall properta delta from the list of pairs

        :return: property difference from transform
        """
        diffs = [p.to_property - p.from_property for p in self.pairs]
        mean = sum(diffs) / float(len(diffs))
        return mean

    def to_grid(self, molecules: bool = False):
        rows = []
        for p in self.pairs:
            if molecules:
                row = [p.from_id, Chem.MolFromSmiles(p.from_smiles), p.from_property,
                       p.to_id, Chem.MolFromSmiles(p.to_smiles), p.to_property]
            else:
                row = [p.from_id, p.from_smiles, p.from_property,
                       p.to_id, p.to_smiles, p.to_property]
            rows.append(row)
        return rows


def transform(smiles: str, properties: List[str], database_file: str) -> TransformResult:
    """
    A rewrite of :func:`mmpdb.do_analysis.transform_command` to return in-memory transform results rather than
    print them out.

    :param smiles: input structure
    :param properties: properties to estimate effect on
    :param database_file: MMP database
    :return: result of the transform
    """

    # MMPDB 2 stuff:
    try:
        args_in = ['transform', '--smiles', smiles, database_file]
        for p in properties:
            args_in.extend(['--property', p])
        args, _ = mmp.parser.parse_known_args(args_in)
        parser = args.subparser
        min_radius = args.min_radius
        assert min_radius in list("012345"), min_radius
        min_radius = int(min_radius)
        min_pairs = int(args.min_pairs)
        min_variable_size = args.min_variable_size
        min_constant_size = args.min_constant_size
        explain = command_support.get_explain(args.explain)
        dataset = dbutils.open_dataset_from_args_or_exit(args)
        property_names = command_support.get_property_names_or_error(parser, args, dataset)

        if args.substructure:
            substructure_pat = Chem.MolFromSmarts(args.substructure)
            if substructure_pat is None:
                parser.error("Cannot parse --substructure %r" % (args.substructure,))
        else:
            substructure_pat = None

        # evaluate --where, --score, and --rule-selection-cutoffs.
        rule_selection_function = analysis_algorithms.get_rule_selection_function_from_args(
            parser, args)

        transform_tool = analysis_algorithms.get_transform_tool(dataset, rule_selection_function)
        transform_record = transform_tool.fragment_transform_smiles(args.smiles)
        if transform_record.errmsg:
            parser.error("Unable to fragment --smiles %r: %s"
                         % (args.smiles, transform_record.errmsg))

        if args.jobs > 1:
            pool = multiprocessing.Pool(processes=args.jobs)
        else:
            pool = None
        try:
            result = transform_tool.transform(
                transform_record.fragments, property_names,
                min_radius=min_radius,
                min_pairs=min_pairs,
                min_variable_size=min_variable_size,
                min_constant_size=min_constant_size,
                substructure_pat=substructure_pat,
                pool=pool,
                explain=explain,
            )
        except analysis_algorithms.EvalError as err:
            sys.stderr.write("ERROR: %s\nExiting.\n" % (err,))
            raise SystemExit(1)

        return result
    except AttributeError:
        # MMPDB3
        db = dbutils.open_database(database_file)
        dataset = db.get_dataset()
        transform_tool = get_transform_tool(dataset)
        transform_record = transform_tool.fragment_transform_smiles(smiles)
        result = transform_tool.transform(transform_record.fragmentations, properties)
        return result

def create_combined_structure(constant_smiles: str, variable_smiles: str,
                              attachment_order: Optional[List[int]] = None) -> Mol:
    """
    From constant and variable substructures create combined structure with labelled atoms indicating bond cuts
    and renderer_highlight molecular property to show variable substructure.

    :param constant_smiles: constant portion of molecule
    :param variable_smiles:  variable portion of molecule
    :param attachment_order: attachment order of variable atoms.  If None, the attachment points in the variable_smiles should be numbered

    :return: Combined structure as :class:`rdkit.Chem.rdchem.Mol`
    """
    constant_mol = Chem.MolFromSmiles(constant_smiles)
    variable_mol = Chem.MolFromSmiles(variable_smiles)

    variable_attachment_no = 0
    for variableNo, atom in enumerate(variable_mol.GetAtoms()):
        atom.SetBoolProp('variable', True)
        atom.SetIntProp('variableNo', variableNo)
        if atom.GetAtomicNum() == 0:
            atom.SetProp("attachment", "variable")
            if atom.HasProp('molAtomMapNumber'):
                variable_index = atom.GetIntProp('molAtomMapNumber') - 1
            else:
                variable_index = variable_attachment_no
            if attachment_order:
                atom.SetIntProp("attachmentIdx", attachment_order[variable_index])
            else:
                atom.SetIntProp("attachmentIdx", variable_index)
            variable_attachment_no += 1

    constant_attachment_no = 0
    for atom in constant_mol.GetAtoms():
        atom.SetBoolProp('variable', False)
        if atom.GetAtomicNum() == 0:
            atom.SetProp("attachment", "constant")
            atom.SetIntProp("attachmentIdx", constant_attachment_no)
            constant_attachment_no += 1

    if attachment_order:
        assert variable_attachment_no == len(attachment_order)
        assert constant_attachment_no == len(attachment_order)

    combine_mol = Chem.CombineMols(constant_mol, variable_mol)
    Chem.SanitizeMol(combine_mol)
    mol = Chem.RWMol(combine_mol)
    n_cuts = constant_attachment_no
    mol.SetIntProp("nCuts", n_cuts)

    for index in range(0, variable_attachment_no):
        _join_mols(mol, index)

    Chem.SanitizeMol(mol)
    Chem.AssignStereochemistry(mol)

    def highlight(atom: Atom) -> bool:
        return atom.HasProp('variable') and atom.GetBoolProp('variable')

    def highlight_bond(bond: Bond, highlight_atoms: List[int]) -> bool:
        start = bond.GetBeginAtomIdx()
        end = bond.GetEndAtomIdx()
        return start in highlight_atoms and end in highlight_atoms

    cuts_bond_idxs = [_bond_idx_for_cut(mol, cut_no) for cut_no in range(1, n_cuts + 1)]

    mol.SetProp("bond_cuts", ' '.join([str(b + 1) for b in cuts_bond_idxs]))
    highlight_atoms = [a.GetIdx() for a in mol.GetAtoms() if highlight(a)]
    highlight_bonds = [b.GetIdx() for b in mol.GetBonds() if highlight_bond(b, highlight_atoms)]

    prop_text = "COLOR #ff0000\nATOMS {}\nBONDS {}".format(' '.join([str(a + 1) for a in highlight_atoms]),
                                                           ' '.join([str(b + 1) for b in highlight_bonds]))
    mol.SetProp('Renderer_Highlight', prop_text)

    return Mol(mol)


def _parity_shell(values):
    """
    Determines the parity of integers in a list

    :param values:
    :return: parity
    """

    # from http://www.dalkescientific.com/writings/diary/archive/2016/08/15/fragment_parity_calculation.html
    # Simple Shell sort; while O(N^2), we only deal with at most 4 values
    values = list(values)
    N = len(values)
    num_swaps = 0
    for i in range(N - 1):
        for j in range(i + 1, N):
            if values[i] > values[j]:
                values[i], values[j] = values[j], values[i]
                num_swaps += 1
    return num_swaps % 2


def _atom_ordering_before(atom: Atom, old_atom: Atom, new_atom: Atom) -> List[Atom]:
    """
    Determine the order of neighbor atoms in the bond table, substituting new_atom for old_atom

    :param atom: center atom
    :param old_atom:
    :param new_atom:
    :return: Ordering
    """
    env_before = [b.GetOtherAtom(atom) for b in atom.GetBonds()]
    index, _ = next((i, a) for i, a in enumerate(env_before) if a.GetIdx() == old_atom.GetIdx())
    env_before[index] = new_atom
    return env_before


def _check_chirality_after(atom: Atom, env_before: List[Atom]) -> None:
    """
    Determine if a chiral inversion has occurred round an atom after creating combined structure

    :param atom:
    :param env_before: Ordering prior to merging
    """
    ids_after = [b.GetOtherAtomIdx(atom.GetIdx()) for b in atom.GetBonds()]
    ids_before = [a.GetIdx() for a in env_before]
    parity_before = _parity_shell(ids_before)
    parity_after = _parity_shell(ids_after)

    if parity_before != parity_after:
        atom.InvertChirality()


def _set_mapping_label(atom: Atom, mapping_no: int) -> None:
    """
    Add labels to atoms to indicate a bond the variable substructure
    :param atom:
    :param mapping_no:
    :return:
    """
    atom.SetBoolProp('cutNo{}'.format(mapping_no), True)
    if atom.HasProp('molAtomMapNumber'):
        mapping_no = int('{}{}'.format(atom.GetIntProp('molAtomMapNumber'), mapping_no))
    atom.SetAtomMapNum(mapping_no)
    atom.SetIntProp('molAtomMapNumber', mapping_no)


def _find_neighbour_stereo_bond(mol: Mol, atom: Atom, old_atom: Atom, new_atom: Atom) -> Optional[
    Tuple[Bond, List[Atom]]]:
    """
    Determines if any trans/cis stereo bonds are affected by the fragment joining.  Returns the affected bond and defining atoms

    :param mol: molecule
    :param atom: center atom
    :param old_atom: old neighbor atom
    :param new_atom: new neighbor atom
    :return: Tuple of any affected bond and atoms, or None
    """

    for bond in atom.GetBonds():
        stereo = bond.GetStereo()
        if stereo != BondStereo.STEREONONE:
            stereo_atoms = [mol.GetAtomWithIdx(i) for i in bond.GetStereoAtoms()]
            try:
                index, _ = next((i, a) for i, a in enumerate(stereo_atoms) if a.GetIdx() == old_atom.GetIdx())
                stereo_atoms[index] = new_atom
            except StopIteration:
                # note need to reset stereo bond information even if deleted atom is not in the list of stereo atoms
                pass
            return bond, stereo_atoms
    return None


def _join_mols(mol: RWMol, attachment_no: int) -> None:
    """
    Joins two disconnected fragment in the molecule by connecting wildcard atoms for the attachment points

    :param mol:
    :param attachment_no:
    """
    constant_atom = None
    variable_atom = None
    constant_neighbour = None
    variable_neighbour = None

    for atom in mol.GetAtoms():

        if atom.GetAtomicNum() == 0:

            if atom.HasProp("attachment") and atom.GetProp("attachment") == "constant" and atom.GetIntProp(
                    "attachmentIdx") == attachment_no:
                constant_atom = atom
                neighbours = atom.GetNeighbors()
                assert len(neighbours) == 1
                constant_neighbour = neighbours[0]

            if atom.HasProp("attachment") and atom.GetProp("attachment") == "variable" and atom.GetIntProp(
                    "attachmentIdx") == attachment_no:
                variable_atom = atom
                neighbours = atom.GetNeighbors()
                assert len(neighbours) == 1
                variable_neighbour = neighbours[0]

    assert constant_atom
    assert variable_atom
    assert constant_neighbour
    assert variable_neighbour

    variable_env_before = _atom_ordering_before(variable_neighbour, variable_atom, constant_neighbour)
    constant_env_before = _atom_ordering_before(constant_neighbour, constant_atom, variable_neighbour)

    stereo_bonds = {}
    t = _find_neighbour_stereo_bond(mol, variable_neighbour, variable_atom, constant_neighbour)
    if t:
        stereo_bond, stereo_atoms = t
        stereo_bonds[stereo_bond] = stereo_atoms
    t = _find_neighbour_stereo_bond(mol, constant_neighbour, constant_atom, variable_neighbour)
    if t:
        stereo_bond, stereo_atoms = t
        stereo_bonds[stereo_bond] = stereo_atoms

    mol.RemoveAtom(constant_atom.GetIdx())
    mol.RemoveAtom(variable_atom.GetIdx())
    mol.AddBond(constant_neighbour.GetIdx(), variable_neighbour.GetIdx(), BondType.SINGLE)

    # see Andrew Dalke's blog for how to handle chiral centers on fragmentation
    #  http://www.dalkescientific.com/writings/diary/archive/2016/08/14/fragment_chiral_molecules.html#fragment_chiral
    _check_chirality_after(variable_neighbour, variable_env_before)
    _check_chirality_after(constant_neighbour, constant_env_before)

    for stereo_bond in stereo_bonds:
        atom_indices = [a.GetIdx() for a in stereo_bonds[stereo_bond]]
        assert len(atom_indices) == 2
        stereo_bond.SetStereoAtoms(atom_indices[0], atom_indices[1])

    _set_mapping_label(variable_neighbour, attachment_no + 1)
    _set_mapping_label(constant_neighbour, attachment_no + 1)


def align_combined_molecules(query: Mol, mapped_query_mol: Mol, mapped_product_mol: Mol,
                             constant_smiles: str) -> None:
    """
    Align the combined query and product molecules to the initial query coordinates for the largest fragment in the
    constant smiles

    :param query:
    :param mapped_query_mol:
    :param mapped_product_mol:
    :param constant_smiles:
    """
    if mapped_query_mol.HasSubstructMatch(query, useChirality=True):
        AllChem.GenerateDepictionMatching2DStructure(mapped_query_mol, query)
    else:
        AllChem.Compute2DCoords(mapped_query_mol)

    template_smiles = constant_smiles.split('.')[0]
    template_mol = Chem.MolFromSmarts(template_smiles)
    map = mapped_query_mol.GetSubstructMatch(template_mol, useChirality=True)
    if mapped_product_mol.HasSubstructMatch(template_mol, useChirality=True) \
            and map:
        template_conformer = Conformer(template_mol.GetNumAtoms())
        query_conformer = mapped_query_mol.GetConformer(0)
        for template_no, query_no in enumerate(map):
            query_point = query_conformer.GetAtomPosition(query_no)
            # print('x {} y {} z {}'.format(query_point.x, query_point.y, query_point.z))
            template_conformer.SetAtomPosition(template_no, Point3D(query_point.x, query_point.y, query_point.z))
        template_mol.AddConformer(template_conformer)
        AllChem.GenerateDepictionMatching2DStructure(mapped_product_mol, template_mol)
    else:
        AllChem.Compute2DCoords(mapped_product_mol)


def result_to_mmp_transform(query: Mol, result: TransformResult) -> MmpTransform:
    """
    Map the data structures returned from an MMPDB query search to our own data structures.

    :param query:
    :param result:
    :return:
    """

    if not query.GetNumConformers():
        AllChem.Compute2DCoords(query)
    else:
        rescale_bond_lengths(query)

    property_ids, properties = zip(*result.property_info_list)
    products = []
    last_from_smiles = None
    last_to_smiles = None
    last_constant_smiles = None

    for product_no, tp in enumerate(result.transform_products):
        product_smiles = tp.smiles
        transform_products = []
        has_product_with_rule = False
        for rule in tp.property_rules:

            if not rule or not rule.property_rule:
                transform_products.append(None)
                continue

            has_product_with_rule = True
            reversed = True if rule.property_rule and rule.is_reversed == 1 else False
            from_smiles = add_attachment_order_to_smiles(rule.variable_smiles)
            # for multiple properties the previous product will often be the same as this one.  If it is reuse it
            previous_matches = last_from_smiles == from_smiles and last_to_smiles == rule.to_smiles \
                               and last_constant_smiles == rule.constant_smiles

            if not previous_matches:

                # create combined structure with labeled variable atoms and attachment points
                mapped_query_mol = create_combined_structure(rule.constant_smiles, from_smiles)
                mapped_product_mol = create_combined_structure(rule.constant_smiles, rule.to_smiles)

                # align around constant core for depiction
                align_combined_molecules(query, mapped_query_mol, mapped_product_mol, rule.constant_smiles)

                # identify a reasonable minimal difference between query and product. This is usually some subset of
                # the from_smiles pattern. We can look for differences in just the transform or in the whole molecule.
                # whole molecule MCSs are better, but slower. There are some edge cases where it is not possible to get
                # a reasonable difference without looking at the whole molecule.
                whole_molecule_mcs = False
                if whole_molecule_mcs:
                    mapping = mcs_with_constant_smiles(mapped_query_mol, mapped_product_mol, rule.constant_smiles)
                    query_mapping, _ = list(zip(*mapping))
                    query_missing = [a for a in range(mapped_query_mol.GetNumAtoms()) if a not in query_mapping]
                else:
                    from_mol = Chem.MolFromSmiles(from_smiles)
                    to_mol = Chem.MolFromSmiles(rule.to_smiles)
                    fragment_mapping = mcs_match_fragments(from_mol, to_mol)
                    query_mapping, _ = list(zip(*fragment_mapping)) if fragment_mapping else [(), ()]
                    query_missing = [atom.GetIdx() for atom in mapped_query_mol.GetAtoms()
                                     if
                                     atom.HasProp('variableNo') and atom.GetIntProp(
                                         'variableNo') not in query_mapping]
                    for atom in from_mol.GetAtoms():
                        if atom.GetAtomicNum() == 0 and atom.GetIdx() not in query_mapping:
                            cut_no = atom.GetAtomMapNum()
                            wildcard_match = _atom_matching_wildcard(mapped_query_mol, cut_no)
                            query_missing.append(wildcard_match)

                for missing_idx in query_missing:
                    mapped_query_mol.GetAtomWithIdx(missing_idx).SetBoolProp('missing', True)

                # find missing atoms in original query
                reference_missing = map_labelled_atoms_to_reference(query, mapped_query_mol)

                # this removes annoying substructure highlighting in IPythonConsole
                #  (which monkey patches substructure search)
                if hasattr(mapped_product_mol, '__sssAtoms'):
                    delattr(mapped_product_mol, '__sssAtoms')
                if hasattr(mapped_query_mol, '__sssAtoms'):
                    setattr(mapped_query_mol, '__sssAtoms', query_missing)

                last_from_smiles = from_smiles
                last_constant_smiles = rule.constant_smiles
                last_to_smiles = rule.to_smiles

            transform_product = MmpTransformProduct(from_smiles=from_smiles, to_smiles=rule.to_smiles,
                                                    n_pairs=rule.count, reversed=reversed, radius=rule.radius,
                                                    value=rule.avg, rule_environment_id=rule.rule_environment_id,
                                                    query=mapped_query_mol,
                                                    product=mapped_product_mol,
                                                    constant_smiles=rule.constant_smiles,
                                                    query_reference_missing=reference_missing)
            transform_products.append(transform_product)

        if has_product_with_rule:
            product = MmpProduct(query=query, product=product_smiles, transform_products=transform_products)
            products.append(product)
    return MmpTransform(query=query, property_names=properties, products=products, property_ids=property_ids)


def environment_rule_pairs(rule_environment_id: int, reverse: bool, property_name, database_file: str) \
        -> MmpEnvironment:
    """
    Recover the database of molecule pairs that are used for an environment rule

    :param rule_environment_id:
    :param reverse:
    :param property_name:
    :param database_file:
    :return:
    """
    db = open_database(database_file)
    cursor = db.get_cursor()

    query = 'select compound1_id, compound2_id from pair where rule_environment_id = ?'
    cursor.execute(query, [rule_environment_id])
    pairs = []
    for compound1_id, compound2_id in cursor.fetchall():
        if reverse:
            compound1_id, compound2_id = compound2_id, compound1_id

        query = """ select public_id, clean_smiles, value 
                      from compound c, compound_property cp, property_name pn
                     where c.id = ?
                       and pn.name = ?
                       and c.id = cp.compound_id 
                       and cp.property_name_id = pn.id"""
        cursor.execute(query, (compound1_id, property_name))
        (compound1_name, compound1_smiles, compound1_value) = cursor.fetchone()
        cursor.execute(query, (compound2_id, property_name))
        (compound2_name, compound2_smiles, compound2_value) = cursor.fetchone()
        pair = MmpPair(from_id=compound1_name, from_smiles=compound1_smiles, from_property=compound1_value,
                       to_id=compound2_name, to_smiles=compound2_smiles, to_property=compound2_value)
        pairs.append(pair)

    return MmpEnvironment(pairs=pairs)


def _atom_has_cut_no(atom: Atom, cut_no: int) -> bool:
    """
    Return true if this atom is labelled with the cut point

    :param atom:
    :param cut_no:
    :return:
    """
    prop = 'cutNo{}'.format(cut_no)
    if atom.HasProp(prop) and atom.GetBoolProp(prop):
        return True
    return False


def _bond_idx_for_cut(mol: Mol, cut_no: int) -> int:
    """
    Find the bond index for a given cut

    :param mol:
    :param cut_no:
    :return:
    """
    for atom in mol.GetAtoms():
        if _atom_has_cut_no(atom, cut_no):
            for neighbor in atom.GetNeighbors():
                if _atom_has_cut_no(neighbor, cut_no):
                    bond = mol.GetBondBetweenAtoms(atom.GetIdx(), neighbor.GetIdx())
                    return bond.GetIdx()
            raise RuntimeError('Unable to find neighbor atom in cut')
    raise ValueError('Unable to find cut {} in molecule'.format(cut_no))


def _atom_matching_wildcard(mol: Mol, cut_no: int) -> int:
    """
    Find the index of the atom in a combined molecule that matches the transform dummy atom for the given
    cut (this the the atom in the constant part of the molecule that matches the cut)

    :param mol:
    :param cut_no:
    :return:
    """
    return next(atom.GetIdx() for atom in mol.GetAtoms()
                if _atom_has_cut_no(atom, cut_no) and not atom.GetBoolProp('variable'))


def add_attachment_order_to_smiles(smiles: str, attachment_order: Optional[List[int]] = None) -> str:
    """
    Adds the specified mapping numbers to the transform smiles

    :param smiles:
    :param attachment_order:
    :return:
    """
    mol = Chem.MolFromSmiles(smiles)
    atoms = [a for a in mol.GetAtoms() if a.GetAtomicNum() == 0]
    if not attachment_order:
        attachment_order = range(len(atoms))
    assert (len(atoms) == len(attachment_order))
    for attachment_no, atom in zip(attachment_order, atoms):
        atom.SetIntProp('molAtomMapNumber', attachment_no + 1)
    return Chem.MolToSmiles(mol, True)


def database_property_names(database_file: str) -> List[str]:
    """
    Returns all the property names for a MMPDB database

    :param database_file:
    :return:
    """
    dbinfo = DBFile(database_file)
    database = dbinfo.open_database()
    dataset = database.get_dataset()
    return dataset.get_property_names()
