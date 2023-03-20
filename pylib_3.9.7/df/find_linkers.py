#!/usr/bin/env python

import argparse
import concurrent.futures as cf
import csv
import gzip
import json
import random
import sqlite3
import sys

from collections import defaultdict
from itertools import repeat
from pathlib import Path
from psutil import cpu_count
from typing import IO, Optional, Union

from rdkit import rdBase, Chem
from rdkit.Chem import rdmolops, rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D

from float_range import FloatRange
from int_range import IntRange


# Find linkers in a file of structures using method described in
# https://www.sciencedirect.com/science/article/abs/pii/S0968089623000421#b0255


class Linker:
    """
    Class for holding a linker, which is a short piece of molecule
    between 2 rings.  Includes SMILES for the left and right hand
    pieces off the linker and the whole of the parent molecule.
    """

    def __init__(self, name: str, linker_smiles: str, mol_smiles: str,
                 left_smiles: str, right_smiles: str,
                 num_atoms: int, linker_length: int,
                 num_donors: int, num_acceptors: int):
        self.name = name
        self.linker_smiles = linker_smiles
        self._num_atoms = num_atoms
        self.mol_smiles = mol_smiles
        self.left_smiles = left_smiles
        self.right_smiles = right_smiles

        self._reversed_linker = None
        self._reversed_left_smiles = None
        self._reversed_right_smiles = None

        self._symmetrical = None
        self._path_length = linker_length
        self._num_donors = num_donors
        self._num_acceptors = num_acceptors

    def __str__(self) -> str:
        retval = (f'Linker {self.name} SMILES : {self.linker_smiles} LeftSMILES'
                  f' : {self.left_smiles} RightSMILES : {self.right_smiles}'
                  f' Symmetrical : {self.symmetrical}')
        return retval

    def __eq__(self, other) -> bool:
        if not isinstance(other, Linker):
            return NotImplemented
        if self.linker_smiles != other.linker_smiles:
            return False

        if (self.left_smiles == other.left_smiles
                and self.right_smiles == other.right_smiles):
            return True
        if self.symmetrical:
            olsmi, orsmi = self.reversed_sides
            if olsmi == other.left_smiles and orsmi == other.right_smiles:
                return True
        return False

    @property
    def path_length(self) -> int:
        return self._path_length

    @property
    def num_atoms(self) -> int:
        return self._num_atoms

    @property
    def symmetrical(self) -> bool:
        if self._symmetrical is None:
            self._assign_symmetrical()
        return self._symmetrical

    @property
    def num_donors(self) -> int:
        return self._num_donors

    @property
    def num_acceptors(self) -> int:
        return self._num_acceptors

    @property
    def reversed_sides(self) -> tuple[str, str]:
        """
        Make SMILES strings of the reversed sides.
        Swap left and right for linker1 and compare to linker2.  It's
        almost always the case that changing the map number doesn't
        change the canonical SMILES (apart from the obvious), but
        from time to time it does e.g.
        O=C1C2CCCCC2C(=O)N1CC(=O)N1CCCSC1=[*:2] goes to
        O=C1C2CCCCC2C(=O)N1CC(=O)N1CCCSC1=[*:1] canonicalises to
        O=C(CN1C(=O)C2CCCCC2C1=O)N1CCCSC1=[*:1]
        """
        if self._reversed_left_smiles is None:
            orsmi = self.left_smiles.replace('[*:1]', '[*:2]')
            self._reversed_right_smiles = Chem.MolToSmiles(Chem.MolFromSmiles(orsmi))
            olsmi = self.right_smiles.replace('[*:2]', '[*:1]')
            self._reversed_left_smiles = Chem.MolToSmiles(Chem.MolFromSmiles(olsmi))

        return self._reversed_left_smiles, self._reversed_right_smiles

    @property
    def reversed_linker(self) -> str:
        """
        Make SMILES string of reversed linker i.e. [*:1] becomes [*:2]
        and vice versa.
        """
        if self._reversed_linker is None:
            other_smi = self.linker_smiles.replace('[*:1]', '[*:XX]')
            other_smi = other_smi.replace('[*:2]', '[*:1]')
            other_smi = other_smi.replace('[*:XX]', '[*:2]')
            self._reversed_linker = Chem.MolToSmiles(Chem.MolFromSmiles(other_smi))
        return self._reversed_linker

    def _assign_symmetrical(self):
        """
        Returns True if linker is symmetrical, which is decided by making a
        new mol with the :1 and :2 atom maps reversed and seeing if it has
        the same canonical SMILES.  Assumes linker_smi is already
        canonical.
        """
        other_smi = self.reversed_linker
        check_smi = Chem.MolToSmiles(Chem.MolFromSmiles(other_smi))
        self._symmetrical = self.linker_smiles == check_smi


def create_mol_supplier(infile) -> Union[Chem.ForwardSDMolSupplier, Chem.SmilesMolSupplier, None]:
    """

    Args:
        infile (str): must end .smi, .sdf or .sdf.gz

    Returns:
        ForwardSDMolSupplier or None
    """
    inpath = Path(infile)
    sfx = inpath.suffix
    gzipped = False
    if sfx == '.gz':
        suffixes = inpath.suffixes
        gzipped = True
        sfx = suffixes[-2]

    if sfx != '.smi' and sfx != '.sdf' and sfx != '.mol':
        print(f'ERROR : input must be a SMILES, SDF, MOL or gzipped SDF'
              f' or MOL, you gave {infile}  Aborting.')
        return None

    if sfx == '.smi':
        return Chem.SmilesMolSupplier(infile, titleLine=False)
    else:
        try:
            if gzipped:
                inf = gzip.open(infile)
                return Chem.ForwardSDMolSupplier(inf)
            else:
                return Chem.ForwardSDMolSupplier(infile)
        except (OSError, FileNotFoundError):
            print(f'ERROR : failed to open file {infile}.  Not Found.')
            return None


def parse_args():
    """
        Use argparse module to parse the command line arguments.

        Returns:
            (Namespace object): Parsed arguments in argparse Namespace
                                object.
        """
    parser = argparse.ArgumentParser(description='Find linkers in molecules.')
    parser.add_argument('-I', '--input-file', dest='infile',
                        required=True,
                        help='Name of file containing conformations.  Must be'
                             ' SMILES, SDF or gzipped SDF.')
    parser.add_argument('-R', '--random-fraction', dest='rand_frac',
                        type=FloatRange(0, 1.0), default=1.0,
                        help='Fraction of input file to analyse at random.'
                             '  Default=%(default)s.')
    parser.add_argument('--max-heavies', dest='max_heavies',
                        type=IntRange(1), default=8,
                        help='Maximum number of heavy atoms in linker,'
                             ' excluding the 2 ring atoms being linked.'
                             '  Default=%(default)s.')
    parser.add_argument('--max-bonds', dest='max_bonds',
                        type=IntRange(2), default=5,
                        help='Maximum number of bonds in shorted path between'
                             ' ring atoms being linked.  Default=%(default)s.')
    parser.add_argument('--min-count', dest='min_count',
                        type=IntRange(1), default=2,
                        help='Minimum number of molecules a linker must be in'
                             ' to be reported.  Default=%(default)s.')
    parser.add_argument('-O', '--output-file', dest='output_files',
                        required=True, action='append',
                        help='Name of output file.  Should have extension'
                             ' .csv, .html, .json, .db or .sdf.  Linkers will'
                             ' be output in descending order of frequency.'
                             '  Multiple instances allowed, so all files can'
                             ' be output in one run. Extension .db is for a'
                             ' SQLite database, which also contains'
                             ' information about the molecules that the'
                             ' linkers were derived from.')
    parser.add_argument('--image-dir', dest='image_dir',
                        default='LinkerImages',
                        help='If writing an HTML file for the output, where to'
                             ' put the SVGs of the linkers.')
    parser.add_argument('--num-procs', dest='num_procs',
                        default=1, type=IntRange(1, cpu_count() - 1),
                        help='Number of processors if running in parallel.')

    args = parser.parse_args()
    for output_file in args.output_files:
        out_sfx = Path(output_file).suffix
        if out_sfx not in ['.csv', '.html', '.json', '.db', '.sdf']:
            print(f'ERROR bad suffix for output file. You gave {out_sfx} not one'
                  f' of ".csv", ".html" or ".json".')
            args = None
            break
    return args


def linker_subs_ok(linker_frags: list[list[Chem.Mol]]) -> bool:
    """
    Substituents can't be more than 1 atom unless they're a ring. So
    *C1COC(*)CC1 is ok but *C1COC(*)C(CC)C1 isn't.
    Args:
        linker_frags:

    Returns:
        True or False
    """
    # if more than 1 frag, only the one with the dummies should
    # be larger than 1 atom
    for linker_frag in linker_frags:
        if len(linker_frag) > 1:
            for i, frag in enumerate(linker_frag):
                dummy_frag = False
                for atom in frag.GetAtoms():
                    if not atom.GetAtomicNum():
                        dummy_frag = True
                        break
                if not dummy_frag and frag.GetNumAtoms() > 1:
                    return False

    return True


def shortest_path_between_dummies(linker: Chem.Mol) -> list[int]:
    """
    Returns the shortest path between 2 dummies in the molecule.
    Raises ValueError if there aren't exactly 2 dummies.
    Args:
        linker: molecule of interest

    Returns:
        list at atom indices describing the shortest path between
        dummies.
    """
    dummies = []
    for at in linker.GetAtoms():
        if at.GetAtomicNum() == 0:
            dummies.append(at.GetIdx())
    if len(dummies) != 2:
        raise ValueError(f'Molecule has {len(dummies)} not 2.')
    shortest_path = rdmolops.GetShortestPath(linker, dummies[0], dummies[1])
    return list(shortest_path)


def split_linker(mol: Chem.Mol) -> list[list[Chem.Mol]]:
    """
    Split the molecule into 2 pieces by breaking the first non-ring
    bond it finds where one end has degree > 2 that doesn't split
    it into 2 mols with a dummy in each.
    Returns the fragments as a list of mols.
    Args:
        mol:

    Returns:

    """
    # print(f'split_linker : {Chem.MolToSmiles(mol)}')
    frag_sets = []
    for bond in mol.GetBonds():
        if not bond.IsInRing():
            beg_atom = bond.GetBeginAtom()
            if not beg_atom.GetAtomicNum():
                continue
            end_atom = bond.GetEndAtom()
            if not end_atom.GetAtomicNum():
                continue
            if beg_atom.GetDegree() > 2 or end_atom.GetDegree() > 2:
                # print(f'splitting on {beg_atom.GetIdx()} -> {end_atom.GetIdx()}'
                #       f'  {beg_atom.GetAtomicNum()} and {end_atom.GetAtomicNum()}'
                #       f' :: {bond.GetBondType()}')
                fragged_mol = rdmolops.FragmentOnBonds(mol, [bond.GetIdx()],
                                                       addDummies=False)
                frags = rdmolops.GetMolFrags(fragged_mol, asMols=True,
                                             sanitizeFrags=False)
                for frag in frags:
                    num_dummies = 0
                    for atom in frag.GetAtoms():
                        if not atom.GetAtomicNum():
                            num_dummies += 1
                    if num_dummies and num_dummies == 2:
                        frag_sets.append(list(frags))
    if not frag_sets:
        frag_sets = [[Chem.Mol(mol)]]
    return frag_sets


def linker_is_ok(linker: Chem.Mol, max_heavies: int) -> bool:
    """
    Returns True if the linker is acceptable, False otherwise.
    That means:
    a. No more than max_heavies, excluding dummies
    b. No substituent more tnan 1 atom
    c. If more than 3 atoms excluding dummies, must contain
       one of a ring bond, a multiple bond or a non-carbon
       atom
    The criterion for maximum number of bonds between ring atoms
    should already have been enforced.
    """
    if linker.GetNumAtoms() > max_heavies + 2:
        return False

    linker_frags = split_linker(linker)
    # print(f'linker : {Chem.MolToSmiles(linker)}  fragged :', end='')
    # for linker_frag in linker_frags:
    #     print(' :: ', end='')
    #     for f in linker_frag:
    #         print(f' {Chem.MolToSmiles(f)}', end='')
    # print('')

    # if it didn't fragment, then there are no substituents to
    # worry about.
    if linker_frags and len(linker_frags[0]) > 1:
        if linker_subs_ok(linker_frags):
            return True

    # Now, if it doesn't have any substituents and it's more than 5
    # atoms including dummies it must contain at least a ring atom,
    # a double bond or a non-carbon.
    if len(linker_frags) == 1 and len(linker_frags[0]) == 1:
        if linker_frags[0][0].GetNumAtoms() < 6:
            return True
        for atom in linker_frags[0][0].GetAtoms():
            if atom.GetAtomicNum() and atom.GetAtomicNum() != 6:
                return True
        for bond in linker_frags[0][0].GetBonds():
            if (bond.GetBondType() != Chem.rdchem.BondType.SINGLE
                    or bond.IsInRing()):
                return True

    return False


def split_molecule(mol: Chem.Mol, bond1: Chem.Bond,
                   bond2: Chem.Bond) -> tuple[Chem.Mol, Chem.Mol, Chem.Mol]:
    """
    Split the molecule at the 2 bonds, returning the left molecule,
    the linker and the right, in that order.  Dummies marked with
    atom map numbers 1 and 2.  If there are already dummies marked
    like that, it will be a mess and probably fail.
    """
    split_mol = rdmolops.FragmentOnBonds(mol, [bond1.GetIdx(), bond2.GetIdx()],
                                         dummyLabels=[(1, 1), (2, 2)])
    split_frags = rdmolops.GetMolFrags(split_mol, asMols=True)
    left_mol = None
    right_mol = None
    linker = None
    for sf in split_frags:
        # print(f'split frag : {Chem.MolToSmiles(sf)}')
        dummies = []
        for atom in sf.GetAtoms():
            if not atom.GetAtomicNum() and atom.GetIsotope():
                dummies.append(atom.GetIsotope())
                atom.SetAtomMapNum(atom.GetIsotope())
                atom.SetIsotope(0)
        # print(f'{Chem.MolToSmiles(sf)} : {dummies}')
        if len(dummies) == 1:
            if dummies[0] == 1:
                left_mol = sf
            elif dummies[0] == 2:
                right_mol = sf
        elif len(dummies) == 2:
            linker = sf

    # print(f'left, linker, right : {Chem.MolToSmiles(left_mol)},'
    #       f' {Chem.MolToSmiles(linker)}, {Chem.MolToSmiles(right_mol)}')
    return left_mol, linker, right_mol


def find_linker_bonds(mol: Chem.Mol, max_length: int) -> list[tuple[Chem.Bond, Chem.Bond, int]]:
    """
    Find all the pairs of bonds that can be the start and end of a
    linker.
    This means one end of each is a ring atom not in the same SSSR ring,
    no more than max_length bonds apart, and neither is a ring bond
    itself.
    Returns a list of the 2 bonds and the distance between the ring
    atoms.
    """
    dist_mat = rdmolops.GetDistanceMatrix(mol)
    # The SSSRs in GetRingInfo are useful but not the final answer.
    # For example, in the 3 ring fused system C12COCCC2CC3CCCNC3C1 the
    # O and N atoms (indices 2 and 11) are recorded as not being in
    # the same ring even though they are both in the outer 14 atom
    # ring.
    sssr_info = mol.GetRingInfo()
    linker_bonds = []
    for at1 in mol.GetAtoms():
        if not at1.IsInRing():
            continue
        at1_bonds = [b for b in at1.GetBonds() if not b.IsInRing()]
        if not at1_bonds:
            continue
        at1_idx = at1.GetIdx()
        for at2 in mol.GetAtoms():
            if not at2.IsInRing():
                continue
            at2_idx = at2.GetIdx()
            if at2_idx <= at1_idx:
                continue
            length = dist_mat[at1_idx][at2_idx]
            if length == 1 or length > max_length:
                continue
            if sssr_info.AreAtomsInSameRing(at1_idx, at2_idx):
                continue
            at2_bonds = [b for b in at2.GetBonds() if not b.IsInRing()]
            if not at2_bonds:
                continue
            for b1 in at1_bonds:
                for b2 in at2_bonds:
                    # ring atoms must be further apart than the
                    # non-ring ends so they're pointing towards
                    # each other
                    e_at1 = b1.GetOtherAtomIdx(at1_idx)
                    e_at2 = b2.GetOtherAtomIdx(at2_idx)
                    if length > dist_mat[e_at1][e_at2]:
                        linker_bonds.append((b1, b2, length))

    return linker_bonds


def count_donors_and_acceptors(mol: Chem.Mol) -> tuple[int, int]:
    """
    Count the number of donors and acceptors in the molecule.  Since
    this is intended to be used primarily for linkers, where we lack
    a detailed environment for the atoms, it just classes all O and N
    atoms as acceptors, and all O and N with at least one H as donors.
    This makes them equivalent to the Lipinski definitions.
    Args:
        mol: the molecule of interest

    Returns:
        tuple of donor count, acceptor count
    """
    num_donors = 0
    num_acceptors = 0
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 7 or atom.GetAtomicNum() == 8:
            num_acceptors += 1
            if atom.GetTotalNumHs():
                num_donors += 1

    return num_donors, num_acceptors


def find_linkers(mol_rec: tuple[str, str], max_heavies: int = 8,
                 max_length: int = 5) -> tuple[str, list[Linker]]:
    """
    Return the SMILES strings of all linkers found in the molecule,
    along with the left hand (attached to the [*:1]) and right hand
    (attached to the [*:2]) pieces.
    Args:
        mol_rec: molecule of interest as a tuple of SMILES string, name
        max_heavies: maximum number of heavy atoms in a linker
        max_length: maximum number of bonds in shortest path
                    between R group atoms.

    Returns:
        list of dicts containing the linkers and details about them.
    """
    # print(f'find_linkers for {mol_rec[0]} : {mol_rec[1]}')
    mol = Chem.MolFromSmiles(mol_rec[0])
    if mol is None or not mol:
        return mol_rec[1], []

    linker_bonds = find_linker_bonds(mol, max_length)

    all_linkers = []
    for bond_pair in linker_bonds:
        # print(f'Splitting on {bond1.GetBeginAtomIdx()} -> {bond1.GetEndAtomIdx()}'
        #       f'  and {bond2.GetBeginAtomIdx()} -> {bond2.GetEndAtomIdx()}')
        left, linker, right = split_molecule(mol, bond_pair[0], bond_pair[1])
        ok = linker_is_ok(linker, max_heavies)
        num_donors, num_acceptors = count_donors_and_acceptors(linker)
        if ok:
            # print(f'Linker {Chem.MolToSmiles(linker)} ok.')
            mol_smiles = Chem.MolToSmiles(mol)
            lnk = Linker(name=mol_rec[1],
                         linker_smiles=Chem.MolToSmiles(linker),
                         mol_smiles=mol_smiles,
                         left_smiles=Chem.MolToSmiles(left),
                         right_smiles=Chem.MolToSmiles(right),
                         num_atoms=linker.GetNumAtoms(),
                         linker_length=bond_pair[2],
                         num_donors=num_donors,
                         num_acceptors=num_acceptors)
            if lnk not in all_linkers:
                all_linkers.append(lnk)
            if not lnk.symmetrical:
                # print(f'not symmetrical : {lnk.linker_smiles}')
                rlsmi, rrsmi = lnk.reversed_sides
                # print(f'orig     : {lnk.linker_smiles} {lnk.left_smiles} {lnk.right_smiles}')
                # print(f'reversed : {lnk.reversed_linker} {rlsmi} {rrsmi}')
                rlnk = Linker(name=mol_rec[1],
                              linker_smiles=lnk.reversed_linker,
                              mol_smiles=mol_smiles,
                              left_smiles=rlsmi,
                              right_smiles=rrsmi,
                              num_atoms=lnk.num_atoms,
                              linker_length=bond_pair[2],
                              num_donors=num_donors,
                              num_acceptors=num_acceptors)
                # print(f'final    : {lnk.linker_smiles} {lnk.left_smiles} {lnk.right_smiles}')
                if rlnk not in all_linkers:
                    all_linkers.append(rlnk)

    # print(f'{Chem.MolToSmiles(mol)} gave {len(all_linkers)}')
    return mol_rec[1], all_linkers


def get_mol_name(mol: Chem.Mol, mol_num: int, stem: str) -> str:
    """
    Return a name for the molecule.  If there's a _Name prop, use that,
    if not make something from the stem an molecule number as passed
    in.
    Args:
        mol:
        mol_num:
        stem:

    Returns:
    name
    """
    try:
        mol_name = mol.GetProp('_Name')
    except KeyError:
        mol_name = f'{stem}{mol_num}'

    return mol_name


def serial_find_all_linkers(input_mols: list[tuple[str, str]], max_heavies: int = 8,
                            max_length: int = 5, silent: bool = False) -> Optional[dict[list[Linker]]]:
    """
    Do the linker search in serial
    Args:
        input_mols:
        max_heavies:
        max_length:
        silent: if True, doesn't use the progress bar

    Returns:

    """
    all_linkers = defaultdict(list)

    for i, mol_rec in enumerate(input_mols):
        mol_name, linkers = find_linkers(mol_rec, max_heavies, max_length)
        for lnk in linkers:
            all_linkers[lnk.linker_smiles].append(lnk)

    return all_linkers


def parallel_find_all_linkers(input_mols: list[tuple[str, str]], num_procs: int,
                              max_heavies: int = 8, max_length: int = 5,
                              silent: bool = False) -> Optional[dict[str, list[Linker]]]:
    """
    Do the linker search in parallel using num_procs processors
    Args:
        input_mols:
        max_heavies:
        max_length:
        num_procs:
        silent: if True, doesn't use the progress bar

    Returns:

    """
    all_linkers = defaultdict(list)
    # chunk_size = 100000
    # next_start = 0
    # num_mols = len(input_mols)
    # with tqdm(total=len(input_mols)) as pbar:
    #     while next_start < num_mols:
    #         with cf.ProcessPoolExecutor(max_workers=cpu_count() - 1) as pool:
    #             futures = []
    #             print(f'Submitting {next_start} to {next_start + chunk_size}')
    #             for i, mol in enumerate(input_mols[next_start:next_start + chunk_size]):
    #                 fut = pool.submit(find_linkers, mol, next_start + i,
    #                                   max_heavies, max_length)
    #                 futures.append(fut)
    #
    #             for fut in cf.as_completed(futures):
    #                 mol_name, linkers = fut.result()
    #                 for linker in linkers:
    #                     all_linkers[linker.linker_smiles].append(linker)
    #                 pbar.update(1)
    #         next_start += chunk_size

    if len(input_mols) > 20:
        chunk_size = 1
    else:
        chunk_size = 20
    with cf.ProcessPoolExecutor(max_workers=cpu_count() - 1) as pool:
        for mol_name, linkers in cf.map(find_linkers, input_mols,
                                             repeat(max_heavies),
                                             repeat(max_length),
                                             max_workers=num_procs,
                                             chunksize=chunk_size,
                                             disable=silent):
            for linker in linkers:
                all_linkers[linker.linker_smiles].append(linker)
    return all_linkers


def find_all_linkers(input_mols: list[tuple[str, str]], max_heavies: int = 8,
                     max_length: int = 5, num_procs=1,
                     silent: bool = False) -> Optional[dict[str, list[Linker]]]:
    """
    Analyse the structures in the file and return linkers found.
    Linkers will have maximum of max_heavies non-H atoms, and the
    topological, through-bond distance between the 2 R groups is no
    more than max_length.
    The atoms marking the ends are ring atoms.

    Args:
        input_mols: the molecules to process, as tuples of SMILES and
                    name
        max_heavies: maximum number of heavy atoms in a linker
        max_length: maximum number of bonds in shortest path
                    between R group atoms.
        num_procs: number of processors for parallel runs
        silent: suppresses progress reports.

    Returns:
        List of Mols that are the linkers or None if there was a
        problem.
    """

    if not silent:
        print(f'Processing {len(input_mols)} molecules.')
    if num_procs == 1:
        all_linkers = serial_find_all_linkers(input_mols, max_heavies,
                                              max_length, silent)
    else:
        all_linkers = parallel_find_all_linkers(input_mols, num_procs,
                                                max_heavies, max_length,
                                                silent)

    return all_linkers


def sort_linkers(linkers: dict[str, list[Linker]]) -> dict[str, list[Linker]]:
    """
    Sort the linkers in descending size of frequency, assumed to be the
    length of the associated list.
    Args:
        linkers:

    Returns:
        same info, different order
    """
    freqs = []
    for smi, mols in linkers.items():
        freqs.append((smi, len(mols)))
    freqs.sort(key=lambda fr: fr[1], reverse=True)

    sorted_linkers = {}
    for f in freqs:
        sorted_linkers[f[0]] = linkers[f[0]]

    return sorted_linkers


def make_image(mol: Chem.Mol, image_dir: Union[str, Path]) -> Path:
    """
    Make an SVG of the molecule in image_dir, named by the mol name
    which is assumed to be unique, and returns name of file.
    image_dir is created if it doesn't exist.
    Args:
        mol:
        image_dir:

    Returns:
        name of image file created.
    """
    image_path = Path(image_dir)
    if not image_path.exists():
        image_path.mkdir()

    image_width = -1
    image_height = -1
    drawer = rdMolDraw2D.MolDraw2DSVG(image_width, image_height)
    rdDepictor.Compute2DCoords(mol)
    rdMolDraw2D.DrawMoleculeACS1996(drawer, mol)
    drawer.FinishDrawing()
    img_file = image_path / f'{mol.GetProp("_Name")}.svg'
    with open(img_file, 'w') as f:
        f.write(drawer.GetDrawingText())

    return img_file


def start_html_file(html_file: str) -> IO:
    htmlf = open(html_file, 'w')

    htmlf.write(f'''<!DOCTYPE html>
    <html lang="en">
      <head>
        <meta charset="UTF-8" />
        <meta name="viewport" content="width=device-width,
              initial-scale=1.0" />
        <meta http-equiv="X-UA-Compatible" content="ie=edge" />
        <link href="index.css" rel="stylesheet" />
        <style>
            .fixTableHead {{
            overflow-y: auto;
            height: 95vh;
            }}
            .fixTableHead thead th {{
            position: sticky;
            top: 0;
            }}
            table {{
            border-collapse: collapse;		
            width: 100%;
            }}
            th,
            td {{
            padding: 8px 15px;
            border: 2px solid #529432;
            }}
            th {{
            background: #ABDD93;
            }}
            .plotCell {{
            width: 40%;
            height: auto;
            }}
            .plotImage {{
            width: 100%;
            height: auto;
            object-fit: scale-down
            }}
        </style>
      </head>
      <body>
        <div class="fixTableHead" style='float:right; width:100%;'>
            <table border="1">
              <thead>
                  <tr>
                    <th>Name</th>
                    <th>Linker</th>
                    <th>Occurrences</th>
                  </tr>
              </thead>
              <tbody>\n''')
    return htmlf


def finish_html_file(htmlf: IO) -> None:
    """
    Write the end of the table and close the file.
    Args:
        htmlf:

    Returns:

    """
    htmlf.write('''        </tbody>
          </table>
          </div>
        </body>
      </html>\n''')
    htmlf.close()


def write_html_linkers(linkers: dict[str, list[str]],
                       html_file: Union[str, Path],
                       image_dir: Path) -> bool:
    """
    Write linkers to HTML file containing images.
    Args:
        linkers:
        html_file:
        image_dir: where to put the images

    Returns:
        bool on success
    """
    try:
        htmlf = start_html_file(html_file)
    except FileNotFoundError:
        print(f'Failed to open {html_file} for writing.')
        return False

    i = 1
    for smi, mols in linkers.items():
        mol = Chem.MolFromSmiles(smi)
        mol_name = f'Linker_{i}'
        mol.SetProp('_Name', mol_name)
        image_file = make_image(mol, image_dir)
        htmlf.write('<tr>')
        htmlf.write(f'<td>{mol_name}</td>')
        htmlf.write(f'<td><img src={image_file} /></td>')
        htmlf.write(f'<td>{len(mols)}</td>')
        htmlf.write('</tr>')
        i += 1

    finish_html_file(htmlf)
    return True


def write_csv_linkers(linkers: dict[str, list[Linker]],
                      csv_file: Union[str, Path]) -> bool:
    """
    Write linkers to CSV file
    Args:
        linkers:
        csv_file:

    Returns:
        bool on success
    """
    try:
        with open(csv_file, 'w', newline='') as fo:
            csvw = csv.writer(fo)
            csvw.writerow(['Name', 'SMILES', 'Frequency'])
            i = 1
            for smi, mols in linkers.items():
                mol_name = f'Linker_{i}'
                csvw.writerow([mol_name, smi, len(mols)])
                i += 1
    except FileNotFoundError:
        print(f'Failed to open {csv_file} for writing.')
        return False

    return True


def create_linkers_table(conn: sqlite3.Connection,
                         linkers: dict[str, list[Linker]]) -> None:
    """
    Create and populate the linkers table, containing the basic
    information about each linker - SMILES, label and
    number of times it occurred in the input database.  The SMILES is
    the primary key for the table

    Args:
        conn: connection to the database
        linkers: the linkers data

    Returns:
        None
    """
    curs = conn.cursor()
    curs.execute('CREATE TABLE linkers'
                 ' (name TEXT PRIMARY KEY,'
                 ' smiles TEXT,'
                 ' frequency INTEGER)')

    i = 1
    linker_rows = []
    for smi, mols in linkers.items():
        l_name = f'Linker_{i}'
        linker_rows.append((l_name, smi, len(mols)))
        i += 1

    curs.executemany('INSERT INTO linkers'
                     ' (name, smiles, frequency)'
                     ' VALUES (?,?,?)', linker_rows)
    conn.commit()


def create_linker_molecules_table(conn: sqlite3.Connection,
                                  linkers: dict[str, list[str]]) -> None:
    """
    Make and populate the table giving the molecule details for each
    molecule that each linker was found in.  The SMILES is a foreign
    key to the linkers table.
    Args:
        conn:
        linkers:

    Returns:

    """
    print('Creating linker_molecules table')
    curs = conn.cursor()
    curs.execute('CREATE TABLE linker_molecules '
                 '(mol_name TEXT,'
                 ' mol_smiles TEXT,'
                 ' linker_smiles TEXT,'
                 ' left_smiles TEXT,'
                 ' right_smiles TEXT,'
                 ' symmetrical BOOLEAN,'
                 ' FOREIGN KEY (linker_smiles) REFERENCES linkers(smiles))')
    for smi, linkers in linkers.items():
        mol_rows = []
        for linker in linkers:
            mol_rows.append((linker.name, linker.mol_smiles, linker.linker_smiles,
                             linker.left_smiles, linker.right_smiles))
        curs.executemany('INSERT INTO linker_molecules'
                         ' (mol_name, mol_smiles, linker_smiles, left_smiles,'
                         '  right_smiles)'
                         ' VALUES (?,?,?,?,?)', mol_rows)
        conn.commit()


def write_sqlite_linkers(linkers: dict[str, list[str]],
                         db_file: Union[str, Path]) -> bool:
    """
    Write linkers to SQLite file.  If file exists, it is overwritten.
    Args:
        linkers:
        db_file:

    Returns:
        bool on success
    """
    if db_file.exists():
        print(f'WARNING - {db_file} will be overwritten.')
        db_file.unlink()

    try:
        conn = sqlite3.connect(db_file)
    except sqlite3.OperationalError:
        print(f'Failed to open database file SQLite {db_file} for writing.')
        return False

    create_linkers_table(conn, linkers)
    create_linker_molecules_table(conn, linkers)
    conn.close()

    return True


def write_json_linkers(linkers: dict[str, list[Linker]],
                       json_file: Union[str, Path]) -> bool:
    """
    Write linkers to JSON file
    Args:
        linkers:
        json_file:

    Returns:
        bool on success
    """
    json_list = []
    i = 1
    for smi, mols in linkers.items():
        mol_name = f'Linker_{i}'
        json_list.append({
            'SMILES': smi,
            'Name': mol_name,
            'Frequency': len(mols)
        })
        i += 1

    try:
        with open(json_file, 'w', newline='') as fo:
            fo.write(json.dumps(json_list, indent=2))
            fo.write('\n')
    except FileNotFoundError:
        print(f'Failed to open {json_file} for writing.')
        return False

    return True


def write_sdf_linkers(linkers: dict[str, list[str]],
                      sdf_file: Union[str, Path]) -> bool:
    """
    Write linkers to HTML file containing images.
    Args:
        linkers:
        sdf_file:

    Returns:
        bool on success
    """
    try:
        writer = Chem.SDWriter(str(sdf_file))
    except (FileNotFoundError, OSError):
        print(f'Failed to open {sdf_file} for writing.')
        return False

    i = 1
    for smi, mols in linkers.items():
        mol = Chem.MolFromSmiles(smi)
        mol_name = f'Linker_{i}'
        mol.SetProp('_Name', mol_name)
        mol.SetProp('Frequency', f'{len(mols)}')
        mol.SetProp('SMILES', smi)
        writer.write(mol)
        i += 1

    return True


def remove_atom_maps(mol: Chem.Mol) -> None:
    for atom in mol.GetAtoms():
        atom.SetAtomMapNum(0)


def read_molecules(infile: str, rand_frac: float) -> Optional[list[tuple[str, str]]]:
    """
    Read all molecules up front, so we can process them in parallel.
    Save them as SMILES strings to reduce the memory use at the expense
    of a bit of extra processing.  When running in parallel, memory is
    pretty critical.
    Args:
        infile: SMILES, SDF or gzipped SDF.
        rand_frac: random fraction of input file to process

    Returns:
        list of SMILES and names, the latter created if missing.
    """
    molecules = []
    suppl = create_mol_supplier(infile)
    if suppl is None:
        return None

    for i, mol in tqdm(enumerate(suppl)):
        if not mol or not mol.GetNumAtoms():
            continue
        if random.random() > rand_frac:
            continue
        mol_name = get_mol_name(mol, i, 'Str_')
        mol.SetProp('_Name', mol_name)
        remove_atom_maps(mol)
        molecules.append((Chem.MolToSmiles(mol), mol_name))

    print(f'Read {len(molecules)} from {infile}')
    return molecules


def write_linkers_files(linkers: dict[str, list[Linker]],
                        output_files: list[str],
                        image_dir: Union[str, Path]) -> bool:
    """
    Write one or more output files as named, with format given by
    suffix.  Suffix should already be enforced as one of '.csv',
    '.html' or '.json'.
    Args:
        linkers: dictionary of linkers, keyed on SMILES string,
                 with list of names of molecules containing that
                 linker
        output_files: one or more output file names
        image_dir: where to put the images for the HTML table.


    Returns:
        bool on success
    """

    for output_file in output_files:
        output_path = Path(output_file)
        if output_path.suffix == '.html':
            res = write_html_linkers(linkers, output_path, image_dir)
            if not res:
                return False
        elif output_path.suffix == '.csv':
            res = write_csv_linkers(linkers, output_path)
            if not res:
                return False
        elif output_path.suffix == '.json':
            res = write_json_linkers(linkers, output_path)
            if not res:
                return False
        elif output_path.suffix == '.db':
            res = write_sqlite_linkers(linkers, output_path)
            if not res:
                return False
        elif output_path.suffix == '.sdf':
            res = write_sdf_linkers(linkers, output_path)
            if not res:
                return False

    return True


def trim_linkers(linkers: dict[str, list[Linker]], min_count: int) -> dict[str, list[Linker]]:
    """
    Remove linkers that are in fewer than min_count molecules and sorts
    them in the dict into descending order of size.
    Args:
        linkers:
        min_count:

    Returns:
        trimmed list
    """
    trimmed_linkers = {}
    for smi, linkers in linkers.items():
        if len(linkers) >= min_count:
            trimmed_linkers[smi] = linkers
    sorted_linkers = sort_linkers(trimmed_linkers)
    return sorted_linkers


def main():
    print(f'Using RDKit version {rdBase.rdkitVersion}.')
    # so that properties, such as the _Name, are pickled when passed
    # into the multiprocessing bit.
    Chem.SetDefaultPickleProperties(Chem.PropertyPickleOptions.AllProps)

    args = parse_args()
    if args is None:
        return False

    input_mols = read_molecules(args.infile, args.rand_frac)
    if input_mols is None:
        return False

    linkers = find_all_linkers(input_mols, args.max_heavies, args.max_bonds,
                               args.num_procs)
    if linkers is None:
        return False
    linkers = trim_linkers(linkers, args.min_count)
    write_linkers_files(linkers, args.output_files, args.image_dir)

    return True


if __name__ == '__main__':
    sys.exit(not main())
