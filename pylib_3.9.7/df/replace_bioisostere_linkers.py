#!/usr/bin/env python

# Use a bioisosteric linkers database to replace linkers in input
# molecules

import argparse
import re
import sqlite3
import sys

from pathlib import Path
from typing import Optional

from rdkit import rdBase, Chem
from rdkit.Chem import rdmolops

import df.find_linkers as fl
from df.int_range import IntRange


def parse_args(cli_args: list[str]):
    """
        Use argparse module to parse the command line arguments.

        Returns:
            (Namespace object): Parsed arguments in argparse Namespace
                                object.
        """
    parser = argparse.ArgumentParser(description='Replace linkers with'
                                                 ' bioisosteres.')
    parser.add_argument('-D', '--database-file', dest='db_file',
                        required=True,
                        help='Name of SQLite file containing the bioisosteric'
                             ' linkers.')
    in_args = parser.add_mutually_exclusive_group(required=True)
    in_args.add_argument('-S', '--query-smiles', dest='query_smiles',
                         help='Single SMILES string of molecule to be'
                              ' processed.  Will probably need to be in'
                              ' quotes.')
    in_args.add_argument('-I', '--input-file', dest='input_file',
                         help='Name of molecule file containing input'
                              ' molecules to be processed.')
    parser.add_argument('-O', '--output-file', dest='out_file',
                        required=True,
                        help='Name of output file.  Should be an SDF or'
                             ' SMILES file with extension .sdf or .smi.')
    parser.add_argument('--max-heavies', dest='max_heavies',
                        type=IntRange(1), default=8,
                        help='Maximum number of heavy atoms in linker,'
                             ' excluding the 2 ring atoms being linked.'
                             '  Default=%(default)s.')
    parser.add_argument('--max-bonds', dest='max_bonds',
                        type=IntRange(2), default=5,
                        help='Maximum number of bonds in shorted path between'
                             ' ring atoms being linked.  Default=%(default)s.')
    parser.add_argument('--plus-delta-length', dest='plus_length',
                        type=IntRange(-1, 5), default=-1,
                        help='Positive delta length for linker compared with'
                             ' that found in the query structure.  Default=-1'
                             ' means no maximum.')
    parser.add_argument('--minus-delta-length', dest='minus_length',
                        type=IntRange(-1, 5), default=-1,
                        help='Negative delta length for linker compared with'
                             ' that found in the query structure.  Default=-1'
                             ' means no minimum.')
    parser.add_argument('--match-donors', dest='match_donors',
                        action='store_true',
                        help='If True and the query linker has an hbond donor,'
                             ' the replacement must too, and not if not.  If'
                             ' False, it will take either.')
    parser.add_argument('--match-acceptors', dest='match_acceptors',
                        action='store_true',
                        help='If True and the query linker has an hbond'
                             ' acceptor, the replacement must too, and not if'
                             ' not.  If False, it will take either.')

    args = parser.parse_args(cli_args)
    return args


def check_db_file(db_file: str):
    """
    Makes sure the db_file is a valid bioisostere one.
    Args:
        db_file:

    Returns:
        bool

    Raises:
         FileNotFoundError if db_file doesn't exist.
         ValueError is db_file doesn't contain the relevant tables.
    """
    # sqlite3 helpfully creates a db file if it doesn't exist!
    if not Path(db_file).exists():
        raise FileNotFoundError(f'{db_file} not available for reading.')

    conn = sqlite3.connect(db_file)
    for table in ['bioisosteres', 'linkers']:
        res = conn.execute('SELECT COUNT(name)'
                           ' FROM sqlite_schema'
                           ' WHERE type="table" AND name=?',
                           (table,)).fetchone()
        if res[0] == 0:
            raise ValueError(f'Invalid database {db_file}: no table "{table}".')


def make_new_smiles(mol_smi: str, linker_smi: str, bios: list[str]) -> list[str]:
    """
    Take the fragmented SMILES string, the SMILES of the current
    linker of interest and a list of bioisostere SMILES and make SMILES
    strings of molecules where the linker is replaced by the
    bioisosteres.  The new SMILES strings are not zipped up, because
    there may be more substitutions to make.
    Args:
        mol_smi: fragmented SMILES string
        linker_smi: SMILES string of current linker of interest
        bios: list of SMILES strings to be attached to left and right
              SMILES

    Returns:
        list of new SMILES strings, still in fragmented form.
    """
    params = rdmolops.MolzipParams()
    params.label = rdmolops.MolzipLabel.AtomMapNumber
    new_smis = []
    for b in bios:
        new_smi = mol_smi.replace(linker_smi, b)
        new_smis.append(new_smi)
        # print(f'"{new_smis[-1]}",')

    return new_smis


def split_input_smiles(query_smiles: str, linker_smis: list[str],
                       max_heavies: int, max_bonds: int) -> str:
    """
    Take the molecule and split on any linkers to produce a
    fragmented SMILES with linkers and pieces, with the dummy map
    numbers adjusted so that multiple linkers don't have the same
    values.
    e.g. split c1ccccc1CCc1cnccc1OCOc1ccccc1 into
    [*:1]c1ccccc1.[*:1]CC[*:2].[*:2]c1cnccc1[*:3].[*:3]OCO[*:4].[*:4]c1ccccc1
    """
    new_smi = query_smiles
    _, linkers = fl.find_linkers((new_smi, ''), max_heavies=max_heavies, max_length=max_bonds)
    # print(f'depth = {len(linker_smis)} : num linkers {len(linkers)}')
    if not linkers:
        return query_smiles

    i = len(linkers) + 1
    new1 = f'[*:{2 * i - 1}]'
    new2 = f'[*:{2 * i}]'
    new_linker = linkers[0].linker_smiles.replace('[*:1]', new1).replace('[*:2]', new2)
    linker_smis.append(new_linker)
    new_left_smi = linkers[0].left_smiles.replace('[*:1]', new1)
    new_right_smi = linkers[0].right_smiles.replace('[*:2]', new2)
    # print(new_linker)
    # print(new_left_smi)
    # print(new_right_smi)
    new_left_smi = split_input_smiles(new_left_smi, linker_smis, max_heavies, max_bonds)
    new_right_smi = split_input_smiles(new_right_smi, linker_smis, max_heavies, max_bonds)

    new_smi = f'{new_left_smi}.{new_linker}.{new_right_smi}'
    return new_smi


def replace_linkers(query_smiles: str, db_file: str,
                    max_heavies: int, max_bonds: int,
                    plus_length: int, minus_length: int,
                    match_donors: bool, match_acceptors: bool) -> list[Chem.Mol]:
    """
    Take the query SMILES string, find any linkers in the structure and
    returns new molecules with the linkers replaced by ones plucked
    from the bioisostere database.  Assumes query_smiles is valid.
    New linkers must have a length (shortest distance in bonds between
    dummies) of l-minus_length to l+plus_length where l is the length
    of each linker in query_smiles.  If match_donors is True, then any
    replacement linkers must have a donor if the query linker does, and
    not if not, and likewise for the match_acceptors.  If either is
    False, it doesn't care.
    Args:
        query_smiles:
        db_file:
        max_heavies:
        max_bonds:
        plus_length:
        minus_length:
        match_donors:
        match_acceptors:

    Returns:
        new molecules
    """
    linker_smis = []
    split_smi = split_input_smiles(query_smiles, linker_smis, max_heavies,
                                   max_bonds)
    print(f'split_smi : {split_smi}')
    print(f'linker_smis : {linker_smis}')
    for l in linker_smis:
        print(f'linker_smi : {l}')

    new_smis = [split_smi]
    for lsmi in linker_smis:
        bios = fetch_bioisosteres(lsmi, db_file, plus_length, minus_length,
                                  match_donors, match_acceptors)
        # If there aren't any bioisosteres for this linker, leave it
        # as is.
        # print(f'linker {lsmi} has {len(bios)} replacements : {bios}')
        if not bios:
            bios = [lsmi]
        next_new_smis = []
        for nm in new_smis:
            next_new_smis.extend(make_new_smiles(nm, lsmi, bios))
        new_smis = next_new_smis

    new_mols = []
    for new_smi in new_smis:
        mol = Chem.MolFromSmiles(new_smi)
        if mol is None:
            print(f'Warning: SMILES {query_smiles} gave invalid product {new_smi}.')
            continue
        zip_mol = rdmolops.molzip(mol)
        new_mols.append(zip_mol)

    return new_mols


def bulk_replace_linkers(mol_file: str, db_file: str,
                         max_heavies: int, max_bonds: int,
                         plus_length: int, minus_length: int,
                         match_donors: bool,
                         match_acceptors: bool) -> Optional[list[Chem.Mol]]:
    """
    Take the structures in the mol file and process them with
    replace_linkers.  Returns None if file can't be read.
    New linkers must have a length (shortest distance in bonds between
    dummies) of l-minus_length to l+plus_length where l is the length
    of each linker in query_smiles.  If match_donors is True, then any
    replacement linkers must have a donor if the query linker does, and
    not if not, and likewise for the match_acceptors.  If either is
    False, it doesn't care.
    Args:
        mol_file:
        db_file:
        max_heavies:
        max_bonds:
        plus_length:
        minus_length:
        match_donors:
        match_acceptors:

    Returns:
        new molecules
    """
    suppl = fl.create_mol_supplier(mol_file)
    if suppl is None:
        return None

    all_new_mols = []
    for mol in suppl:
        if not mol or not mol.GetNumAtoms():
            continue
        print(f'Doing {mol.GetProp("_Name")}')
        mol_smi = Chem.MolToSmiles(mol)
        new_mols = replace_linkers(mol_smi, db_file, max_heavies, max_bonds,
                                   plus_length, minus_length, match_donors,
                                   match_acceptors)
        all_new_mols.extend(new_mols)

    return all_new_mols


def trim_linkers_by_length(conn: sqlite3.Connection, query_linker: str,
                           linkers: list[str], plus_length: int,
                           minus_length) -> list[str]:
    """
    Remove any linkers that have length outside the deltas given.
    """
    sql1 = """SELECT DISTINCT path_length, linker_smiles
     FROM linkers WHERE linker_smiles = ?"""

    row = conn.execute(sql1, (query_linker,)).fetchone()
    if plus_length == -1:
        max_length = 1000
    else:
        max_length = row[0] + plus_length
    if minus_length == -1:
        min_length = 0
    else:
        min_length = row[0] - minus_length

    sql2 = f"""SELECT DISTINCT path_length, linker_smiles
     FROM linkers WHERE
      linker_smiles IN ({','.join(['?' for _ in range(len(linkers))])})"""
    new_linkers = []
    for row in conn.execute(sql2, linkers):
        if min_length <= row[0] <= max_length:
            new_linkers.append(row[1])
    return new_linkers


def trim_linkers_by_hbonding(conn: sqlite3.Connection, query_linker: str,
                             linkers: list[str], match_donors: bool,
                             match_acceptors: bool) -> list[str]:
    """
    Remove any linkers that fail the hbonding tests.
    """
    sql1 = """SELECT DISTINCT num_donors, num_acceptors, linker_smiles
     FROM linkers WHERE linker_smiles = ?"""

    row = conn.execute(sql1, (query_linker,)).fetchone()
    num_donors = row[0]
    num_acceptors = row[1]

    if match_donors:
        check_donors = True
        if row[0]:
            must_have_donor = True
        else:
            must_have_donor = False
    else:
        check_donors = False
    if match_acceptors:
        check_acceptors = True
        if row[1]:
            must_have_acceptor = True
        else:
            must_have_acceptor = False
    else:
        check_acceptors = False

    print(f'trim_linkers_by_hbonding : {row[0]}  {row[1]} : {match_donors}'
          f' {num_donors} : {match_acceptors} {num_acceptors}')
    sql2 = f"""SELECT DISTINCT num_donors, num_acceptors, linker_smiles
     FROM linkers WHERE
      linker_smiles IN ({','.join(['?' for _ in range(len(linkers))])})"""
    new_linkers = []
    for row in conn.execute(sql2, linkers):
        if match_donors:
            if not (num_donors and row[0]) or (not num_donors and not row[0]):
                print(f'next row {row[0]} and {row[1]} skipping')
                continue
        if match_acceptors:
            if not (num_acceptors and row[1])\
                    or (not num_acceptors and not row[1]):
                print(f'next row {row[0]} and {row[1]} skipping')
                continue
        print(f'next row {row[0]} and {row[1]} keeping')
        new_linkers.append(row[2])
    return new_linkers


def fetch_bioisosteres(linker_smi: str, db_file: str,
                       plus_length: int, minus_length: int,
                       match_donors: bool, match_acceptors: bool) -> Optional[list[str]]:
    """
    Pull all bioisosteres from db_file that include the linker_smi
    on one side or the other.  Returns SMILES strings of the
    swaps.
    Bioisosteres must have a length (shortest distance in bonds between
    dummies) of l-minus_length to l+plus_length where l is the length
    of each linker in query_smiles.  If match_donors is True, then any
    replacement linkers must have a donor if the query linker does, and
    not if not, and likewise for the match_acceptors.  If either is
    False, it doesn't care.
    Args:
        linker_smi: assumed to contain 2 dummies for the link points
                    e.g. [*:1]C[*:2] and be in canonical SMILES to
                    match the linkers in the db_file.
        db_file: name of valid SQLite3 database.
        plus_length:
        minus_length:
        match_donors:
        match_acceptors:

    Returns:
        SMILES strings in a list

    Raises: FileNotFoundError if db_file doesn't exist.
    """
    check_db_file(db_file)
    # print(f'fetch_bioisosteres : {linker_smi}')
    # the linker_smi won't necessarily have the dummies with map
    # numbers 1 and 2, but they need to be so for looking up.
    dummy_regex = re.compile(r'\[\*:(\d+)\]')
    dummy_matches = sorted([int(dm) for dm in dummy_regex.findall(linker_smi)])
    # print(f'searching {linker_smi} : dummy matches : {dummy_matches}')
    new1 = f'[*:{dummy_matches[0]}]'
    new2 = f'[*:{dummy_matches[1]}]'
    mended_smi = linker_smi.replace(new1, '[*:1]').replace(new2, '[*:2]')

    conn = sqlite3.connect(db_file)
    sql = """SELECT linker1_smiles, linker2_smiles FROM bioisosteres
    WHERE linker1_smiles = ? or linker2_smiles = ?"""
    linkers = conn.execute(sql, (mended_smi, mended_smi)).fetchall()
    bios = []
    for linker in linkers:
        if linker[0] == mended_smi:
            bios.append(linker[1])
        else:
            bios.append(linker[0])

    bios = trim_linkers_by_length(conn, mended_smi, bios, plus_length,
                                  minus_length)

    if match_donors or match_donors:
        bios = trim_linkers_by_hbonding(conn, mended_smi, bios, match_donors,
                                        match_acceptors)

    bios = [b.replace('[*:1]', new1).replace('[*:2]', new2) for b in bios]
    # print(bios)
    return bios


def write_mols(mols: list[Chem.Mol], out_file: str) -> bool:
    """
    Writes the molecules to the named file.  Returns bool on success.
    Args:
        mols:
        out_file:

    Returns:

    """
    out_path = Path(out_file)
    if out_path.suffix == '.smi':
        writer = Chem.SmilesWriter(out_file, includeHeader=False)
    elif out_path.suffix == '.sdf':
        writer = Chem.SDWriter(out_file)
    else:
        print(f'ERROR : unrecognised extension for file {out_file}. Nothing'
              f' written.')
        return False

    for mol in mols:
        writer.write(mol)

    return True


def check_smiles(smiles: str) -> bool:
    """
    Return False if molecule can't be parsed successfully.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None or not mol.GetNumAtoms():
        return False
    else:
        return True


def main(cli_args):
    print(f'Using RDKit version {rdBase.rdkitVersion}.')
    args = parse_args(cli_args)
    if args is None:
        return False

    if args.query_smiles is not None:
        if not check_smiles(args.query_smiles):
            return False
        new_mols = replace_linkers(args.query_smiles, args.db_file,
                                   args.max_heavies, args.max_bonds,
                                   args.plus_length, args.minus_length,
                                   args.match_donors, args.match_acceptors)
    else:
        new_mols = bulk_replace_linkers(args.input_file, args.db_file,
                                        args.max_heavies, args.max_bonds,
                                        args.plus_length, args.minus_length)
    if new_mols is None or not new_mols:
        return False

    if not write_mols(new_mols, args.out_file):
        return False


if __name__ == '__main__':
    sys.exit(not main(sys.argv[1:]))
