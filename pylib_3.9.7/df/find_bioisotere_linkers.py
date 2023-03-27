#!/usr/bin/env python
#
# Process a ChEMBL database and pull out the bioisosteric linkers, as
# defined in
# https://www.sciencedirect.com/science/article/abs/pii/S0968089623000421#b0255

import argparse
import concurrent.futures as cf
import copy
import itertools
import json
import sqlite3
import sys

import df.find_linkers as fl

from os import close, cpu_count
from pathlib import Path
from tempfile import mkstemp
from typing import Any, IO, Optional, Union

from rdkit import rdBase, Chem
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D

from df.int_range import IntRange


class Bioisostere:
    """
    Class for holding a bioisostere, defined as 2 linkers (short pieces
    of molecule between 2 rings) where the rest of the 2 molecules are
    the same. I.e. if you take linker 1 out of molecule 1 and replace
    linker 2 in molecule 2 with it, you would get molecule 1.  Keeps a
    list of the molecules that gave rise to it.
    """

    def __init__(self, linker1: fl.Linker, linker2: fl.Linker,
                 doc: str):
        """
        Takes 2 linkers and the ID of the document they came from.
        The latter doesn't matter except that we want ot keep track
        of how many different documents the examples came from.  Ertl
        only counts them when they come from at least 2 documents.
        Args:
            linker1:
            linker2:
            doc:
        """
        # SMILES strings of the 2 linkers, in defined ordered for
        # better duplicate matching.
        if linker1.linker_smiles < linker2.linker_smiles:
            self._linker1 = linker1.linker_smiles
            self._linker2 = linker2.linker_smiles
        else:
            self._linker1 = linker2.linker_smiles
            self._linker2 = linker1.linker_smiles
        self._rev_linker1, self._rev_linker2 = self.reversed_linkers()
        # the examples are the tuples of linkers that produced the
        # bioisostere, so we can go back and probe the molecular
        # environments they came from, for example, and the doc ID.
        # Keep a separate set of the docs that the examples come from
        # for convenience.
        self._docs = set()
        self._examples = []
        self.add_example(linker1, linker2, doc)

    def __eq__(self, other: Any) -> bool:
        if not isinstance(other, Bioisostere):
            return NotImplemented

        if (self._linker1 == other._linker1
                and self._linker2 == other._linker2):
            return True
        if (self._rev_linker1 == other._linker1
                and self._rev_linker2 == other._linker2):
            return True

        return False

    def __str__(self) -> str:
        return (f'Linker 1 : {self._linker1}  Linker 2 : {self._linker2}'
                f' Num Examples : {self.num_examples}')

    @property
    def linker1_smiles(self) -> str:
        return self._linker1

    @property
    def linker2_smiles(self) -> str:
        return self._linker2

    @property
    def num_examples(self) -> int:
        return len(self._examples)

    @property
    def examples(self) -> list[tuple[fl.Linker, fl.Linker, str]]:
        return self._examples

    @property
    def num_docs(self) -> int:
        return len(self._docs)

    def add_example(self, linker1: fl.Linker, linker2: fl.Linker,
                    doc: str, dont_sort: bool = False) -> None:
        """
        Add the 2 linkers and document ID to the example set if they're
        not in already.  Always keep them so that linker1.linker_smiles
        sorts smaller.
        Args:
            linker1:
            linker2:
            doc:
            dont_sort: by default, the examples are sorted by doc. This
                       turns that off, for when it's called by
                       add_examples, which will sort them at the end.

        Returns:
            None
        """

        if linker1.linker_smiles > linker2.linker_smiles:
            linker1, linker2 = linker2, linker1
        found_it = False
        for ex in self._examples:
            if (linker1.mol_smiles == ex[0].mol_smiles
                    and linker2.mol_smiles == ex[1].mol_smiles
                    and doc == ex[2]):
                found_it = True
        if not found_it:
            self._examples.append((linker1, linker2, doc))
            self._docs.add(doc)
            if not dont_sort:
                self._sort_examples()

    def add_examples(self, new_exs: list[tuple[fl.Linker, fl.Linker, str]]) -> None:
        """
        Add a whole bunch of new examples, probably from another
        Bioisostere that is being merged in.
        Args:
            new_exs:

        Returns:

        """
        for new_ex in new_exs:
            self.add_example(new_ex[0], new_ex[1], new_ex[2], dont_sort=True)
        self._sort_examples()

    def reversed_linkers(self) -> tuple[str, str]:
        """
        Make SMILES strings of the reversed linkers.
        Swap left and right for linker1 and compare to linker2.  It's
        almost always the case that changing the map number doesn't
        change the canonical SMILES (apart from the obvious), but
        from time to time it does e.g.
        O=C1C2CCCCC2C(=O)N1CC(=O)N1CCCSC1=[*:2] goes to
        O=C1C2CCCCC2C(=O)N1CC(=O)N1CCCSC1=[*:1] canonicalises to
        O=C(CN1C(=O)C2CCCCC2C1=O)N1CCCSC1=[*:1]
        """
        ol1 = self._linker1.replace('[*:1]', '[*:XX]')
        ol1 = ol1.replace('[*:2]', '[*:1]')
        ol1 = ol1.replace('[*:XX]', '[*:2]')
        ol1 = Chem.MolToSmiles(Chem.MolFromSmiles(ol1))
        ol2 = self._linker2.replace('[*:1]', '[*:XX]')
        ol2 = ol2.replace('[*:2]', '[*:1]')
        ol2 = ol2.replace('[*:XX]', '[*:2]')
        ol2 = Chem.MolToSmiles(Chem.MolFromSmiles(ol2))

        if ol1 < ol2:
            return ol1, ol2
        else:
            return ol2, ol1

    def _sort_examples(self) -> None:
        self._examples.sort(key=lambda e: e[2])


def parse_args(cli_args: list[str]):
    """
        Use argparse module to parse the command line arguments.

        Returns:
            (Namespace object): Parsed arguments in argparse Namespace
                                object.
        """
    parser = argparse.ArgumentParser(description='Find bioisosteric linkers in'
                                                 ' ChEMBL database.')
    in_args = parser.add_mutually_exclusive_group(required=True)
    in_args.add_argument('-C', '--chembl-database-file', dest='chembl_file',
                         help='Name of SQLite ChEMBL database file.')
    in_args.add_argument('--input-series-file', dest='in_series',
                         help='Name of a JSON file containing series from'
                              ' previous run.')
    parser.add_argument('-O', '--output-db-file', dest='out_db_file',
                        required=True,
                        help='Name of SQLite database for results.')
    parser.add_argument('--output-series-file', dest='out_series',
                        help='Name of JSON file to write the series to.')
    parser.add_argument('--keep-herg', dest='keep_herg',
                        action='store_true',
                        help='By default, HERG assays are removed.  This'
                             ' retains them.')
    parser.add_argument('--keep-p450', dest='keep_p450',
                        action='store_true',
                        help='By default, P450 assays are removed.  This'
                             ' retains them.')
    parser.add_argument('--limit', dest='limit',
                        default=-1, type=IntRange(-1),
                        help='For testing, the limit on the number of assays'
                             ' in the initial query.  Default=-1 means no'
                             ' limit.')
    parser.add_argument('--max-heavies', dest='max_heavies',
                        type=IntRange(1), default=8,
                        help='Maximum number of heavy atoms in linker,'
                             ' excluding the 2 ring atoms being linked.'
                             '  Default=%(default)s.')
    parser.add_argument('--max-bonds', dest='max_bonds',
                        type=IntRange(2), default=5,
                        help='Maximum number of bonds in shortest path between'
                             ' ring atoms being linked.  Default=%(default)s.')
    parser.add_argument('--min-examples', dest='min_examples',
                        default=2, type=IntRange(1),
                        help='Minimum number of examples of a bioisostere'
                             ' in different documents/assays for it to be'
                             ' reported.  Default=%(default)s.')
    parser.add_argument('--num-procs', dest='num_procs',
                        default=1, type=IntRange(1, cpu_count() - 1),
                        help='Number of processors if running in parallel.')
    parser.add_argument('--image-dir', dest='image_dir',
                        default='BioisostereImages',
                        help='If writing an HTML file for the output, where to'
                             ' put the SVGs of the linkers.')
    parser.add_argument('--html-file', dest='html_file',
                        help='Name of file for HTML table of bioisosteric'
                             ' linkers.')
    args = parser.parse_args(cli_args)
    return args


def extract_herg_assays(dbfile: str) -> set[str]:
    """
    Pull the names of assays with TID 165, which is hERG,
    CHEMBL240
    Args:
        dbfile:

    Returns:

    """
    sql = """SELECT DISTINCT chembl_id FROM assays WHERE tid=165"""
    con = sqlite3.connect(dbfile)
    herg_assays = [assay[0] for assay in con.execute(sql)]
    # print(herg_assays)
    return set(herg_assays)


def extract_p450_assays(dbfile: str) -> set[str]:
    """
    Pull the names of assays with P450 in the name.
    Args:
        dbfile:

    Returns:

    """
    sql = """SELECT chembl_id FROM assays
     WHERE assays.tid IN
      (SELECT tid FROM target_dictionary
      WHERE pref_name LIKE '%P450%')"""
    con = sqlite3.connect(dbfile)
    p450_assays = [assay[0] for assay in con.execute(sql)]
    # print(herg_assays)
    return set(p450_assays)


def assay_id_to_chembl_id(dbfile: str) -> dict[str, str]:
    """
    The activities table contains assay_id values, which are internal
    to the database.  Users expect the ChEMBL ID, which is what you
    find in the web interface.  This returns the map of the former
    to the latter.
    Args:
        dbfile:

    Returns:

    """
    sql = 'SELECT assay_id, chembl_id FROM assays'
    con = sqlite3.connect(dbfile)
    assay_map = {}
    for row in con.execute(sql):
        assay_map[row[0]] = row[1]

    return assay_map


def extract_series(dbfile: str, limit: int) -> list[dict[str, Union[str, list[str]]]]:
    """
    Pull the series out of the ChEMBL database.  As defined by Ertl,
    a series is one where compounds have activities <10uM in the
    same assay_id and doc_id i.e. in the same assay as reported in the
    same publication and there are at least 3 of them.

    Args:
        dbfile:
        limit: number of assays considered in first query, for
               testing.  -1 means no limit.

    Returns:
        list of lists of SMILES strings, each inner list being a
        series
    """
    # someone more adept in SQL could probably do this all in one go.
    assay_map = assay_id_to_chembl_id(dbfile)

    sql = """SELECT doc_id, assay_id, value, canonical_smiles, molecule_dictionary.chembl_id
FROM activities, compound_structures, molecule_dictionary
WHERE activities.molregno = compound_structures.molregno
AND molecule_dictionary.molregno = activities.molregno
AND units == 'uM' and value < 10
ORDER BY doc_id, assay_id, molecule_dictionary.chembl_id
    """
    if limit != -1:
        sql += f'\nLIMIT {limit}'
    con = sqlite3.connect(dbfile)
    series = []
    last_doc = None
    last_assay = None
    for (doc_id, assay_id, acty, smi, chemblid) in con.execute(sql):
        # print(f'doc: {doc_id} assay: {assay_id} activity: {acty}'
        #       f' SMILES: {smi}  ChEMBLID: {chemblid}')
        if doc_id == last_doc and assay_id == last_assay:
            series[-1]['doc'] = doc_id
            series[-1]['assay'] = assay_map[assay_id]
            series[-1]['SMILES'].append(smi)
            series[-1]['ChEMBLID'].append(chemblid)
        else:
            last_doc = doc_id
            last_assay = assay_id
            series.append({'doc': str(doc_id),
                           'assay': assay_id,
                           'SMILES': [smi],
                           'ChEMBLID': [chemblid]})

    out_series = [ser for ser in series if len(ser['SMILES']) > 2]
    # print(f'initial assays : {sorted([ser["assay"] for ser in series])}')

    return out_series


def summarise_series(series: list[dict]) -> None:
    mean_size = 0
    for ser in series:
        # print(f"{ser['doc']}  {ser['assay']} {len(ser['SMILES'])}")
        mean_size += len(ser['SMILES'])
    num_sers = len(series)
    mean_size /= num_sers
    max = len(series[0]['SMILES'])
    min = len(series[-1]['SMILES'])
    median = len(series[num_sers // 2]['SMILES'])
    print(f'Number of series : {len(series)}  median size : {median}'
          f'  mean size : {mean_size:6.2f} max : {max} ({series[0]["assay"]})  min : {min}')
    # for i, ser in enumerate(series[:20]):
    #     print(f'Series {i} size : {len(ser["SMILES"])}')


def add_bioisosteres_if_different(new_bio: Bioisostere,
                                  bioisosteres: list[Bioisostere]) -> None:
    """
    If the bioisostere is new, add a new Bioisostere, otherwise add it
    to the existing one
    Args:
        new_bio:
        bioisosteres:

    Returns:

    """
    found_it = False
    for bios in bioisosteres:
        if bios == new_bio:
            bios.add_examples(new_bio.examples)
            found_it = True
            break

    if not found_it:
        bioisosteres.append(new_bio)


def make_series(db_file: str, keep_herg: bool, keep_p450: bool,
                limit: int) -> list[dict[str, Union[str, list[str]]]]:
    """
    Process the ChEMBL SQLite database to pull out series of compounds.
    As defined by Ertl, a series is one where compounds have
    activities <10uM in the same assay_id and doc_id i.e. in the same
    assay as reported in the same publication and there are at least
    3 of them.
    Args:
        db_file: Name of a ChEMBL SQLite database file.
        keep_herg: if True, leaves in the HERG assays
        keep_p450: if True, leaves in the P450 assays.
        limit: number of assays considered in first query, for
               testing.  -1 means no limit.

    Returns:
        a list of dicts of information about the series.
    """
    series = extract_series(db_file, limit)
    series.sort(key=lambda s: len(s['SMILES']), reverse=True)
    print('Initial series : ', end='')
    summarise_series(series)

    if keep_herg:
        series_no_herg = series
    else:
        herg_assays = extract_herg_assays(db_file)
        series_no_herg = [ser for ser in series if ser['assay'] not in herg_assays]
        print('After removing HERG results : ', end='')
        summarise_series(series_no_herg)

    if keep_p450:
        series_no_p450 = series_no_herg
    else:
        p450_assays = extract_p450_assays(db_file)
        series_no_p450 = [ser for ser in series_no_herg if ser['assay'] not in p450_assays]
        print('After removing P450 results : ', end='')
        summarise_series(series_no_p450)

    return series_no_p450


def find_bioisosteres_in_series(series: dict[str, Union[str, list[str]]],
                                max_heavies: int, max_bonds: int) -> list[Bioisostere]:
    """
    Process the series to find bioisosteres, which are linkers that
    come from molecules in the series where the left and right sides
    of the original molecules were the same.
    Args:
        series: the molecule series, as a dict that must include SMILES
                and ChEMBLID
        max_heavies:
        max_bonds:

    Returns:
            The list of biaisosteres
    """
    molecules = [(smi, name) for smi, name in zip(series['SMILES'], series['ChEMBLID'])]
    linkers = fl.find_all_linkers(molecules, max_heavies=max_heavies,
                                  max_length=max_bonds)
    linker_smis = list(linkers.keys())
    bioisosteres = []
    for i in range(len(linker_smis) - 1):
        # print(f'{i} of {len(linker_smis)} : {len(linkers[linker_smis[i]])}')
        for linker_i in linkers[linker_smis[i]]:
            for j in range(i + 1, len(linker_smis)):
                # print(f'  {i} ::  {j} of {len(linker_smis)} : {len(linkers[linker_smis[j]])}')
                for linker_j in linkers[linker_smis[j]]:
                    if linker_i == linker_j:
                        continue
                    linker_j_mutant = copy.copy(linker_j)
                    linker_j_mutant._linker = linker_i._linker
                    linker_j_mutant._linker_smiles = None
                    if linker_i == linker_j_mutant:
                        new_bio = Bioisostere(linker_i, linker_j, series['doc'])
                        add_bioisosteres_if_different(new_bio, bioisosteres)
    return bioisosteres


def write_series_json(series: list[dict[str, Union[str, list[str]]]], json_file: str) -> bool:
    """
    Write the series into the given JSON file. Returns bool on success.
    """
    try:
        with open(json_file, 'w') as f:
            f.write(json.dumps(series, indent=4))
            f.write('\n')
    except FileNotFoundError:
        print(f'ERROR : failed to write file {json_file}.')
        return False

    return True


def read_series_json(json_file: str) -> Optional[list[dict[str, Union[str, list[str]]]]]:
    """
    Read the series from the given JSON file. Returns None if file not
    found.
    """
    try:
        with open(json_file, 'r') as f:
            series = json.loads(f.read())
    except FileNotFoundError:
        print(f'ERROR : failed to read file {json_file}.')
        return None

    return series


def extract_bioisosteres(series: list[dict[str, Union[str, list[str]]]],
                         max_heavies: int, max_bonds: int,
                         min_examples: int = 2,
                         num_procs: int = 1) -> list[Bioisostere]:
    """
    Pull out all the bioisosteres from the series passed in.
    Args:
        series:
        max_heavies:
        max_bonds:
        min_examples:
        num_procs

    Returns:

    """

    def merge_bioisosteres(new_bios, old_bios):
        if old_bios is None:
            old_bios = new_bios
        else:
            for nb in new_bios:
                add_bioisosteres_if_different(nb, old_bios)
        return old_bios

    if num_procs == 1:
        bioisosteres = None
        for ser in series:
            next_bios = \
                find_bioisosteres_in_series(ser, max_heavies=max_heavies,
                                            max_bonds=max_bonds)
            bioisosteres = merge_bioisosteres(next_bios, bioisosteres)
    else:
        bioisosteres = None
        with cf.ProcessPoolExecutor(max_workers=num_procs) as pool:
            for next_bios in pool.map(find_bioisosteres_in_series,
                                      series,
                                      itertools.repeat(max_heavies),
                                      itertools.repeat(max_bonds),
                                      chunksize=1):
                bioisosteres = merge_bioisosteres(next_bios, bioisosteres)

    trimmed_bios = []
    max_ex = -1
    for bios in bioisosteres:
        if bios.num_docs >= min_examples:
            trimmed_bios.append(bios)
        if bios.num_docs > max_ex:
            max_ex = bios.num_docs

    if not trimmed_bios:
        print(f'No bioisosteres with at least {min_examples} examples in'
              f' different documents/assays.  The maximum was {max_ex}.')
    else:
        trimmed_bios.sort(key=lambda b: (b.num_docs, b.num_examples), reverse=True)
    return trimmed_bios


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
                    <th>Index</th>
                    <th>Linker 1</th>
                    <th>Example 1</th>
                    <th>Linker 2</th>
                    <th>Example 2</th>
                    <th>Doc ID</th>
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


def make_image(mol_smiles: str, image_dir: Union[str, Path],
               mol_name: Optional[str] = None,
               linker_atoms: Optional[list[int]] = None) -> str:
    """
    Make an SVG of the molecule in image_dir, named by the mol name
    if given, or a random name if not, and returns name of file.
    image_dir is created if it doesn't exist.  If mol_name is given,
    it is used as a legend in the image, and if linker_smiles is given,
    it will be used to highlight the linker atoms.  To do that, the
    atom maps are changed to ring atoms.

    Args:
        mol_smiles:
        image_dir:
        mol_name:
        linker_atoms:

    Returns:
        name of image file created.
    """
    image_path = Path(image_dir)
    if not image_path.exists():
        image_path.mkdir()

    mol = Chem.MolFromSmiles(mol_smiles)
    if mol_name is None:
        fh, img_file = mkstemp(suffix='.svg', dir=image_dir)
        close(fh)
        mol_name = ''
    else:
        img_file = image_path / f'{mol_name}.svg'

    if linker_atoms is not None:
        match_bonds = []
        for i in range(len(linker_atoms) - 1):
            for j in range(i + 1, len(linker_atoms)):
                bond = mol.GetBondBetweenAtoms(linker_atoms[i], linker_atoms[j])
                if bond is not None:
                    match_bonds.append(bond.GetIdx())
    else:
        match_bonds = None

    image_width = -1
    image_height = -1
    drawer = rdMolDraw2D.MolDraw2DSVG(image_width, image_height)
    rdDepictor.Compute2DCoords(mol)
    rdMolDraw2D.DrawMoleculeACS1996(drawer, mol, legend=mol_name,
                                    highlightAtoms=linker_atoms,
                                    highlightBonds=match_bonds)
    drawer.FinishDrawing()
    with open(img_file, 'w') as f:
        f.write(drawer.GetDrawingText())

    return img_file


def write_bioisostere_table(biososteres: list[Bioisostere],
                            html_file: str = 'bioisosteres.html',
                            images_dir: str = 'BioisostereImages') -> None:
    htmlf = start_html_file(html_file)
    for i, bios in enumerate(biososteres):
        l1img = make_image(bios.linker1_smiles, images_dir)
        l2img = make_image(bios.linker2_smiles, images_dir)
        for ex in bios.examples:
            l1eximg = make_image(ex[0].mol_smiles, images_dir, ex[0].name)
            l2eximg = make_image(ex[1].mol_smiles, images_dir, ex[1].name)
            htmlf.write('<tr>')
            htmlf.write(f'<td>{i} : {bios.num_docs}</td>')
            htmlf.write(f'<td><img src={l1img} /></td>')
            htmlf.write(f'<td><img src={l1eximg}/></td>')
            htmlf.write(f'<td><img src={l2img}/></td>')
            htmlf.write(f'<td><img src={l2eximg}/></td>')
            htmlf.write(f'<td>{ex[2]}</td>')
            htmlf.write(f'</tr>')

    finish_html_file(htmlf)


def create_database_tables(conn: sqlite3.Connection) -> None:
    """
    Set up the output database tables.
    The bioisosteres table is the key one, holding the bioisosteres as
    pairs of linker SMILES, as well as counts of where they came from.
    The linkers table holds all the linkers used in bioisosteres,
    and refers back to it using bioisostere_id, so you can find the
    two linkers that the bioisostere refers to by searching this
    by bioisostere_id.  There are multiple linker records wih the
    same linker_smiles so that you can find out the environments of
    the linkers in the molecules they came from.
    """
    curs = conn.cursor()
    sql = """CREATE TABLE bioisosteres(
    bioisostere_id INTEGER PRIMARY KEY,
    linker1_smiles TEXT NON NULL,
    linker2_smiles TEXT NON NULL,
    num_docs INTEGER,
    num_examples INTEGER)
    """
    curs.execute(sql)

    sql = """CREATE TABLE linkers(
    linker_id INTEGER PRIMARY KEY,
    name TEXT NON NULL,
    linker_smiles TEXT NON NULL,
    linker_atoms TEXT NON NULL,
    mol_smiles TEXT NON NULL,
    left_smiles TEXT NON NULL,
    right_smiles TEXT NON NULL,
    path_length INTEGER,
    num_donors INTEGER,
    num_acceptors INTEGER,
    bioisostere_id INTEGER NOT NULL,
    FOREIGN KEY (bioisostere_id) REFERENCES bioisosteres(bioisostere_id)
    )
    """
    curs.execute(sql)

    conn.commit()


def add_bioisostere_to_table(conn: sqlite3.Connection, bios: Bioisostere) -> int:
    """
    Add just the biosisostere details itself to the biosisostere table.
    Returns its ID.
    """
    sql = """INSERT INTO bioisosteres
    (linker1_smiles, linker2_smiles, num_docs, num_examples)
    VALUES (?,?,?,?)
    """
    curs = conn.cursor()
    curs.execute(sql, (bios.linker1_smiles, bios.linker2_smiles,
                       bios.num_docs, bios.num_examples))
    conn.commit()
    curs.execute('SELECT last_insert_rowid()')
    id = int(curs.fetchone()[0])
    return id


def add_linker_to_table(curs: sqlite3.Cursor, linker: fl.Linker,
                        bios_id: int) -> int:
    """
    Add the linker to the linkers table, returning its primary key.
    """
    sql = """INSERT INTO linkers
    (name, linker_smiles, linker_atoms, mol_smiles, left_smiles, right_smiles,
     path_length, num_donors, num_acceptors, bioisostere_id)
    VALUES (?,?,?,?,?,?,?,?,?,?)
    """
    linker_atoms = ' '.join([str(la) for la in linker.linker_atoms])
    curs.execute(sql, (linker.name, linker.linker_smiles, linker_atoms,
                       linker.mol_smiles, linker.left_smiles,
                       linker.right_smiles, linker.path_length,
                       linker.num_donors, linker.num_acceptors,
                       bios_id))
    curs.execute('SELECT last_insert_rowid()')
    id = int(curs.fetchone()[0])
    return id


def add_linkers_to_table(conn: sqlite3.Connection, bios: Bioisostere,
                         bios_id: int) -> None:
    """
    Add the examples from the bioisostere, primary key bios_id, to the
    linkers and examples tables.
    """
    curs = conn.cursor()
    for ex in bios.examples:
        _ = add_linker_to_table(curs, ex[0], bios_id)
        _ = add_linker_to_table(curs, ex[1], bios_id)
    conn.commit()


def write_bioisosteres_db(bioisosteres: list[Bioisostere], db_file: str) -> None:
    """
    Write the bioisosteres to a SQLite database.
    Args:
        bioisosteres:
        db_file:

    Returns:

    """
    db_path = Path(db_file)
    if db_path.exists():
        print(f'WARNING : existing file {db_file} will be overwritten.')
        db_path.unlink()
    conn = sqlite3.connect(db_file)
    create_database_tables(conn)
    for bios in bioisosteres:
        bios_id = add_bioisostere_to_table(conn, bios)
        add_linkers_to_table(conn, bios, bios_id)
    conn.close()


def main(cli_args):
    print(f'Using RDKit version {rdBase.rdkitVersion}.')

    args = parse_args(cli_args)
    if args is None:
        return False

    if args.chembl_file is not None:
        series = make_series(args.chembl_file, args.keep_herg, args.keep_p450,
                             args.limit)
    else:
        series = read_series_json(args.in_series)
    if series is None:
        return False
    summarise_series(series)

    if args.out_series is not None:
        write_series_json(series, args.out_series)

    bioisosteres = extract_bioisosteres(series, args.max_heavies,
                                        args.max_bonds, args.min_examples,
                                        args.num_procs)

    print(f'Number of bioisosteres : {len(bioisosteres)}')
    for bios in bioisosteres:
        print(bios)
        for ex in bios.examples:
            print(f'  {ex[0].name} {ex[0].mol_smiles} ::'
                  f' {ex[1].name} {ex[1].mol_smiles}')

    if args.html_file is not None:
        write_bioisostere_table(bioisosteres, html_file=args.html_file,
                                images_dir=args.image_dir)
    write_bioisosteres_db(bioisosteres, args.out_db_file)

    return True


if __name__ == '__main__':
    sys.exit(not main(sys.argv[1:]))
