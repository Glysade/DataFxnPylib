import itertools
import os
import re
import uuid
from io import TextIOWrapper
from typing import List, Dict, Optional, Tuple

from Bio import SeqIO
from Bio.Application import ApplicationError
from Bio.SeqRecord import SeqRecord

from ruse.bio.applications import SeqretFeatureCommandline, ShowDbCommandLine
from ruse.bio.bio_util import SequenceType, add_annotation_comment


class SequenceRetrieval:
    usa_db_matcher = re.compile(r"^([\w\_]+)\:")


class SwissIOWrapper(TextIOWrapper):
    def __init__(self, handle: TextIOWrapper):
        super().__init__(handle)
        self.handle = handle
        self.read_pe = False

    def readline(self, size=-1) -> str:
        line = self.handle.readline(size)
        if self.read_pe and line.startswith(';'):
            self.read_pe = False
            return self.readline(size)
        self.read_pe = line.startswith('PE')
        return line


def emboss_database_type(database_name: str) -> Optional[SequenceType]:
    databases = _available_databases()
    if not databases:
        return None
    if database_name in databases['Protein']:
        return SequenceType.PROTEIN
    if database_name in databases['Nucleotide']:
        return SequenceType.DNA
    return None


def retrieve_emboss_records(database: str, ids: List[str]) -> List[SeqRecord]:
    usas = ["{}:{}".format(database, id) for id in ids]
    sequences = retrieve_sequences(usas)
    assert len(sequences) == len(usas)
    return sequences


def _available_databases() -> Optional[Dict[str, List[str]]]:
    if not os.getenv('EMBOSSRC'):
        return None;
    file = '/tmp/showdb.txt'
    read_file = False
    if os.path.exists(file):
        mtime = os.path.getmtime(file)
        embossrc = os.path.join(os.getenv('EMBOSSRC'), ".embossrc")
        emboss_mtime = os.path.getmtime(embossrc)
        if emboss_mtime > mtime:
            read_file = True
    if not os.path.exists(file) or read_file:
        showdb_args = {'protein': True, 'nucleic': True, 'type': True, 'outfile': file}
        cline = ShowDbCommandLine(**showdb_args)
        cline()

    databases = {'Protein': [], 'Nucleotide': []}
    matcher = re.compile(r"^([\w\_]+)\s+(\w+)")
    with open(file, 'r') as fh:
        for line in fh.readlines():
            if not line.startswith('#'):
                m = matcher.match(line)
                if m:
                    db = m.group(1)
                    t = m.group(2)
                    databases.setdefault(t, []).append(db)

    return databases


def check_usas(usas: List[str], sequence_type: Optional[SequenceType] = None) -> Tuple[
    Optional[SequenceType], List[str]]:
    databases = _available_databases()
    if not databases:
        return None, []
    protein_databases = databases['Protein']
    dna_databases = databases['Nucleotide']
    expanded_usas = []
    for usa in usas:
        m = SequenceRetrieval.usa_db_matcher.match(usa)
        if not m:
            if not sequence_type:
                raise ValueError("No database specified in {} and sequence type missing".format(usa))
            if sequence_type == SequenceType.PROTEIN:
                queries = ['{}:{}'.format(db, usa) for db in protein_databases]
                expanded_usas.extend(queries)
            else:
                assert sequence_type == SequenceType.DNA
                queries = ['{}:{}'.format(db, usa) for db in dna_databases]
                expanded_usas.extend(queries)
        else:
            expanded_usas.append(usa)
            db = m.group(1)
            if db not in protein_databases and db not in dna_databases:
                raise ValueError("Unknown database {}".format(db))
            if not sequence_type:
                if db in protein_databases:
                    sequence_type = SequenceType.PROTEIN
                if db in dna_databases:
                    sequence_type = SequenceType.DNA
            else:
                if db in dna_databases and sequence_type == SequenceType.PROTEIN:
                    raise ValueError("Mixed protein and dna databases USA queries")
                if db in protein_databases and sequence_type == SequenceType.DNA:
                    raise ValueError("Mixed protein and dna databases USA queries")
    return sequence_type, expanded_usas


def retrieve_sequences(usas: List[str], max_sequences: int = 10000, sequence_type: Optional[SequenceType] = None) -> \
        Optional[list]:
    (sequence_type, usas) = check_usas(usas, sequence_type)
    if not sequence_type:
        return None
    basename = str(uuid.uuid4())
    if len(usas) > 1:
        list_file = '{}_usa.list'.format(basename)
        with open(list_file, 'w') as fh:
            fh.writelines(usa + '\n' for usa in usas)
        query = 'list:{}'.format(list_file)
        database = "unknown"
    else:
        list_file = None
        query = usas[0]
        m = SequenceRetrieval.usa_db_matcher.match(query)
        assert m
        database = m.group(1)

    if sequence_type == SequenceType.DNA:
        emboss_format = "embl"
        bio_format = "embl"
    else:
        emboss_format = "swiss"
        bio_format = "swiss"

    out_file = "{}_out.seqs".format(basename)
    args = {'sequence': query, 'outseq': out_file, 'osformat': emboss_format, 'feature': True}
    cline = SeqretFeatureCommandline(**args)

    # print("emboss command line is {}".format(cline))
    try:
        stdout, stderr = cline()
    except ApplicationError as e:
        # This is the message when the query has no match against the database- I don't really think this is an error
        # as such
        if "Bad value for '-sequence'" in e.stderr:
            return []
        raise e

    if list_file:
        os.remove(list_file)

    with open(out_file, 'r') as fh:
        swiss_fh = SwissIOWrapper(fh)
        seq_iter = itertools.islice(SeqIO.parse(swiss_fh, bio_format), max_sequences)
        sequences = list(seq_iter)
    for sequence in sequences:
        add_annotation_comment(sequence, "Seqret database: {}".format(database))
    os.remove(out_file)
    return sequences
