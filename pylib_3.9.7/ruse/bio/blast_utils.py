"""
==============
blast-utils.py
==============

Copyright (C) 2017-2022 Glysade, LLC

A set of utility functions for handing NCBI Blast search and Entrez queries
"""

import os
import re
import subprocess
import uuid
from collections import namedtuple
from glob import glob
from typing import List, Dict

from Bio import SeqIO, Entrez
from Bio.Data.IUPACData import ambiguous_dna_values
from Bio.SeqRecord import SeqRecord

from ruse.bio.bio_util import entrez_email
from ruse.bio.blast_search import full_path_to_blast_exe
from ruse.bio.sequence_retrieval import emboss_database_type, retrieve_emboss_records
from ruse.util.frozen import Frozen

Entrez.email = entrez_email()


def truncate_full_length(seq: SeqRecord) -> SeqRecord:
    maximum_length = 10_000
    if len(seq) <= maximum_length:
        return seq
    else:
        return seq[0: maximum_length]


def retrieve_entrez_records(db: str, target_ids: List[str]) -> List[SeqRecord]:
    """
    Uses :func:`Bio.Entrez.efetch` to retrieve reference sequences from NCBI

    :param db: Etnrez database ("nt" or "nr")
    :param target_ids: List of Entrez/NCBI identifiers
    :return:  a List of Genbank entries as :class:`Bio.SeqRecord.SeqRecord` objects
    """

    with Entrez.efetch(db=db, rettype="gb", retmode="text",
                       id=target_ids) as handle:
        records = []
        try:
            for entry in SeqIO.parse(handle, 'gb'):
                records.append(entry)
        except TypeError:
            if len(records) != len(target_ids):
                raise ValueError('Unable to retrieve all records from Entrez');
        # this line that I used prior to Biopython 1.79 is throwing a type error
        # Interestingly the new code above doesn't do that.
        # records = list(SeqIO.parse(handle, "gb"))
    return records


def blast_databases_list() -> List[str]:
    dir = os.environ['BLASTDB']
    assert dir, "BLASTDB environment variable not set!"
    files = glob(os.path.join(dir, '*'))
    files = [os.path.basename(f) for f in files]
    files = [f.split('.')[0] for f in files]
    dbs = list(set(files))
    return dbs

def _create_dna_equivalence_lookup() -> Dict[str, Dict[str, bool]]:
    """
    Create an equivalence lookup for ambiguous dna bases

    :return: a Dictionary of Dictionaries that is a lookup to determine if two ambiguous dna bases are equivalent
    """
    lookup = {}
    for n1, n1matches in ambiguous_dna_values.items():
        for n2, n2matches in ambiguous_dna_values.items():
            if n1 not in lookup:
                lookup[n1] = {}
            if set(n1matches) & set(n2matches):
                lookup[n1][n2] = True
            else:
                lookup[n1][n2] = False
    return lookup


class SequenceMatcher(object):
    """
    A class to facilitate matching ambiguous dna sequences
    """
    dna_equivalence = _create_dna_equivalence_lookup()

    @classmethod
    def match_sequences(cls, seq: str, sub: str) -> int:
        """
        A class method to determine the the first position that an ambiguous dna subsequence can match a dna sequence
        :param seq: The main sequence
        :param sub: Subsequence
        :return: The position of start of the first match of the subsequence to the sequence, or -1 if no such match exists
        """
        seq = seq.replace('-', '')
        sub = sub.replace('-', '')

        for seq_pos in range(0, len(seq)-len(sub)+1):
            if cls._match_portion(seq_pos, seq, sub):
                return seq_pos
        return -1

    @classmethod
    def _match_portion(cls, start: int, seq: str, sub: str) -> bool:
        for sub_pos in range(0, len(sub)):
            if not cls.dna_equivalence[seq[start+sub_pos]][sub[sub_pos]]:
                return False
        return True


class BlastRecord(namedtuple('BlastRecord', ['id', 'record'])):
    """
    Stores a local blast database sequence record and the id used to retrieve it

    Attributes:
        * id: identifier that was used to extract the local record
        * record: sequence as :class:`Bio.SeqRecord.SeqRecord`
    """


class BlastRecords(Frozen):
    """
    Class to facilitate retrieval of sequence records from a local database
    Uses the external program *blastdbcmd* to export records from the local database

    Attributes:
        * ids: list of ids retrieved from the local database
        * records: retrieved records as list of :class:`Bio.SeqRecord.SeqRecord`
        * record_map: Dictionary mapping ids to records
        * database: The local database
        * error: Any blastdbcmd error
    """

    def __init__(self):
        """
        Empty constructor
        """
        self.ids = []  # type: List[str]
        self.records = []  # type: List[BlastRecord]
        self.record_map = {}  # type: Dict[str, BlastRecord]
        self.database = ""
        self.error = ""

    def retrieve_from_local_database(self, database: str, ids: List[str]) -> None:
        if emboss_database_type(database):
            sequences = retrieve_emboss_records(database, ids)
            records = [BlastRecord(sequence.name, sequence) for sequence in sequences]
        else:
            records = self.retrieve_from_local_blast_database(database, ids)
        self.records = records
        self.record_map = {record.id: record for record in self.records}
        for id in ids:
            assert id in self.record_map
        self.ids = ids

    def retrieve_from_local_blast_database(self, database: str, ids: List[str]) -> List[BlastRecord]:
        """
        Retrieves a list of target sequences from a local blast database

        :param database: The name of the local database
        :param ids: The ids to retrieve
        """

        self.database = database
        file_id = str(uuid.uuid4())
        file_in = "{}_in.txt".format(file_id)
        file_out = "{}_out.fasta".format(file_id)
        with open(file_in, 'w') as fh:
            for id in ids:
                fh.write("{}\n".format(id))

        blastdbcmd = full_path_to_blast_exe('blastdbcmd')
        args = [blastdbcmd, '-db', self.database, '-entry_batch', file_in, '-outfmt', '%f', '-out', file_out]
        try:
            subprocess.run(args, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        except subprocess.CalledProcessError as ex:
            self.error = "Blastdbcmd error {}".format(ex)
            return
        else:
            with open(file_out) as fh:
                sequences = list(SeqIO.parse(fh, 'fasta'))
        finally:
            os.remove(file_out)
            os.remove(file_in)

        assert len(sequences) == len(ids)
        records = [BlastRecord(id, s) for id, s in zip(ids, sequences)]
        return records

