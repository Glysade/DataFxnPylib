"""
===========
bio_util.py
===========

Copyright (C) 2017-2022 Glysade, LLC

A number of utility functions to augment Biopython functionality
"""

import os
from enum import Enum
from typing import List, Optional

from Bio.Seq import Seq, UndefinedSequenceError
from Bio import Entrez, SeqIO
from Bio.SeqFeature import SeqFeature
from Bio.SeqRecord import SeqRecord

from ruse.util.data_table import ColumnInfo, ColumnDataType

# IUPAC ambiguous RNA/DNA with gap
DNA_ALPHABET = 'GATUCRYWSMKHBVDN-'
# IUPAC extended protein with gap and stop (no 3 letter codes)
# Also incomplete codon '_'
PROTEIN_ALPHABET = 'ACDEFGHIKLMNPQRSTVWYBXZJUO-*_'


class SequenceType(Enum):
    PROTEIN = 1
    DNA = 2


def is_dna(seq: str) -> bool:
    """
    Return true if seq could be a dna sequence.  Note that there are many protein sequences that could be dna sequences

    :param seq: the sequence to test
    :return: True if the sequence contains all DNA alphabet
    """

    for b in seq:
        if b.upper() not in DNA_ALPHABET:
            return False
    return True


def is_dna_record(rec: SeqRecord):
    if 'molecule_type' in rec.annotations:
        molecule_type = rec.annotations['molecule_type'].upper()
        if 'DNA' in molecule_type:
            return True
        if 'RNA' in molecule_type:
            return True
        return False
    return is_dna(rec.seq)


def is_protein(seq: str) -> bool:
    """
    Return true if seq could be a protein sequence, note any dna sequence could also be a protein sequence, so
    not is_dna may be a better test

    :param seq: The sequence to test
    :return: True if the sequence contains all Protein alphabet
    """

    for b in seq:
        if b.upper() not in PROTEIN_ALPHABET:
            return False
    return True


def entrez_email() -> str:
    """
    Returns the email to be passed to NCBI when retrieving information from NCBI using the Entrez interface.
    This is set using the ENTREZ_EMAIL environment variable

    :return: Email to pass to NCBI
    """
    if 'ENTREZ_EMAIL' not in os.environ:
        raise ValueError("Environment variable ENTREZ_EMAIL not set.")
    else:
        return os.environ['ENTREZ_EMAIL']


def extract_feature(record: SeqRecord, feature_id: str = None, qualifier_name: str = None, qualifier_value: str = None,
                    type: str = None) -> SeqFeature:
    """
    Extracts a feature from a sequence record
    Can search by type, feature_id or both qualifier_name and qualifier_value
    Matching is case insensitive

    :param record: :class'`SeqRecord` object to extract feature from
    :param feature_id: feature id string criteria
    :param qualifier_name: qualifier name criteria
    :param qualifier_value: qualifier value criteria
    :param type: feature type restriction
    :return: matching :class:`SeqFeature`
    """
    if not record:
        return None

    if qualifier_name or qualifier_value:
        assert qualifier_value, "qualifier_name set but no value given"
        assert qualifier_name, "qualifier_value set but no name given"

    assert qualifier_name or feature_id or type, "no criteria given"

    def match_feature(feature: SeqFeature) -> bool:
        if not feature.location:
            return False
        if type:
            if feature.type.lower() != type.lower():
                return False
        if feature_id:
            if feature.id.lower() != feature_id.lower():
                return False
        if qualifier_name:
            lookup = {k.lower(): v for k, v in feature.qualifiers.items()}
            if qualifier_name.lower() not in lookup:
                return False
            if qualifier_value.lower() not in [n.lower() for n in lookup[qualifier_name.lower()]]:
                return False
        return True

    feature = next((f for f in record.features if match_feature(f)), None)
    if feature:
        out_record = feature.extract(record)
        out_record.annotations = record.annotations.copy()
        return out_record
    return None


def sub_sequence(record: SeqRecord, start: int, end: int) -> Optional[SeqRecord]:
    """
    Extract a sequence range from a :class:`SeqRecord` while copying annotations

    :param record: sequence to extract range from
    :param start: start of sequence range
    :param end: end of sequence range
    :return: sequence range
    """
    assert end > start
    if not record:
        return None
    sub_record = record[start:end]
    sub_record.annotations = record.annotations.copy()
    return sub_record


def sequences_to_file(file: str, sequences: List[SeqRecord], format: str = 'fasta') -> None:
    """
    Writes sequences to a file

    :param file:
    :param sequences:
    :param format:
    """
    with open(file, 'w') as fh:
        SeqIO.write(sequences, fh, format)


def is_genbank_column(column_def: ColumnDataType) -> bool:
    if column_def['dataType'] == 'binary' and column_def['properties']['ContentType'] == 'chemical/x-genbank':
        return True
    else:
        return False


def add_annotation_comment(record: SeqRecord, comment: str) -> None:
    if 'comment' in record.annotations:
        current_comment = record.annotations['comment']
        if isinstance(current_comment, str):
            record.annotations['comment'] = [current_comment, comment]
        elif isinstance(current_comment, (tuple, list)):
            comments = list(current_comment)
            comments.append(comment)
            record.annotations = comments
        else:
            raise ValueError("Unknown type for seqrecord comment annotation")
    else:
        record.annotations['comment'] = [comment]


def ok_sequence(rec: SeqRecord) -> bool:
    if not rec or not rec.seq:
        return False
    s = str(rec.seq)
    if not is_dna(s) and not is_protein(s):
        return False
    not_all_allowed = set('_-')
    if set(s) <= not_all_allowed:
        return False
    return True


def is_defined_sequence(seq: Seq) -> bool:
    try:
        bytes(seq)
    except UndefinedSequenceError:
        return False
    return True


Entrez.email = entrez_email()
