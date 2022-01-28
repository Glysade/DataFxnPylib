import os
import uuid
from collections import deque
from typing import List, Optional, NamedTuple

from Bio import SeqIO
from Bio.GenBank.Scanner import GenBankScanner
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord

from ruse.bio.applications import RestrictCommandLine


class RestrictionSite(NamedTuple):
    seq_name: str
    start: int
    end: int
    strand: str
    enzyme_name: str
    restriction_site: str
    prime5: int
    prime3: int
    frag5: Optional[int]
    frag3: Optional[int]
    primerev5: Optional[int]
    primerev3: Optional[int]
    fragrev5: Optional[int]
    fragrev3: Optional[int]


def enzyme_cleavage_sites(sequences: List[SeqRecord], enzymes: Optional[List[str]] = None, sitelen: int = 4) -> List[
    SeqRecord]:

    def line_to_site(line):
        args = line.strip().split('\t')
        if args[0] == 'SeqName':
            return None
        frag5 = int(args[8]) if args[8] != '.' else None
        frag3 = int(args[9]) if args[9] != '.' else None
        primerev5 = int(args[10]) if args[10] != '.' else None
        primerev3 = int(args[11]) if args[11] != '.' else None
        fragrev5 = int(args[12]) if args[12] != '.' else None
        fragrev3 = int(args[13]) if args[13] != '.' else None
        return RestrictionSite(args[0], int(args[1]), int(args[2]), args[3], args[4], args[5], args[6], args[7],
                               frag5=frag5, frag3=frag3, primerev5=primerev5, primerev3=primerev3, fragrev5=fragrev5,
                               fragrev3=fragrev3)

    def annotated_seq(in_seq: SeqRecord, site_iter):
        id = in_seq.id
        features = []
        while len(site_iter) > 0:
            if site_iter[0] is None:
                site_iter.popleft()
                continue
            if site_iter[0].seq_name != id:
                break
            site = site_iter.popleft()
            strand = -1 if site.strand == '-' else +1
            notes = [
                "enzyme_name: {}".format(site.enzyme_name),
                "restriction_site: {}".format(site.restriction_site),
                "5prime: {}".format(site.prime5),
                "3prime: {}".format(site.prime3)
            ]
            if site.frag5:
                notes.append("5frag: {}".format(site.frag5))
            if site.frag3:
                notes.append("3frag: {}".format(site.frag3))
            if site.primerev5:
                notes.append("5primerev: {}".format(site.primerev5))
            if site.primerev3:
                notes.append("3primerev: {}".format(site.primerev3))
            if site.fragrev5:
                notes.append("5fragrev: {}".format(site.fragrev5))
            if site.fragrev3:
                notes.append("3fragrev: {}".format(site.fragrev3))
            site_feature = SeqFeature(FeatureLocation(site.start - 1, site.end, strand=strand), type='misc_feature',
                                      qualifiers={'note': notes})
            features.append(site_feature)
        return SeqRecord(in_seq.seq, id=id, name=in_seq.name, description=in_seq.description, features=features,
                         annotations={'molecule_type': 'DNA'})

    basename = str(uuid.uuid4())
    in_file = '{}_restrict_in.fasta'.format(basename)
    out_file = '{}_restrict_out.txt'.format(basename)
    with open(in_file, 'w') as fh:
        SeqIO.write(sequences, fh, 'fasta')

    if enzymes:
        enzymes_arg = ','.join(enzymes)
    else:
        enzymes_arg = 'all'
    cline = RestrictCommandLine(sequence=in_file, enzymes=enzymes_arg, sitelen=str(sitelen), rformat='excel',
                                outfile=out_file)
    print("restrict command line is {}".format(cline))
    stdout, stderr = cline()

    out_sequences = []
    with open(out_file, 'r') as fh:
        lines = fh.readlines()
        site_iter = deque((line_to_site(l) for l in lines))
        for in_seq in sequences:
            out_seq = annotated_seq(in_seq, site_iter)
            out_sequences.append(out_seq)

    if os.path.exists(in_file):
        os.remove(in_file)
    if os.path.exists(out_file):
        os.remove(out_file)

    return out_sequences

