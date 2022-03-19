import os.path
import uuid
from typing import List, Optional

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from df.bio_helper import query_from_request, sequences_to_column
from df.data_transfer import DataFunction, DataFunctionRequest, DataFunctionResponse, string_input_field, ColumnData, \
    TableData, DataType
from ruse.bio.antibody import numbering_and_regions_from_sequence, ANTIBODY_NUMBERING_COLUMN_PROPERTY
from ruse.bio.applications import IgblastpCommandLine
from ruse.bio.blast_parse import build_common_alignments
from ruse.bio.igblast_parse import IgResult, parse_igblastp_results, annotate_igblast_regions
from ruse.bio.sequence_align import copy_features_to_gapped_sequence


def process_ig_search_results(query_sequence: SeqRecord, results: IgResult, table_name: str) -> DataFunctionResponse:
    ungapped = SeqRecord(query_sequence.seq.ungap(), query_sequence.id, query_sequence.name, query_sequence.description)
    alignment: List[SeqRecord] = build_common_alignments(ungapped, results.hits)
    annotated_query = annotate_igblast_regions(results, query_sequence)
    copy_features_to_gapped_sequence(annotated_query, alignment[0])

    pairs: List[Optional[str]] = [None]
    target_ids: List[Optional[str]] = [None]
    evalues: List[Optional[float]] = [None]
    bits: List[Optional[float]] = [None]
    chain: List[Optional[str]] = [None]
    for hit in results.hits:
        pairs.append(f'{hit.query}|{hit.target}')
        target_ids.append(hit.target_id)
        evalues.append(hit.evalue)
        bits.append(hit.bits)
        chain.append(hit.chain)

    alignment_column = sequences_to_column(alignment, 'Aligned Sequence', True)
    antibody_numbering = numbering_and_regions_from_sequence(alignment[0])
    if antibody_numbering:
        alignment_column.properties[ANTIBODY_NUMBERING_COLUMN_PROPERTY] = antibody_numbering.to_column_json()
    pair_column = ColumnData(name='Aligned Sequence Pairs', dataType=DataType.STRING,
                             contentType='chemical/x-sequence-pair', values=pairs)
    chain_column = ColumnData(name='Chain', dataType=DataType.STRING, values=chain)
    target_id_column = ColumnData(name='Target Id', dataType=DataType.STRING, values=target_ids)
    e_value_column = ColumnData(name='EValue', dataType=DataType.FLOAT, values=evalues)
    bit_column = ColumnData(name='Bits', dataType=DataType.FLOAT, values=bits)
    columns = [alignment_column, pair_column, chain_column, target_id_column, e_value_column, bit_column]

    output_table = TableData(tableName=table_name, columns=columns)
    response = DataFunctionResponse(outputTables=[output_table])
    return response


class IgBlastpLocalTextSearch(DataFunction):

    def execute(self, request: DataFunctionRequest) -> DataFunctionResponse:
        query_sequence = query_from_request(request)
        ig_data = string_input_field(request, 'igBlastDataDirectory')
        os.environ['IGDATA'] = ig_data

        id = str(uuid.uuid4())
        v_database = string_input_field(request, 'germlineDbV')
        v_database = os.path.join(os.environ['IGDATA'], 'database', v_database)
        organism = string_input_field(request, 'organism')
        domain_system = string_input_field(request, 'domainSystem')
        format = '"7 std qseq sseq btop"'

        query_file = f'query_{id}.fasta'
        with open(query_file, 'w') as fh:
            ungapped = SeqRecord(query_sequence.seq.ungap(), query_sequence.id)
            SeqIO.write([ungapped], fh, 'fasta')

        out_file = f'igblastp_{id}.out'
        command = IgblastpCommandLine(germline_db_V=v_database, organism=organism, domain_system=domain_system,
                                      outfmt=format, query=query_file, out=out_file)
        stdout, stderr = command()

        if stderr or stdout:
            raise RuntimeError(f'IGBLASTP failed stderr: {stderr} stdout: {stdout}')

        results: IgResult = parse_igblastp_results(out_file)
        response = process_ig_search_results(query_sequence, results, 'Igblastp Results')
        return response
