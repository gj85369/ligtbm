import json
import click
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastpCommandline
from Bio.Application import ApplicationError


class BLASTError(Exception):
    pass


def parse_titles(titles_line):
    pdbid_descr_to_chains = {}
    unparsed_titles = []
    for title in titles_line.split('>'):
        title = title.strip()
        if len(title.split('|')) != 5:
            unparsed_titles.append(title)
            continue
        gi, jobid, pdb, pdbid, chain_descr = title.split('|')
        if gi =='gi' and jobid.isnumeric() and pdb == 'pdb' and len(pdbid) == 4:
            chain = chain_descr[0]
            descr = chain_descr.replace('{chain} Chain {chain}, '.format(chain=chain), '')
            pdbid_descr = (pdbid, descr)
            if pdbid_descr not in pdbid_descr_to_chains:
                pdbid_descr_to_chains[pdbid_descr] = ''
            pdbid_descr_to_chains[pdbid_descr] += chain
        else:
            unparsed_titles.append(title)
            continue
    pdbid_descr_chains_list = [pdbid_descr + (chains,) for pdbid_descr, chains in pdbid_descr_to_chains.items()]
    return pdbid_descr_chains_list, unparsed_titles


def run_blastp(query_seq_fasta_file, evalue_threshold=1e-20, identity_threshold=0.9, database='pdb', num_threads=1):
    blastp_cline = NcbiblastpCommandline(query=query_seq_fasta_file, db=database, out='blast_result.xml',
                                         evalue=evalue_threshold, outfmt=5, num_threads=num_threads)
    try:
        blastp_cline()
    except ApplicationError as exc:
        error_msg = exc.stderr.strip()
        raise BLASTError(error_msg)
    with open('blast_result.xml') as f:
        blast_record = NCBIXML.read(f)
    all_unparsed_titles = []
    pdbid_to_blast_alignments = {}
    for alignment in blast_record.alignments:
        pdbid_descr_chains_list, unparsed_titles = parse_titles(alignment.title)
        all_unparsed_titles += unparsed_titles
        for hsp in alignment.hsps:
            hsp_dict = {'score': hsp.score,
                        'bits': hsp.bits,
                        'expect': hsp.expect,
                        'num_alignments': hsp.num_alignments,
                        'identities': hsp.identities,
                        'positives': hsp.positives,
                        'gaps': hsp.gaps,
                        'strand': hsp.strand,
                        'frame': hsp.frame,
                        'query': hsp.query,
                        'query_start': hsp.query_start,
                        'match': hsp.match,
                        'sbjct': hsp.sbjct,
                        'sbjct_start': hsp.sbjct_start,
                        'align_length': hsp.align_length,
                        # additional metrics
                        'identity': float(hsp.identities) / hsp.align_length,
                        'sbjct_length': alignment.length}

            if hsp_dict['identity'] >= identity_threshold:
                for pdbid, descr, chains in pdbid_descr_chains_list:
                    if pdbid not in pdbid_to_blast_alignments:
                        pdbid_to_blast_alignments[pdbid] = []
                    hsp_dict.update({'sbjct_pdbid': pdbid,
                                     'sbjct_chains': chains,
                                     'sbjct_description': descr})
                    pdbid_to_blast_alignments[pdbid].append(hsp_dict.copy())
    return pdbid_to_blast_alignments, all_unparsed_titles


@click.command()
@click.argument('query_seq_fasta_file', type=click.Path(exists=True))
@click.argument('output_json', type=click.Path())
@click.option('--evalue_threshold', type=float, default=1e-20, help='Expectation value threshold.')
@click.option('--identity_threshold', type=float, default=0.9, help='Sequence identity threshold.')
@click.option('--database', type=str, default='pdb', help='The database to use for BLAST protein search.')
@click.option('--num_threads', type=int, default=1, help='Number of threads to use for parallel processing.')
def _run_blastp_click(query_seq_fasta_file, output_json, evalue_threshold, identity_threshold, database, num_threads):
    '''
        BLAST-based script for database queries. Search for sequences similar to the one from
        QUERY_SEQ_FASTA_FILE in DATABASE. The result will be written in json format in OUTPUT_JSON.
    '''
    blast_alignments, unparsed_titles = run_blastp(query_seq_fasta_file, evalue_threshold, identity_threshold, database, num_threads)
    print('Unparsed titles: {}'.format(unparsed_titles))
    with open(output_json, 'w') as f:
        json.dump(blast_alignments, f, indent=4)


if __name__ == '__main__':
    _run_blastp_click()
