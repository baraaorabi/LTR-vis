#!/usr/bin/env python3
import argparse
from pyfasta import Fasta
from Bio.Seq import Seq

def parse_args():
    parser = argparse.ArgumentParser(
        description="Get reads that map to a specific gene")
    # parser.add_argument("-p",
    #                     "--paf",
    #                     type=str,
    #                     required=True,
    #                     help="Path to PAF file")
    # parser.add_argument("-r",
    #                     "--reads",
    #                     type=str,
    #                     required=True,
    #                     help="Path to reads FASTQ file")
    parser.add_argument("-g",
                        "--gene",
                        type=str,
                        required=True,
                        help="Gene name")
    parser.add_argument("-ar",
                        "--gtf",
                        type=str,
                        required=True,
                        help="Path to reference annotation GTF file")
    parser.add_argument("-gr",
                        "--genome",
                        type=str,
                        required=True,
                        help="Path to reference genome FASTA file")
    # parser.add_argument("-ro",
    #                     "--reads-out",
    #                     type=str,
    #                     required=True,
    #                     help="Path to output reads FASTQ file")
    # parser.add_argument("-go",
    #                     "--gene-out",
    #                     type=str,
    #                     required=True,
    #                     help="Path to output reads FASTA file")
    args = parser.parse_args()
    return args

def get_gene_info(gtf_path, gene_name):
    ginfo = dict(name='', id='', transcript_ids=list(), start=-1, end=-1, chrom='', ori='', seq='')
    for line in open(gtf_path):
        if line[0] == '#':
            continue
        line = line.rstrip().split('\t')
        if line[2] != 'gene':
            continue
        attributes = dict()
        for att in line[8].split(';'):
            if len(att) < 2:
                continue
            att = att.split('"')
            k = att[0].strip('" ')
            v = att[1].strip('" ')
            attributes[k] = v
        if attributes['gene_name'] != gene_name and attributes['gene_id'] != gene_name:
            continue
        if ginfo['name'] != '':
            raise Exception('Ambigious gene name/id: {} vs {}'.format(attributes['gene_name'], ginfo['name']))
        ginfo['name'] = attributes['gene_name']
        ginfo['id'] = attributes['gene_id']
        ginfo['chrom'] = line[0]
        ginfo['start'] = int(line[3])-1
        ginfo['end'] = int(line[4])-1
        ginfo['ori'] = line[6]
    for line in open(gtf_path):
        if line[0] == '#':
            continue
        line = line.rstrip().split('\t')
        if line[2] != 'transcript':
            continue
        attributes = dict()
        for att in line[8].split(';'):
            if len(att) < 2:
                continue
            att = att.split('"')
            k = att[0].strip('" ')
            v = att[1].strip('" ')
            attributes[k] = v
        if attributes['gene_id'] != ginfo['id']:
            continue
        ginfo['transcript_ids'].append(attributes['transcript_id'])
    return ginfo

def get_seq(fasta, contig, start, end, ori):
    ref = Fasta(fasta)
    if contig in ref:
        seq = ref[contig][start:end]
    elif 'chr{}'.format(contig) in ref:
        seq = ref['chr{}'.format(contig)][start:end]
    elif len(contig)>3 and contig[:3] == 'chr':
        seq = ref[contig[3:]][start:end]
    else:
        raise Exception('Cannot find {} contig in ref ({})'.format(contig, sorted(ref.keys())) )
    if ori == '+':
        return
    elif ori == '-':
        seq = str(Seq(str(seq)).reverse_complement())
    else:
        raise Exception('Unexpected orientation: {}'.format(ori))
    return seq

def main():
    args = parse_args()
    ginfo = get_gene_info(gtf_path=args.gtf, gene_name=args.gene)
    ginfo['seq'] = get_seq(fasta=args.genome, contig=ginfo['chrom'], start=ginfo['start'], end=ginfo['end'], ori=ginfo['ori'])
    print(ginfo)

if __name__ == "__main__":
    main()
