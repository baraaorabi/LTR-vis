#!/usr/bin/env python3
import argparse
from collections import namedtuple
from Bio.Seq import Seq

def parse_args():
    parser = argparse.ArgumentParser(
        description="Get specific gene and output its sequence")
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
    parser.add_argument("-go",
                        "--gene-out",
                        type=str,
                        required=True,
                        help="Path to output reads FASTA file")
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

def get_fasta_seq(fasta, contig, start, end, ori):
    contig_seq = list()
    flag = False
    for line in open(fasta):
        if line[0] == '>':
            if flag:
                break
            line = line.split()
            fasta_contig = line[0][1:]
            if contig in [fasta_contig, 'chr{}'.format(fasta_contig)] or fasta_contig in [contig, 'chr{}'.format(contig)]:
                flag = True
        elif flag:
            contig_seq.append(line.rstrip())
    contig_seq = ''.join(contig_seq)
    assert(start <= end <= len(contig_seq))
    seq = contig_seq[start:end]
    if ori == '+':
        return seq
    elif ori == '-':
        seq = str(Seq(str(seq)).reverse_complement())
        return seq
    else:
        raise Exception('Unexpected orientation: {}'.format(ori))

def output_gene(outpath, ginfo):
    outfile = open(outpath, 'w')
    fields = ['id', 'name', 'seq']
    print('>{}'.format(ginfo['id']), file=outfile)
    print('{}'.format(ginfo['seq']), file=outfile)
    outfile.close()

def main():
    args = parse_args()
    print('Getting {} gene info'.format(args.gene))
    ginfo = get_gene_info(gtf_path=args.gtf, gene_name=args.gene)
    ginfo['seq'] = get_fasta_seq(fasta=args.genome, contig=ginfo['chrom'], start=ginfo['start'], end=ginfo['end'], ori=ginfo['ori'])
    print('Outputting {} gene to {}'.format(args.gene, args.gene_out))
    output_gene(outpath=args.gene_out, ginfo=ginfo)

if __name__ == "__main__":
    main()
