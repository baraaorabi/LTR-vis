#!/usr/bin/env python3
import argparse
from collections import namedtuple
from Bio.Seq import Seq

ReadAln = namedtuple('ReadAln', [
    'query',
    'qstart',
    'qend',
    'target',
    'tstart',
    'tend',
    'ori',
])

def parse_args():
    parser = argparse.ArgumentParser(
        description="Get reads that map to a specific gene")
    parser.add_argument("-p",
                        "--paf",
                        type=str,
                        required=True,
                        help="Path to PAF file")
    parser.add_argument("-r",
                        "--reads",
                        type=str,
                        required=True,
                        help="Path to reads FASTQ file")
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
    parser.add_argument("-ro",
                        "--reads-out",
                        type=str,
                        required=True,
                        help="Path to output reads TSV file")
    # parser.add_argument("-ao",
    #                     "--alns-out",
    #                     type=str,
    #                     required=True,
    #                     help="Path to output read alignments to whole transcriptome TSV file")
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

def get_read_ids_from_paf(paf, ginfo):
    read_ids = set()
    read_oris = dict()
    targets = set()
    targets.add(ginfo['id'])
    for tid in ginfo['transcript_ids']:
        targets.add(tid)
    for line in open(paf):
        line = line.rstrip().split('\t')
        if not line[5].split('.')[0] in targets:
            continue
        if not 'tp:A:P' in line[12:]:
            continue
        rid = line[0]
        read_ids.add(rid)
        read_oris[rid] = line[4]
    # for line in open(paf):
    #     line = line.rstrip().split('\t')
    #     rid = line[0]
    #     if not rid in read_ids:
    #         continue
    #     read_alns.append(
    #         ReadAln(
    #             query  = rid,
    #             qstart = int(line[2]),
    #             qend   = int(line[3]),
    #             ori    = line[4],
    #             target = line[5].split('.')[0],
    #             tstart = int(line[7]),
    #             tend   = int(line[8]),
    #         )
    #     )
    return read_ids,read_oris

def output_reads(outpath, fastq, read_ids, read_oris, gene_ori):
    rid_to_seq = dict()
    flag = False
    for idx,line in enumerate(open(fastq)):
        if idx % 4 == 0:
            rid = line[1:].rstrip().split()[0]
            flag = rid in read_ids
        if not flag:
            continue
        if idx % 4 == 1:
            rid_to_seq[rid] = [line.rstrip()]
        if idx % 4 == 3:
            rid_to_seq[rid].append(line.rstrip())
    assert(idx % 4 == 3)

    outfile = open(outpath,'w')
    for rid,(seq,qual) in rid_to_seq.items():
        if read_oris[rid] == gene_ori:
            seq = str(Seq(str(seq)).reverse_complement())
            print(rid,read_oris[rid], gene_ori)
        print('\n'.join(['@{}'.format(rid), seq, '+', qual]), file=outfile)
    outfile.close()

# def output_alns_tsv(outpath, read_alns):
#     outfile = open(outpath, 'w')
#     print('\t'.join([str(f) for f in ReadAln._fields]), file=outfile)
#     for aln in read_alns:
#         print('\t'.join([str(f) for f in aln]), file=outfile)
#     outfile.close()

def main():
    args = parse_args()
    print('Getting {} gene info'.format(args.gene))
    ginfo = get_gene_info(gtf_path=args.gtf, gene_name=args.gene)
    print('Getting read alignments from {}'.format(args.paf))
    read_ids,read_oris = get_read_ids_from_paf(paf=args.paf, ginfo=ginfo)

    # print('Outputting read transcriptome alignments to {}'.format(args.alns_out))
    # output_alns_tsv(outpath=args.alns_out, read_alns=read_alns)
    print('Outputting reads ({}) to {}'.format(args.reads, args.reads_out))
    output_reads(outpath=args.reads_out, fastq=args.reads, read_ids=read_ids, read_oris=read_oris, gene_ori=ginfo['ori'])

if __name__ == "__main__":
    main()
