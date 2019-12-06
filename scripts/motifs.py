#!/usr/bin/env python3
import argparse
import re

def parse_args():
    parser = argparse.ArgumentParser(
        description="Get reads that map to a specific gene")
    parser.add_argument("-o",
                        "--output",
                        type=str,
                        required=True,
                        help="Path to out TSV file")
    # parser.add_argument("-r",
    #                     "--reads",
    #                     type=str,
    #                     required=True,
    #                     help="Path to reads FASTQ file")
    parser.add_argument("-g",
                        "--gene",
                        type=str,
                        required=True,
                        help="Gene FASTA")
    parser.add_argument("-m",
                        "--motifs",
                        type=str,
                        required=True,
                        help="Gene name")
    args = parser.parse_args()
    return args

# def get_reads(fastq):
#     reads = dict()
#     for idx,l in enumerate(open(fastq)):
#         l = l.rstrip()
#         if idx%4 == 0:
#             q = l[1:]
#         elif idx%4==1:
#             s = l
#         elif idx%4==3:
#             reads[q]=s
#     return reads

def get_gene(fastq):
    for idx,l in enumerate(open(fastq)):
        l = l.rstrip()
        if idx%4 == 0:
            q = l[1:]
        elif idx%4==1:
            s = l
    return q,s

def get_motifs(tsv):
    motifs = dict()
    for l in open(tsv):
        l = l.rstrip().split('\t')
        motifs[l[0]]=l[1]
    return motifs

def get_matches(motifs, reads, outpath):
    outfile = open(outpath, 'w+')
    for mk,mv in motifs.items():
        print(mv)
        for k,v in reads.items():
            for i in re.finditer(mv,v):
                print(k,mk,i.start(),i.end(),sep='\t',file=outfile)
    outfile.close()


def main():
    args = parse_args()
    reads = dict()
    q,s = get_gene(args.gene)
    reads[q]=s
    motifs = get_motifs(args.motifs)
    get_matches(motifs,reads,args.output)

if __name__ == "__main__":
    main()
