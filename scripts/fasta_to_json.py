#!/usr/bin/env python3
import argparse
import json

def parse_args():
    parser = argparse.ArgumentParser(
        description="Turn FASTA into JSON format")
    parser.add_argument("-f",
                        "--fasta",
                        type=str,
                        required=True,
                        help="Path to FASTA file")
    parser.add_argument("-jo",
                        "--json-out",
                        type=str,
                        required=True,
                        help="Path to output JSON file")
    args = parser.parse_args()
    return args

dna_to_int = dict(
    A=0,
    C=1,
    G=2,
    T=3,
    a=0,
    c=1,
    g=2,
    t=3,
    N=3,
    n=3,
)

def parse_fasta(fasta):
    lines = open(fasta).readlines()
    assert(len(lines)==2)

    target = lines[0].rstrip().split()[0][1:]
    seq = [dna_to_int[c] for c in lines[1].rstrip()]
    qual = [60 for _ in range(len(lines[1].rstrip()))]
    return target,seq,qual

def main():
    args = parse_args()
    target,seq,qual = parse_fasta(fasta=args.fasta)
    data = dict(
        name=target,
        seq=seq,
        qual=qual
    )
    outfile = open(args.json_out, 'w')
    print(json.dumps(data, indent=4), file=outfile)
    outfile.close()
    # rid_to_tblocks = get_target_blocks(tlength=tlength, rid_to_aln=rid_to_aln, rid_to_len=rid_to_len)

if __name__ == "__main__":
    main()
