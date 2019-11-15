#!/usr/bin/env python3
import argparse
import json

def parse_args():
    parser = argparse.ArgumentParser(
        description="Turn FASTQ into JSON format")
    parser.add_argument("-q",
                        "--fastq",
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

def parse_fastq(fastq):
    rid_to_data = dict()
    for idx,line in enumerate(open(fastq)):
        if idx % 4 == 0:
            rid = line[1:].rstrip().split()[0]
            seq = list()
            qual = list()
            assert(not rid in rid_to_data)
        if idx % 4 == 1:
            seq = [dna_to_int[c] for c in line.rstrip()]
        if idx % 4 == 3:
            qual = [ord(c)-33 for c in line.rstrip()]
            assert(len(qual)==len(seq))
            rid_to_data[rid] = dict(
                seq=seq,
                qual=qual
            )
    assert(idx % 4 == 3)
    return rid_to_data

def main():
    args = parse_args()
    rid_to_data = parse_fastq(fastq=args.fastq)

    data = list()
    for rid,rdata in rid_to_data.items():
        data.append(dict(
            name=rid,
            seq=rdata['seq'],
            qual=rdata['qual'],
        ))
    outfile = open(args.json_out, 'w')
    print(json.dumps(data, indent=4), file=outfile)
    outfile.close()
    # rid_to_tblocks = get_target_blocks(tlength=tlength, rid_to_aln=rid_to_aln, rid_to_len=rid_to_len)

if __name__ == "__main__":
    main()
