#!/usr/bin/env python3
import argparse
import re
from collections import namedtuple
import json
# Aln = namedtuple('Aln', [
#     'qstart',
#     'qend',
#     'tstart',
#     'tend',
#     'ori',
# ])

def parse_args():
    parser = argparse.ArgumentParser(
        description="Turn reads PAF alignments into JSON format")
    parser.add_argument("-p",
                        "--paf",
                        type=str,
                        required=True,
                        help="Path to PAF file")
    parser.add_argument("-t",
                        "--target",
                        type=str,
                        default='ENSG00',
                        help="Name of target")
    parser.add_argument("-jo",
                        "--json-out",
                        type=str,
                        required=True,
                        help="Path to output JSON file")
    args = parser.parse_args()
    return args

def get_aln_info(paf, target):
    rid_to_blks = dict()
    tlen = -1
    for line in open(paf):
        line = line.rstrip().split('\t')
        if not target in line[5]:
            continue
        rid = line[0]
        assert(not rid in rid_to_blks)
        if tlen == -1:
            tlen = int(line[6])
        else:
            assert(tlen == int(line[6]))
        qlen = int(line[1])
        qstart = int(line[2])
        qend = int(line[3])
        tstart = int(line[7])
        tend = int(line[8])
        ori = line[4]
        assert(ori in ['-','+'])
        cigar = ''
        for f in line[12:]:
            if f[:len('cg:Z:')] == 'cg:Z:':
                cigar = f[len('cg:Z:'):]
                cigar = re.findall(r'(\d+)([A-Z]{1})', cigar)
        rid_to_blks[rid] = list()
        rid_to_blks[rid].append(dict(
            type='d',
            start=0,
            end=tstart,
        ))
        rid_to_blks[rid].append(dict(
            type='i',
            pos=tstart,
            length=qstart,
        ))
        tpos1 = tstart
        qpos1 = qstart
        tpos2 = tstart
        qpos2 = qstart
        for size,op in cigar:
            size = int(size)
            if op == 'M':
                tpos2 += size
                qpos2 += size
            elif op == 'I':
                tpos2 += 0
                qpos2 += size
                if size > 15:
                    rid_to_blks[rid].append(dict(
                        type='m',
                        start=tpos1,
                        end=tpos2,
                        orientation=ori,
                    ))
                    rid_to_blks[rid].append(dict(
                        type='i',
                        pos=tpos2,
                        length=size,
                    ))
                    tpos1 = tpos2
                    qpos1 = qpos2
            elif op == 'D':
                tpos2 += size
                qpos2 += 0
            elif op == 'N':
                rid_to_blks[rid].append(dict(
                    type='m',
                    start=tpos1,
                    end=tpos2,
                    orientation=ori,
                ))
                tpos2 += size
                qpos2 += 0
                rid_to_blks[rid].append(dict(
                    type='d',
                    start=tpos1,
                    end=tpos2,
                ))
                tpos1 = tpos2
                qpos1 = qpos2
        if tpos2 > tpos1:
            rid_to_blks[rid].append(dict(
                type='m',
                start=tpos1,
                end=tpos2,
                orientation=ori,
            ))
        rid_to_blks[rid].append(dict(
            type='i',
            pos=tend,
            length=qlen-qend,
        ))
        rid_to_blks[rid].append(dict(
            type='d',
            start=tend,
            end=tlen,
        ))
    return rid_to_blks,tlen

def main():
    args = parse_args()
    rid_to_blks,tlen = get_aln_info(paf=args.paf, target=args.target)

    data = list()
    for rid,alns in rid_to_blks.items():
        data.append(dict(
            target_length=tlen,
            name=rid,
            data=alns
        ))
    outfile = open(args.json_out, 'w')
    print(json.dumps(data, indent=4), file=outfile)
    outfile.close()
    # rid_to_tblocks = get_target_blocks(tlength=tlength, rid_to_aln=rid_to_aln, rid_to_len=rid_to_len)

if __name__ == "__main__":
    main()
