#!/usr/bin/env python3
import argparse
import re
from collections import namedtuple
import json
from itertools import chain
from itertools import islice

dna_to_int = dict(
    A=1,
    C=2,
    G=3,
    T=4,
    N=4,
)

def str2bool(v):
    if isinstance(v, bool):
       return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')

VERBOSE = False

def parse_args():
    parser = argparse.ArgumentParser(
        description="Turn reads PAF alignments into JSON format")
    parser.add_argument("-q",
                        "--query",
                        type=str,
                        required=True,
                        help="Path to FASTA/Q file for queries")
    parser.add_argument("-t",
                        "--target",
                        type=str,
                        required=True,
                        help="Path to FASTA/Q file for target")
    parser.add_argument("-p",
                        "--paf",
                        type=str,
                        required=True,
                        help="Path to PAF file")
    parser.add_argument("-jo",
                        "--json-out",
                        type=str,
                        required=True,
                        help="Path to output JSON file")
    parser.add_argument("--verbose", type=str2bool, nargs='?',
                            const=True, default=False,
                            help="Verbose")
    args = parser.parse_args()
    return args

def parse_paf(paf):
    queries = dict()
    paf = open(paf)

    line = paf.readline()
    trgt = line.rstrip().split('\t')[5]
    trgt_l = int(line.rstrip().split('\t')[6])

    for line in chain([line], paf):
        line = line.rstrip().split('\t')
        assert(trgt   == line[5])
        assert(trgt_l == int(line[6]))
        quer = line[0]
        assert(not quer in queries)
        quer_l = int(line[1])
        quer_s = int(line[2])
        quer_e = int(line[3])
        quer_o = line[4]
        assert(quer_o in ['-','+'])
        trgt_s = int(line[7])
        trgt_e = int(line[8])
        cigar = list()
        for f in line[12:]:
            if f[:len('cg:Z:')] == 'cg:Z:':
                cigar = f[len('cg:Z:'):]
                cigar = re.findall(r'(\d+)([M|I|D|N|S|H|P|=|X]{1})', cigar)
                cigar = [(int(s),o) for s,o in cigar]
        queries[quer] = dict(
            quer_l=quer_l,
            quer_s=quer_s,
            quer_e=quer_e,
            quer_o=quer_o,
            trgt_s=trgt_s,
            trgt_e=trgt_e,
            cigar=cigar,
        )
    return trgt,trgt_l,queries

def get_target(target, target_name):
    target = open(target, 'r')
    line = target.readline()
    if line[0] == '@':
        is_fa = False
    elif line[0] == '>':
        is_fa = True
    else:
        raise Exception('File {} must be FASTA/Q formatted'.format(target))

    in_target = False
    for idx,line in enumerate(chain([line], target)):
        if is_fa:
            if idx % 2 == 0:
                in_target = target_name in line[1:].rstrip().split('.')[0]
            elif in_target:
                seq = [dna_to_int[c] for c in line.rstrip().upper()]
                break
        else:
            if idx % 4 == 0:
                in_target = target_name == line[1:].rstrip().split('.')[0]
            elif idx % 4 == 1 and in_target:
                seq = [dna_to_int[c] for c in line.rstrip().upper()]
            elif idx % 4 == 3 and in_target:
                qual = [ord(c)-33 for c in line.rstrip()]
                break
    if is_fa:
        qual = [ord(c)-33 for c in '-'*len(seq)]
    assert(len(qual)==len(seq))
    return seq,qual

def get_queries(queries, query_file):
    query_file = open(query_file, 'r')
    line = query_file.readline()
    if line[0] == '@':
        is_fa = False
    elif line[0] == '>':
        is_fa = True
    else:
        raise Exception('File {} must be FASTA/Q formatted'.format(target))
    for idx,line in enumerate(chain([line], query_file)):
        if is_fa:
            if idx % 2 == 0:
                quer = line[1:].rstrip().split('.')[0]
            else:
                seq = [dna_to_int[c] for c in line.rstrip().upper()]
                qual = [ord(c)-33 for c in '-'*len(seq)]
                if not quer in queries:
                    queries[quer] = dict(
                        quer_l=len(seq),
                        quer_s=0,
                        quer_e=0,
                        quer_o='+',
                        trgt_s=0,
                        trgt_e=0,
                        cigar=list(),
                    )
                queries[quer]['seq'] = seq
                queries[quer]['qual'] = qual
                assert(len(seq) == queries[quer]['quer_l'])
        else:
            if idx % 4 == 0:
                quer = line[1:].rstrip().split('.')[0]
            elif idx % 4 == 1:
                seq = [dna_to_int[c] for c in line.rstrip().upper()]
            elif idx % 4 == 3:
                qual = [ord(c)-33 for c in line.rstrip()]
                assert(len(qual)==len(seq))
                if not quer in queries:
                    queries[quer] = dict(
                        quer_l=len(seq),
                        quer_s=0,
                        quer_e=0,
                        quer_o='+',
                        trgt_s=0,
                        trgt_e=0,
                        cigar=list(),
                    )
                queries[quer]['seq'] = seq
                queries[quer]['qual'] = qual
                assert(len(seq) == queries[quer]['quer_l'])

def compute_cigar(queries, trgt_l):
    for k,q in queries.items():
        read_length = 0
        t_to_q  = [-1 for _ in range(trgt_l)]
        intervals = list()
        t_pos = 0
        q_pos = q['quer_s']
        read_length += q_pos
        for pos in range(0, q['trgt_s']):
            t_pos+=1
        if t_pos > 0:
            intervals.append(dict(orientation=q['quer_o'], type='d', size=t_pos, start=0, end=t_pos-1))
        if q_pos > 0:
            intervals.append(dict(orientation=q['quer_o'], type='i', size=q_pos, start=q['trgt_s'], end=q['trgt_s']))
        for size,op in q['cigar']:
            if op == 'M' or op == 'X':
                if len(intervals) == 0 or intervals[-1]['type'] != 'm':
                    intervals.append(dict(orientation=q['quer_o'], type='m',size=0,start=t_pos,end=t_pos-1))
                read_length+= size
                for pos in range(t_pos, t_pos+size):
                    t_to_q[pos] = q_pos
                    q_pos += 1
                    t_pos += 1
                    intervals[-1]['end'] += 1
                    intervals[-1]['size'] += 1
            elif op == 'I':
                read_length+= size
                if size > 5:
                    intervals.append(dict(orientation=q['quer_o'], type='i',size=size,start=t_pos,end=t_pos))
                else:
                    pass
                q_pos += size
            elif op == 'D' or op == 'N':
                if size > 20:
                    intervals.append(dict(orientation=q['quer_o'], type='d',size=size,start=t_pos,end=t_pos+size-1))
                else:
                    if intervals[-1]['type'] != 'm':
                        intervals.append(dict(orientation=q['quer_o'], type='m',size=0,start=t_pos,end=t_pos-1))
                    intervals[-1]['end'] += size
                    intervals[-1]['size'] += size
                t_pos += size
            else:
                raise Exception('Dunno what to do with cigar {}:{}'.format(op,size))
        assert(q['quer_e']==read_length), [q['quer_e'], read_length]
        for i in intervals:
            assert(i['type']=='i' or i['size'] == i['end']-i['start']+1),i
        assert (q_pos ==  q['quer_e'])
        assert (t_pos ==  q['trgt_e']),[t_pos, q['trgt_e']]
        if q['quer_e'] > 0 and q['quer_l'] - q['quer_e'] > 0:
            intervals.append(dict(orientation=q['quer_o'], type='i', size=q['quer_l'] - q['quer_e'], start=q['trgt_e'], end=q['trgt_e']))
        size = (trgt_l) - (q['trgt_e'])
        if size > 0:
            intervals.append(dict(orientation=q['quer_o'], type='d', size=size, start=q['trgt_e'], end=trgt_l-1))

        queries[k]['t_to_q']    = t_to_q
        queries[k]['intervals'] = intervals

def output_json(trgt, trgt_l, trgt_seq, trgt_qual, queries, outpath):
    outfile = open(outpath, 'w')
    data = dict(
        target=dict(
            name=trgt,
            length=trgt_l,
            seq=trgt_seq,
            qual=trgt_qual
        ),
        queries=[
            dict(
                name = k,
                seq  = [q['seq'][q['t_to_q'][p]]  if q['t_to_q'][p] != -1 else 0 for p in range(trgt_l)],
                qual = [q['qual'][q['t_to_q'][p]] if q['t_to_q'][p] != -1 else 0 for p in range(trgt_l)],
                intervals = q['intervals']
            )
            for k,q in queries.items()
        ]
    )
    print(json.dumps(data, indent=4), file=outfile)
    outfile.close()
    if not VERBOSE:
        return
    print(''.join(str(x) for x in data['target']['seq']),data['target']['name'])
    print(''.join(chr(x+33) for x in data['target']['qual']))
    for d in data['queries']:
        # if d['name'] != 'READ2':
        #     continue
        q = queries[d['name']]
        print(''.join(str(x%10) if x != -1 else '!' for x in queries[d['name']]['t_to_q'] ))
        print(''.join(str(x) for x in d['seq']),d['name'],'\n',  queries[d['name']]['cigar'], d['intervals'])
        print(''.join(chr(x+33) for x in d['qual']))


def main():
    args = parse_args()
    global VERBOSE
    VERBOSE = args.verbose
    trgt,trgt_l,queries = parse_paf(paf=args.paf)
    if VERBOSE:
        print(trgt,trgt_l)
    trgt_seq,trgt_qual = get_target(target=args.target, target_name=trgt)
    get_queries(queries=queries, query_file=args.query)
    compute_cigar(queries=queries, trgt_l=trgt_l)
    output_json(trgt=trgt, trgt_l=trgt_l, trgt_seq=trgt_seq, trgt_qual=trgt_qual, queries=queries, outpath=args.json_out)

    # for q in queries:
    #     print('quer_l', queries[q]['quer_l'])
    #     print('quer_s', queries[q]['quer_s'])
    #     print('quer_e', queries[q]['quer_e'])
    #     print('quer_o', queries[q]['quer_o'])
    #     print('trgt_s', queries[q]['trgt_s'])
    #     print('trgt_e', queries[q]['trgt_e'])
    #     for c in queries[q]['cigar']:
    #         print(c)
    #     # print('t_to_q', queries[q]['t_to_q'])
    #     for i in queries[q]['intervals']:
    #         print(i)
            # print(''.join([str(queries[q]['seq'][queries[q]['t_to_q'][p]]) if queries[q]['t_to_q'][p] != -1 else '-' for p in  range(i['start'],i['end']) ]))



if __name__ == "__main__":
    main()
