from Bio import SeqIO
from Bio.Seq import Seq
import argparse
import numpy as np
from tqdm.auto import tqdm
import os

def parse_arguments():
    parser = argparse.ArgumentParser(description='Extract the MUM sequences')
    parser.add_argument('-m', '--mumfile', type=str, help='bumfile to process')
    parser.add_argument('-i', '--index', type=int, default=0, help='The index of the file in the corresponding filelist')
    parser.add_argument('-o', '--output', type=str, help='The name of the output file', default="mums.fa")
    parser.add_argument('-f', '--filelist', type=str, help='filelist file')
    
    args = parser.parse_args()
    if args.filelist == None:
        args.filelist = os.path.splitext(args.mumfile)[0] + '.lengths'
    return args

def from_bums(bumsfile):
    with open(bumsfile, 'rb') as f:
        n_seqs, n_mums = np.fromfile(f, count = 2, dtype=np.uint64)
        mum_lengths = np.fromfile(f, count = n_mums, dtype=np.int64)
        mum_starts = np.fromfile(f, count = n_seqs * n_mums, dtype=np.int64).reshape(n_mums, n_seqs)
        mum_strands = np.fromfile(f, dtype=np.uint8)
        mum_strands = np.unpackbits(mum_strands, count=n_mums * n_seqs).reshape(n_mums, n_seqs).astype(bool)
    return mum_lengths, mum_starts, mum_strands

def main():
    args = parse_arguments()
    if args.filelist == None:
        args.filelist = args.mumfile.replace('.bums', '.lengths')
    line = open(args.filelist, 'r').read().splitlines()[args.index]
    file = line.split()[0]
    print(file)
    seq = Seq(''.join([l for l in open(file, 'r').read().splitlines() if not l.startswith('>')]))
    # seq = Seq('').join([record.seq for record in SeqIO.parse(file, "fasta")])
    assert len(seq) == int(line.split()[-1])
    mum_lengths, mum_starts, mum_strands = from_bums(args.mumfile)
    with open(args.output, 'w') as out:
        outlines = []
        for i in tqdm(range(len(mum_lengths))):
            cur = seq[mum_starts[i, args.index] : mum_starts[i, args.index] + mum_lengths[i]]  
            if mum_strands[i, args.index]:
                outlines.append(str(cur))
            else:
                outlines.append(str(cur.reverse_complement()))
        out.write('\n'.join(outlines))
if __name__ == "__main__":
    main()

