import numpy as np
import argparse
try:
    from utils import get_sequence_lengths
except ImportError:
    from mumemto.utils import get_sequence_lengths


def remove_start_dollar(mums, s1_bv):
    new_mums = []
    l, starts, strands = mums
    dollars = np.where(s1_bv[starts[0] : starts[0]+l])[0]
    if len(dollars) == 0:
        new_mums.append((l, starts, strands))
    else:
        last_pos_left = 0
        last_pos_right = l
        for idx in range(len(dollars)):
            new_l = dollars[idx] - last_pos_left
            if new_l >= 20:
                new_starts = [s + last_pos_left if strand == '+' else s + last_pos_right - new_l for s, strand in zip(starts, strands)]
                new_mums.append((new_l, new_starts, strands))
            last_pos_left = dollars[idx] + 1
            last_pos_right = last_pos_right - new_l - 1
        new_l = l - last_pos_left
        if new_l >= 20:
            new_starts = [s + last_pos_left if strand == '+' else s for s, strand in zip(starts, strands)]
            new_mums.append((new_l, new_starts, strands))
    return new_mums
    
def main():
    parser = argparse.ArgumentParser(description='Merge MUMs files')
    parser.add_argument('--mums_file', '-m', required=True, help='Path to merged MUMs file')
    parser.add_argument('--s1_mums', '-s1', required=True, help='Path to s1 MUMs file')
    parser.add_argument('--s2_mums', '-s2', required=True, help='Path to s2 MUMs file')
    parser.add_argument('--output', '-o', help='Path to output merged MUMs file', default='merged.mums')
    
    args = parser.parse_args()
    
    if not args.mums_file.endswith('.mums'):
        args.mums_file += '.mums'
    if not args.s1_mums.endswith('.mums'):
        args.s1_mums += '.mums'
    if not args.s2_mums.endswith('.mums'):
        args.s2_mums += '.mums'

    with open(args.mums_file, 'r') as f:
        mums = [line.split() for line in f.read().splitlines()]
        mums = [[int(l[0]), [int(i) for i in l[1].split(',')], l[2].split(',')] for l in mums]

    with open(args.s1_mums, 'r') as f:
        s1 = [line.split() for line in f.read().splitlines()]
        s1 = [[int(l[0]), [int(i) for i in l[1].split(',')], l[2].split(',')] for l in s1]
    with open(args.s2_mums, 'r') as f:
        s2 = [line.split() for line in f.read().splitlines()]
        s2 = [[int(l[0]), [int(i) for i in l[1].split(',')], l[2].split(',')] for l in s2]
        
    ### get lengths
    s1_lens, s2_lens = get_sequence_lengths(args.mums_file.replace('.mums', '.lengths'), multilengths=True)

    ### build bitvectors for mum starts in s1 and s2
    s1_starts = np.cumsum(s1_lens)
    s2_starts = np.cumsum(s2_lens)
    s1_bv = np.zeros(sum(s1_lens) + 1, dtype=bool)
    s1_bv[s1_starts - 1] = 1
    s2_bv = np.zeros(sum(s2_lens) + 1, dtype=bool)
    s2_bv[s2_starts - 1] = 1
    breaks = [s1_starts, s2_starts]
    offsets = [np.concatenate(([0], s1_starts)), np.concatenate(([0], s2_starts))]

    ### get thresholds
    with open(args.s1_mums.replace('.mums', '.thresh'), 'rb') as f:
        s1_thresholds = np.fromfile(f, dtype=np.uint16)
    with open(args.s1_mums.replace('.mums', '.thresh_rev'), 'rb') as f:
        s1_thresholds_rev = np.fromfile(f, dtype=np.uint16)
    with open(args.s2_mums.replace('.mums', '.thresh'), 'rb') as f:
        s2_thresholds = np.fromfile(f, dtype=np.uint16)
    with open(args.s2_mums.replace('.mums', '.thresh_rev'), 'rb') as f:
        s2_thresholds_rev = np.fromfile(f, dtype=np.uint16)
        
    ### split grandMUMs that span multiple $ in the concatenated mums into match segments
    dollar_less = []
    for m in mums:
        dollar_less.extend(remove_start_dollar(m, s1_bv))

    ### find which mum a segment belongs to both subsets
    starts = np.array([m[1] for m in dollar_less])
    mum_idx = np.array([np.searchsorted(breaks[idx], starts[:,idx], side='right') for idx in range(2)]).transpose()

    ### main merging algorithm
    new_thresholds = []
    new_thresholds_rev = []
    merged = []
    for idx, (l, starts, strands) in enumerate(dollar_less):
        # first check if it is no longer unique
        mum_id1 = mum_idx[idx, 0]
        # get (left offset to start of match in mum, right offset to end of match in mum)
        offset1 = (starts[0] - offsets[0][mum_id1], offsets[0][mum_id1+1] - starts[0] - l - 1)
        # check that the match is still unique in s1
        thresh1 = s1_thresholds[starts[0]]
        if thresh1 == 0 or l <= thresh1:
            continue
        mum_id2 = mum_idx[idx, 1]
        # get (left offset to start of match in mum, right offset to end of match in mum)
        offset2 = (starts[1] - offsets[1][mum_id2], offsets[1][mum_id2 +1] - starts[1] - l - 1)
        # check that the match is still unique in s2
        thresh2 = s2_thresholds[starts[1]]
        if thresh2 == 0 or l <= thresh2:
            continue        
        
        ### do 1 first, then 2
        new_starts = []
        new_strands = []
        if strands[0] == '+': # technically always the case, but to be thorough
            for s, strand in zip(s1[mum_id1][1], s1[mum_id1][2]):
                new_starts.append(s + offset1[0] if strand == '+' else s + offset1[1])
                new_strands.append(strand)
        # else:
        #     for s, strand in zip(s1[mum_id1][1], s1[mum_id1][2]):
        #         new_starts.append(s + offset1 if strand == '-' else s + lendiff1)
        #         new_strands.append('-' if strand == '+' else '+')
        if strands[1] == '+':
            for s, strand in zip(s2[mum_id2][1], s2[mum_id2][2]):
                new_starts.append(s + offset2[0] if strand == '+' else s + offset2[1])
                new_strands.append(strand)
        else:
            for s, strand in zip(s2[mum_id2][1], s2[mum_id2][2]):
                new_starts.append(s + offset2[0] if strand == '+' else s + offset2[1])
                new_strands.append('-' if strand == '+' else '+')
        merged.append((int(l), tuple([int(x) for x in new_starts]), tuple(new_strands)))
        
        curs1_thresh = (s1_thresholds[starts[0] : starts[0] + l], s1_thresholds_rev[offsets[0][mum_id1] + offset1[1]: offsets[0][mum_id1+1] - 1 - offset1[0]])
        curs2_thresh = (s2_thresholds[starts[1] : starts[1] + l], s2_thresholds_rev[offsets[1][mum_id2] + offset2[1]: offsets[1][mum_id2+1] - 1 - offset2[0]])
        if strands[1] == '-':    
            new_thresholds.extend(np.where((curs1_thresh[0] > 0) & (curs2_thresh[1] > 0), np.maximum(curs1_thresh[0], curs2_thresh[1]), np.zeros_like(curs1_thresh[0])))
            new_thresholds.extend([0])
            new_thresholds_rev.extend(np.where((curs1_thresh[1] > 0) & (curs2_thresh[0] > 0), np.maximum(curs1_thresh[1], curs2_thresh[0]), np.zeros_like(curs1_thresh[1])))
            new_thresholds_rev.extend([0])
        else:
            new_thresholds.extend(np.where((curs1_thresh[0] > 0) & (curs2_thresh[0] > 0), np.maximum(curs1_thresh[0], curs2_thresh[0]), np.zeros_like(curs1_thresh[0])))
            new_thresholds.extend([0])
            new_thresholds_rev.extend(np.where((curs1_thresh[1] > 0) & (curs2_thresh[1] > 0), np.maximum(curs1_thresh[1], curs2_thresh[1]), np.zeros_like(curs1_thresh[1])))
            new_thresholds_rev.extend([0])

    
    ### write output
    # if args.output:
    if not args.output.endswith('.mums'):
        args.output += '.mums'
    with open(args.output, 'w') as f:
        for m in merged:
            f.write('%d\t%s\t%s\n' % (m[0], ','.join(map(str, m[1])), ','.join(m[2])))
    with open(args.output.replace('.mums', '.thresh'), 'wb') as f:
        f.write(np.array(new_thresholds, dtype=np.uint16).tobytes())
    with open(args.output.replace('.mums', '.thresh_rev'), 'wb') as f:
        f.write(np.array(new_thresholds_rev, dtype=np.uint16).tobytes())
    # else:
    #     try:
    #         for m in merged:
    #             sys.stdout.write('%d\t%s\t%s\n' % (m[0], ','.join(map(str, m[1])), ','.join(m[2])))
    #     except BrokenPipeError:
    #         sys.stdout = os.fdopen(sys.stdout.fileno(), 'wb', 0)
    #         sys.exit(0)
    with open(args.output.replace('.mums', '.lengths'), 'w') as out:
        with open(args.s1_mums.replace('.mums', '.lengths'), 'r') as f:
            out.write(f.read().strip() + '\n')
        with open(args.s2_mums.replace('.mums', '.lengths'), 'r') as f:
            out.write(f.read().strip() + '\n')
    
if __name__ == "__main__":
    main()
