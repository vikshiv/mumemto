import numpy as np
import argparse
try:
    from utils import get_sequence_lengths, parse_mums_generator
except ImportError:
    from mumemto.utils import get_sequence_lengths, parse_mums_generator


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
                new_starts = [s + last_pos_left if strand else s + last_pos_right - new_l for s, strand in zip(starts, strands)]
                new_mums.append((new_l, new_starts, strands))
            last_pos_left = dollars[idx] + 1
            last_pos_right = last_pos_right - new_l - 1
        new_l = l - last_pos_left
        if new_l >= 20:
            new_starts = [s + last_pos_left if strand else s for s, strand in zip(starts, strands)]
            new_mums.append((new_l, new_starts, strands))
    return new_mums

def main():
    parser = argparse.ArgumentParser(description='Merge MUMs files')
    parser.add_argument('--merged_mums', '-m', required=True, help='Path to MUMs of MUMs file')
    parser.add_argument('mum_files', metavar='MUM_FILES', nargs='+', help='Paths to MUMs files to merge')
    parser.add_argument('--output', '-o', help='Path to output merged MUMs file', default='merged.mums')
    
    args = parser.parse_args()
    
    assert len(args.mum_files) >= 2, "At least two MUMs files are required for merging"
    
    if not args.merged_mums.endswith('.mums'):
        args.merged_mums += '.mums'
    for i in range(len(args.mum_files)):
        if not args.mum_files[i].endswith('.mums'):
            args.mum_files[i] += '.mums'

    premerge_mums = [list(parse_mums_generator(m)) for m in args.mum_files]
    
    ### get lengths
    mum_lens = get_sequence_lengths(args.merged_mums.replace('.mums', '.lengths'), multilengths=True)

    ### build bitvectors for mum starts in comb mums
    mum_starts = [np.cumsum(lens) for lens in mum_lens]
    set1_bv = np.zeros(sum(mum_lens[0]) + 1, dtype=bool)
    set1_bv[mum_starts[0] - 1] = 1
    mum_offsets = [np.concatenate(([0], starts)) for starts in mum_starts]

    ### get thresholds
    thresholds, rev_thresholds = [], []
    for m in args.mum_files:
        with open(m.replace('.mums', '.thresh'), 'rb') as f:
            thresholds.append(np.fromfile(f, dtype=np.uint16))
        with open(m.replace('.mums', '.thresh_rev'), 'rb') as f:
            rev_thresholds.append(np.fromfile(f, dtype=np.uint16))
    
    NUM_SETS = len(mum_lens)
    
    ### split grandMUMs that span multiple $ in the concatenated mums into match segments
    dollar_less = []
    for m in parse_mums_generator(args.merged_mums):
        dollar_less.extend(remove_start_dollar(m, set1_bv))

    ### find which mum a segment belongs to both subsets
    starts = np.array([m[1] for m in dollar_less])
    mum_idx = np.array([np.searchsorted(mum_starts[idx], starts[:,idx], side='right') for idx in range(NUM_SETS)]).transpose()

    ### main merging algorithm
    new_thresholds = []
    new_thresholds_rev = []
    merged = []
    for idx, (l, starts, strands) in enumerate(dollar_less):
        # first check if it is no longer unique
        offset = []
        thresh_check = True
        for i in range(NUM_SETS):
            mumid = mum_idx[idx, i]
            # get (left offset to start of match in mum, right offset to end of match in mum)
            offset.append((starts[i] - mum_offsets[i][mumid], mum_offsets[i][mumid+1] - starts[i] - l - 1))
            # check that the match is still unique in subset
            thresh = thresholds[i][starts[i]]
            if thresh == 0 or l <= thresh:
                thresh_check = False
                break
        if not thresh_check:
            continue
        
        new_starts = []
        new_strands = []
        for i in range(NUM_SETS):
            mumid = mum_idx[idx, i]
            for s, strand in zip(premerge_mums[i][mumid][1], premerge_mums[i][mumid][2]): # get matching mum
                new_starts.append(s + offset[i][0] if strand else s + offset[i][1])
                # the mum in set i matches in the forward direction
                new_strands.append(strand if strands[i] else not strand)
        merged.append((int(l), tuple([int(x) for x in new_starts]), tuple(new_strands)))
    
        cur_thresh = []
        cur_revthresh = []
        for i in range(NUM_SETS):
            thresh = (thresholds[i][starts[i] : starts[i] + l], rev_thresholds[i][mum_offsets[i][mum_idx[idx, i]] + offset[i][1]: mum_offsets[i][mum_idx[idx, i]+1] - 1 - offset[i][0]])
            if strands[i]:
                cur_thresh.append(thresh[0]); cur_revthresh.append(thresh[1])
            else:
                cur_revthresh.append(thresh[0]); cur_thresh.append(thresh[1])
        cur_thresh = np.array(cur_thresh)
        cur_revthresh = np.array(cur_revthresh)
        new_thresholds.extend(np.where(np.all(cur_thresh > 0, axis=0), np.max(cur_thresh, axis=0), 0))
        new_thresholds.extend([0])
        new_thresholds_rev.extend(np.where(np.all(cur_revthresh > 0, axis=0), np.max(cur_revthresh, axis=0), 0))
        new_thresholds_rev.extend([0])
        
    ### write output
    # if args.output:
    if not args.output.endswith('.mums'):
        args.output += '.mums'
    with open(args.output, 'w') as f:
        for m in merged:
            f.write('%d\t%s\t%s\n' % (m[0], ','.join(map(str, m[1])), ','.join(['+' if x else '-' for x in m[2]])))
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
        for m in args.mum_files:
            with open(m.replace('.mums', '.lengths'), 'r') as f:
                out.write(f.read().strip() + '\n')
    
if __name__ == "__main__":
    main()
