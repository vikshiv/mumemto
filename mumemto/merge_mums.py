import numpy as np
import argparse
try:
    from utils import get_sequence_lengths, parse_mums_generator
except ImportError:
    from mumemto.utils import get_sequence_lengths, parse_mums_generator
import os
import sys
import subprocess
from tqdm.auto import tqdm
import shutil


def parse_arguments(args=None):  
    parser = argparse.ArgumentParser(description='Merge MUMs files')
    parser.add_argument('--merged_mums', '-m', help='Path to MUMs of MUMs file (only for string merging)')
    parser.add_argument('mum_files', metavar='MUM_FILES', nargs='+', help='Paths to MUMs files to merge')
    parser.add_argument('--output', '-o', help='Path to output merged MUMs file', default='merged.mums')
    parser.add_argument('--verbose', '-v', action='store_true', help='Print verbose output')
    
    if args is None:
        args = parser.parse_args()
    else:
        args = parser.parse_args(args)
    
    assert len(args.mum_files) >= 2, "At least two MUMs files are required for merging"
    
    for i in range(len(args.mum_files)):
        if not args.mum_files[i].endswith('.mums'):
            args.mum_files[i] += '.mums'
            
    if args.merged_mums is not None and not args.merged_mums.endswith('.mums'):
        args.merged_mums += '.mums'
        
    if args.merged_mums is not None and not os.path.exists(args.merged_mums):
        print(f"Error: MUMs of MUMs file {args.merged_mums} does not exist. Omit -m to run merge from start.", file=sys.stderr)
        sys.exit(1)
    return args
    
def merge_anchor_lengths(args):
    length_files = [m.replace('.mums', '.lengths') for m in args.mum_files]
    if not args.output.endswith('.mums'):
        args.output += '.mums'
    out = open(args.output.replace('.mums', '.lengths'), 'w')
    with open(length_files[0], 'r') as f:
        anchor_path = f.readline().split()[0]
    for m in length_files:
        with open(m, 'r') as f:
            first_line = f.readline().split()[0]
            if first_line != anchor_path:
                print(f"Error: Cannot perform anchor-merge. Anchor sequence is not identical in each partition. Ensure paths are identical in the first line of each lengths file.", file=sys.stderr)
                sys.exit(1)
    
    first_file = True
    lines = []
    for m in length_files:
        with open(m, 'r') as f:                
            for l in f.read().splitlines():
                l = l.split()
                if first_file or l[0] != anchor_path:
                    lines.append(l)
        first_file = False
    entry_count = np.array([len(l) for l in lines])
    # all complex or simple lengths, just concatenate them
    if np.all(entry_count == 3) or np.all(entry_count == 2):
        out.write('\n'.join([' '.join(l) for l in lines]))
    # mix of simple and complex, convert simple to complex entries
    else:
        new_lines = []
        for l in lines:
            if len(l) == 3:
                new_lines.append(l)
            else:
                new_lines.append([l[0], '*', l[1]])
                new_lines.append([l[0], os.path.basename(l[0]), l[1]])
        out.write('\n'.join([' '.join(l) for l in new_lines]))
    out.close()
    
def merge_lengths(args):
    if not args.output.endswith('.mums'):
        args.output += '.mums'
    out = open(args.output.replace('.mums', '.lengths'), 'w')
    lines = []
    for m in args.mum_files:
        with open(m.replace('.mums', '.lengths'), 'r') as f:                
            for l in f.read().splitlines():
                l = l.split()
                lines.append(l)
    entry_count = np.array([len(l) for l in lines])
    # all complex or simple lengths, just concatenate them
    if np.all(entry_count == 3) or np.all(entry_count == 2):
        out.write('\n'.join([' '.join(l) for l in lines]))
    # mix of simple and complex, convert simple to complex entries
    else:
        new_lines = []
        for l in lines:
            if len(l) == 3:
                new_lines.append(l)
            else:
                new_lines.append([l[0], '*', l[1]])
                new_lines.append([l[0], os.path.basename(l[0]), l[1]])
        out.write('\n'.join([' '.join(l) for l in new_lines]))
    out.close()
    
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

def run_merger(args):
    script_dir = os.path.dirname(os.path.abspath(__file__))
    if not shutil.which('extract_mums') or not shutil.which('mumemto_exec'):
        extract_script = os.path.join(script_dir, '../extract_mums')
        mumemto_script = os.path.join(script_dir, '../mumemto_exec')
        if not os.path.exists(extract_script) or not os.path.exists(mumemto_script):
            print("Error: mumemto not installed. Please use make to install mumemto first", file=sys.stderr)
            sys.exit(1)
    else:
        extract_script = 'extract_mums'
        mumemto_script = 'mumemto_exec'
    pbar = tqdm(args.mum_files, desc="Extracting MUM sequences", disable=not args.verbose)
    for f in pbar:
        cmd = [extract_script, '-m', f]
        result = subprocess.run(cmd)
        if result.returncode == 1:
            pbar.close()
            print("Error: Partial MUMs detected. Aborting merge. Cleaning up...", file=sys.stderr)
            for f in args.mum_files:
                if os.path.exists(f.replace('.mums', '_mums.fa')):
                    os.remove(f.replace('.mums', '_mums.fa'))
            sys.exit(1)
    
    cmd = [mumemto_script] + [f.replace('.mums', '_mums.fa') for f in args.mum_files] + ['-o', args.output + '_temp_merged']
    if args.verbose:
        print(f"Running command: {' '.join(cmd)}", file=sys.stderr)
    subprocess.run(cmd)
    
    args.merged_mums = args.output + '_temp_merged.mums'

def run_anchor_merger(args):
    if args.verbose:
        print("*.athresh files detected, running anchor merging...", file=sys.stderr)
    script_dir = os.path.dirname(os.path.abspath(__file__))
    anchor_merge_script = 'anchor_merge'
    if not shutil.which(anchor_merge_script):
        anchor_merge_script = os.path.realpath(os.path.join(script_dir, '../anchor_merge'))
    cmd = [anchor_merge_script] + args.mum_files + ['-o', args.output]
    if args.verbose:
        cmd.append('-v')
    if args.verbose:
        print(f"Running command: {' '.join(cmd)}", file=sys.stderr)
    subprocess.run(cmd)

def main(args):
    anchor_merge = all([os.path.exists(m.replace('.mums', '.athresh')) for m in args.mum_files])
    if anchor_merge:
        if args.merged_mums is not None:
            print("Error: -m is only for string merging, but anchor-based merging detected. Ignoring -m.", file=sys.stderr)
        merge_anchor_lengths(args)
        run_anchor_merger(args)
        sys.exit(0)

    threshold_exists = all([os.path.exists(m.replace('.mums', '.thresh')) for m in args.mum_files])
    if not threshold_exists:
        print("Error: *.thresh or *.athresh files required for all inputs for merging.", file=sys.stderr)
        sys.exit(1)
    
    cleanup = args.merged_mums is None
    if args.merged_mums is None:
        run_merger(args)
    
    
    premerge_mums = [list(parse_mums_generator(m)) for m in args.mum_files]
    
    ### get lengths
    mum_lens = get_sequence_lengths(args.merged_mums[:-5] + '.lengths', multilengths=True)

    NUM_SETS = len(mum_lens)
    
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
    
    
    assert len(thresholds) == NUM_SETS, "input # of MUM files does not match merged MUM input file"
    assert len(rev_thresholds) == NUM_SETS, "input # of MUM files does not match merged MUM input file"
    assert len(premerge_mums) == NUM_SETS, "input # of MUM files does not match merged MUM input file"
    assert len(mum_offsets) == NUM_SETS, "input # of MUM files does not match merged MUM input file"
    
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
    for idx, (l, starts, strands) in tqdm(enumerate(dollar_less), total=len(dollar_less), desc="Merging MUMs", disable=not args.verbose):
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
    
    # with open(args.output.replace('.mums', '.lengths'), 'w') as out:
    #     for m in args.mum_files:
    #         with open(m.replace('.mums', '.lengths'), 'r') as f:
    #             out.write(f.read().strip() + '\n')
    merge_lengths(args)
    
    if cleanup:
        for f in args.mum_files:
            os.remove(f.replace('.mums', '_mums.fa'))
        os.remove(args.merged_mums)
        os.remove(args.merged_mums[:-5] + '.lengths')
        
if __name__ == "__main__":
    args = parse_arguments()
    main(args)