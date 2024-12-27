import numpy as np

def find_coll_blocks(mums, max_length, dpi=500, size=(6.4, 4.8), max_break=None, verbose=False):
    if max_break is None:
        bp_per_inch = max_length / (dpi * size[0])
        max_break = min(bp_per_inch, 100000)
    if verbose:
        print('max gap within a collinear block:', max_break)
    starts = np.array([m[1] for m in mums])
    mum_orders = starts.transpose().argsort()
    mum_gaps = []
    flips = set([])
    for i in range(mum_orders.shape[0]):
        cur = []
        for l in range(1, mum_orders.shape[1]):
            left, right = mum_orders[i][l-1], mum_orders[i][l]
            if mums[left][2][i] == mums[right][2][i]:
                if mums[left][2][i] == '+':
                    cur.append((left, right))
                else:
                    cur.append((right, left))
                    flips.add((right, left))
        mum_gaps.append(cur)
    common_gaps = set.intersection(*map(set, mum_gaps))
    left, right = zip(*common_gaps)
    left, right = set(list(left)), set(list(right))
    true_collinear_mums = sorted(list(left.intersection(right)))
    right_coll_mums = sorted(list(left.difference(set(true_collinear_mums)))) # have a right pair, but not a left
    left_coll_mums = sorted(list(right.difference(set(true_collinear_mums)))) # have a left pair, but not a right
    large_blocks = list(zip(right_coll_mums, left_coll_mums))
    ### find the longest stretches of collinear mums
    small_blocks = []
    for l, r in large_blocks:
        last = l
        for i in range(l, r):
            lens = np.full(len(mums[i][1]), mums[i][0])
            lens[(mums[i+1][1] < mums[i][1])] = mums[i+1][0] 
            gap_lens = np.abs(mums[i][1] - mums[i+1][1]) - lens
            if gap_lens.max() > max_break and last < i:
                small_blocks.append((last, i))
                last = i + 1
        if last != r:
            small_blocks.append((last, r))
    return large_blocks, small_blocks, mum_gaps

def get_block_order(mums, blocks):
    starts = np.array([m[1] for m in mums])
    mum_orders = starts.transpose().argsort()
    ### get coll_block order
    coll_block_starts = [b[0] for b in blocks]
    coll_block_orders = []
    for i in range(mum_orders.shape[0]):
        poi = np.where(np.isin(mum_orders[i], coll_block_starts))[0]
        values = mum_orders[i][poi]
        poi = np.argsort(values)
        coll_block_orders.append(np.argsort(poi))
    return coll_block_orders

def parse_mums(mumfile, seq_lengths, lenfilter=0, subsample=1):
    def reverse_strand(l, starts, strands):
        new_starts = np.array([p if s == '+' or s == '' else seq_lengths[idx] - p - l for idx, (p, s) in enumerate(zip(starts, strands))])
        return (l, new_starts, strands)
    count = 0
    for l in open(mumfile, 'r').readlines():
        if count % subsample == 0:
            l = l.strip().split()
            if int(l[0]) >= lenfilter:
                yield reverse_strand(int(l[0]), [int(v) if v else -1 for v in l[1].split(',')], tuple(l[2].split(',')))
        count += 1