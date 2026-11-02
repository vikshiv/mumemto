#!/usr/bin/env python3
"""Contig layout untangling for mumemto viz.

Layout advice = per-genome contig order + orientation (±1).
``untangle_layout`` keeps a reference genome fixed and greedily orders/flips
the others by median mapped position / majority strand vs that reference.

Applying a layout:
  - renumbers contig ids so display order is 0,1,2,...
  - on flipped contigs: strand flip and ``start' = L - start - mum_len``
    (L = contig length), so collinear-block order becomes …3,2,1 with −
"""

from __future__ import annotations

import sys
from dataclasses import dataclass
from typing import List, Optional, Sequence, Tuple

import numpy as np

try:
    from utils import MUMdata, find_coll_blocks
except ImportError:
    from mumemto.utils import MUMdata, find_coll_blocks


@dataclass
class ContigLayout:
    """Per-genome contig order (rank) and orientation.

    order[g] = list of old contig ids left→right
    orient[g][old_cid] = +1 or -1
    """
    order: List[List[int]]
    orient: List[np.ndarray]  # length n_contigs_g, values ±1

    @classmethod
    def identity(cls, multilengths: Sequence[Sequence[int]]) -> "ContigLayout":
        order = [list(range(len(ml))) for ml in multilengths]
        orient = [np.ones(len(ml), dtype=np.int8) for ml in multilengths]
        return cls(order=order, orient=orient)

    def multilengths_reordered(self, multilengths: Sequence[Sequence[int]]) -> List[List[int]]:
        return [[int(ml[c]) for c in self.order[g]] for g, ml in enumerate(multilengths)]

    def contig_names_reordered(
        self, contig_names: Sequence[Sequence[str]]
    ) -> List[List[str]]:
        return [[contig_names[g][c] for c in self.order[g]] for g in range(len(self.order))]


def contig_offsets(order: Sequence[int], lengths: Sequence[int]) -> np.ndarray:
    """Left edge of each old contig id when concatenated in ``order``."""
    left = np.zeros(len(lengths), dtype=np.int64)
    pos = 0
    for cid in order:
        left[cid] = pos
        pos += int(lengths[cid])
    return left


def apply_contig_layout(
    mums: MUMdata,
    multilengths: Sequence[Sequence[int]],
    layout: ContigLayout,
) -> Tuple[MUMdata, List[List[int]]]:
    """Apply order/orientation advice to contig-coord MUMs.

    - Remap contig ids so display order is 0,1,2,...
    - On flipped contigs: strand := not strand, and
      ``start' = L - start - mum_len`` (still contig-local).

    That coord flip puts MUMs in reverse order along the contig so
    ``find_coll_blocks`` still sees a valid reverse-strand run (…3,2,1 with −).
    """
    assert mums.contig_ids is not None
    n_seqs = mums.num_seqs
    new_cids = mums.contig_ids.copy()
    new_starts = mums.starts.copy()
    new_strands = mums.strands.copy()

    for g in range(n_seqs):
        lens = np.asarray(multilengths[g], dtype=np.int64)
        # old_cid -> new_cid (rank index)
        old_to_new = np.empty(len(layout.order[g]), dtype=np.int32)
        for new_id, old_cid in enumerate(layout.order[g]):
            old_to_new[old_cid] = new_id

        present = mums.contig_ids[:, g] >= 0
        old = mums.contig_ids[present, g]
        new_cids[present, g] = old_to_new[old]

        flipped = layout.orient[g] < 0
        flip_mask = np.zeros(mums.num_mums, dtype=bool)
        flip_mask[present] = flipped[old]
        if np.any(flip_mask):
            cid_f = mums.contig_ids[flip_mask, g]
            x = mums.starts[flip_mask, g].astype(np.int64)
            ml = mums.lengths[flip_mask].astype(np.int64)
            new_starts[flip_mask, g] = lens[cid_f] - x - ml
            new_strands[flip_mask, g] = ~mums.strands[flip_mask, g]

    out = MUMdata.from_arrays(
        mums.lengths.copy(),
        new_starts,
        new_strands,
        blocks=None,
        contig_ids=new_cids,
    )
    return out, layout.multilengths_reordered(multilengths)


def gapped_genome_length(contig_lengths: Sequence[int], spacer: int = 0) -> int:
    """Display length of one genome with ``spacer`` bp between contigs."""
    n = len(contig_lengths)
    if n == 0:
        return 0
    return int(sum(contig_lengths)) + int(spacer) * max(0, n - 1)


def trim_trailing_empty_contigs(
    mums: MUMdata,
    multilengths: Sequence[Sequence[int]],
    verbose: bool = False,
    min_frac: float = 1e-4,
) -> Tuple[MUMdata, List[List[int]]]:
    """Drop trailing contigs with little/no MUM evidence per genome.

    A trailing contig is removed while its total MUM bp is 0, or below
    ``min_frac`` of that genome's total MUM bp (default 0.01%). Mid-genome
    sparse contigs are kept. Contig ids for the kept prefix stay 0..keep_n-1.
    If a genome has no MUM hits at all, its contig list is left unchanged.
    """
    assert mums.contig_ids is not None
    new_mls: List[List[int]] = []
    new_cids = mums.contig_ids.copy()
    for g, ml in enumerate(multilengths):
        n = len(ml)
        if n == 0:
            new_mls.append([])
            continue
        present = mums.contig_ids[:, g] >= 0
        covered = np.zeros(n, dtype=np.int64)
        if np.any(present):
            c = mums.contig_ids[present, g]
            # sum of MUM lengths per contig
            covered = np.bincount(
                c, weights=mums.lengths[present].astype(np.float64), minlength=n,
            )[:n].astype(np.int64)
        total = int(covered.sum())
        thresh = int(min_frac * total) if total > 0 else 0
        last = n - 1
        # Drop trailing contigs that are empty/sparse in absolute MUM bp or
        # cover <2% of their own length (typical unplaced junk at FASTA end).
        while last > 0:
            L = int(ml[last])
            frac = (covered[last] / L) if L > 0 else 0.0
            if covered[last] <= thresh or frac < 0.02:
                last -= 1
            else:
                break
        if last < 0:
            keep_n = n
        else:
            keep_n = last + 1
        dropped = n - keep_n
        if verbose and dropped:
            print(
                f'trim genome {g}: drop {dropped} trailing empty/sparse contigs '
                f'({n} → {keep_n}, thresh={thresh} bp)',
                file=sys.stderr,
            )
        new_mls.append([int(x) for x in ml[:keep_n]])
        bad = present & (mums.contig_ids[:, g] >= keep_n)
        if np.any(bad):
            new_cids[bad, g] = -1
    out = MUMdata.from_arrays(
        mums.lengths.copy(),
        mums.starts.copy(),
        mums.strands.copy(),
        blocks=None,
        contig_ids=new_cids,
    )
    return out, new_mls


def contig_local_to_concat(
    mums: MUMdata,
    multilengths: Sequence[Sequence[int]],
    keep_contig_ids: bool = False,
    spacer: int = 0,
) -> MUMdata:
    """Convert contig-local starts to concat-global for plotting (ids 0..n-1 order).

    Does not flip coordinates; orientation is already encoded in strands via
    ``apply_contig_layout``.

    If ``spacer`` > 0, insert that many bp of empty space between consecutive
    contigs (``--mode spaced``). Block finding should use spacer=0.
    """
    assert mums.contig_ids is not None
    spacer = int(spacer)
    new_starts = mums.starts.copy()
    for g, ml in enumerate(multilengths):
        left = np.zeros(len(ml), dtype=np.int64)
        pos = 0
        for c, L in enumerate(ml):
            left[c] = pos
            pos += int(L)
            if c < len(ml) - 1:
                pos += spacer
        present = (mums.starts[:, g] >= 0) & (mums.contig_ids[:, g] >= 0)
        cid = mums.contig_ids[present, g]
        new_starts[present, g] = left[cid] + mums.starts[present, g]
    return MUMdata.from_arrays(
        mums.lengths.copy(),
        new_starts,
        mums.strands.copy(),
        blocks=None,
        contig_ids=(mums.contig_ids.copy() if keep_contig_ids else None),
    )


def reorder_mums(mums: MUMdata, order: np.ndarray) -> MUMdata:
    """Return MUMdata with rows permuted by ``order`` (keeps contig_ids)."""
    cids = None if mums.contig_ids is None else mums.contig_ids[order]
    return MUMdata.from_arrays(
        mums.lengths[order],
        mums.starts[order],
        mums.strands[order],
        contig_ids=cids,
    )


def split_blocks_at_contig_breaks(blocks: list, contig_ids: np.ndarray) -> list:
    """Split collinear blocks wherever contig id changes on any genome."""
    out = []
    for l, r in blocks:
        start = l
        for i in range(l, r):
            if not np.array_equal(contig_ids[i], contig_ids[i + 1]):
                if start <= i:
                    out.append((start, i))
                start = i + 1
        if start <= r:
            out.append((start, r))
    return out


def prepare_contig_plot(
    mums: MUMdata,
    multilengths: Sequence[Sequence[int]],
    layout: ContigLayout,
    max_break: Optional[int] = None,
    verbose: bool = False,
    spacer: int = 0,
) -> Tuple[MUMdata, list, List[List[int]]]:
    """Apply layout advice and return concat MUMs + blocks for viz.

    1. Find collinear blocks on identity concat (before flips).
    2. Split at contig boundaries (so reorder cannot stretch a ribbon).
    3. Apply layout: remap contig ids; on flips
       ``start' = L - start - mum_len`` (L = contig length) and strand flip.
    4. Concat with new contig order for plotting (optional inter-contig spacer).

    Blocks are kept from step 1 so differential flips across genomes cannot
    drop MUMs from the plot; the L−offset transform keeps geometry correct.
    """
    assert mums.contig_ids is not None
    if max_break is None:
        max_break = 0

    # --- blocks on identity (no spacer; adjacency is biological) ---
    ident_concat = contig_local_to_concat(mums, multilengths, keep_contig_ids=True)
    order = np.argsort(ident_concat.starts[:, 0], kind='mergesort')
    mums_s = reorder_mums(mums, order)
    ident_sorted = reorder_mums(ident_concat, order)
    blocks = find_coll_blocks(ident_sorted, max_break=max_break, verbose=verbose)
    blocks = split_blocks_at_contig_breaks(blocks, mums_s.contig_ids)
    if verbose:
        print(f'collinear blocks after contig splits: {len(blocks)}', file=sys.stderr)

    # --- apply order/orient (ids + L-offset + strand) ---
    mums_s, ml_new = apply_contig_layout(mums_s, multilengths, layout)
    plot_mums = contig_local_to_concat(
        mums_s, ml_new, keep_contig_ids=False, spacer=spacer,
    )
    plot_mums.blocks = blocks
    return plot_mums, blocks, ml_new


# ---- Untangle (greedy, reference-fixed) ------------------------------------

def untangle_layout(
    mums: MUMdata,
    multilengths: Sequence[Sequence[int]],
    ref: int = 0,
    verbose: bool = False,
) -> ContigLayout:
    """Order/flip contigs by median mapped position / majority strand vs reference.

    Keeps ``ref`` (default genome 0) in identity order/orientation; only other
    genomes are reordered and flipped.
    """
    assert mums.contig_ids is not None
    layout = ContigLayout.identity(multilengths)
    ref_left = contig_offsets(layout.order[ref], multilengths[ref])

    ref_cid = mums.contig_ids[:, ref]
    ref_x = mums.starts[:, ref].astype(np.int64)
    ref_present = (ref_x >= 0) & (ref_cid >= 0)
    ref_mid = np.full(mums.num_mums, np.nan)
    if np.any(ref_present):
        cid = ref_cid[ref_present]
        x = ref_x[ref_present]
        ml = mums.lengths[ref_present].astype(np.int64)
        ref_mid[ref_present] = ref_left[cid] + x + ml / 2.0

    for g in range(mums.num_seqs):
        if g == ref:
            continue
        n_ctg = len(multilengths[g])
        scores = [[] for _ in range(n_ctg)]
        strand_votes = np.zeros(n_ctg, dtype=np.int64)
        for i in range(mums.num_mums):
            cg = mums.contig_ids[i, g]
            if cg < 0 or mums.starts[i, g] < 0 or not ref_present[i]:
                continue
            scores[cg].append(ref_mid[i])
            same = mums.strands[i, g] == mums.strands[i, ref]
            strand_votes[cg] += 1 if same else -1

        medians = []
        for c in range(n_ctg):
            if scores[c]:
                medians.append((float(np.median(scores[c])), c))
            else:
                medians.append((float('inf'), c))
        medians.sort()
        new_order = [c for _, c in medians]
        new_orient = np.ones(n_ctg, dtype=np.int8)
        for c in range(n_ctg):
            if strand_votes[c] < 0:
                new_orient[c] = -1
        layout.order[g] = new_order
        layout.orient[g] = new_orient
        if verbose:
            n_flip = int(np.sum(new_orient < 0))
            print(f'untangle genome {g}: order={new_order} flips={n_flip}', file=sys.stderr)
    return layout
