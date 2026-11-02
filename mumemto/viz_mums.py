#!/usr/bin/env python3

from matplotlib import pyplot as plt
from matplotlib import cm
from matplotlib.collections import PolyCollection
from mpl_toolkits.axes_grid1 import make_axes_locatable
import os, sys
import numpy as np
import argparse
try:
    from utils import (
        MUMdata, find_coll_blocks, get_sequence_lengths, get_seq_paths,
        ensure_contig_coords,
    )
    from untangle import (
        ContigLayout, untangle_layout, prepare_contig_plot,
        gapped_genome_length, trim_trailing_empty_contigs,
    )
except ImportError:
    from mumemto.utils import (
        MUMdata, find_coll_blocks, get_sequence_lengths, get_seq_paths,
        ensure_contig_coords,
    )
    from mumemto.untangle import (
        ContigLayout, untangle_layout, prepare_contig_plot,
        gapped_genome_length, trim_trailing_empty_contigs,
    )

def parse_arguments(args=None):    
    parser = argparse.ArgumentParser(description="Plots a synteny plot of MUMs from mumemto")
    # parser.add_argument('--filelist', '-f', dest='filelist', help='path to filelist from mumemto', required=True)
    # parser.add_argument('--mums', '-m', dest='mumfile', help='path to *.mum file from mumemto', required=True)
    # parser.add_argument('--lengths','-l', dest='lens', help='lengths file, first column is seq length in order of filelist', required=True)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--input-prefix', '-i', dest='prefix', help='prefix for filelist, mums, and lengths files (uses *.mums file over *.bumbl file if both exist)')
    group.add_argument('--mums', '-m', dest='mumfile', help='path to *.mums (or *.bumbl) file from mumemto')
    
    parser.add_argument('--lengths','-l', dest='lens', help='lengths file, first column is seq length in order of filelist')
    
    parser.add_argument('--filelist', '-f', dest='filelist', help='filelist defines the order of seqs in the plot. First col = full path from lengths file')
    parser.add_argument('--labels', '-n', dest='labels', help='labels for each sequence in the plot, one a line, matching the order of the filelist (if given) or the lengths file')
    parser.add_argument('--len-filter','-L', dest='lenfilter', help='only plot MUMs longer than threshold', default=0, type=int)
    parser.add_argument('--subsample','-s', dest='subsample', help='subsample every Nth mum', default=1, type=int)
    parser.add_argument('--center','-c', dest='center', action='store_true', help='center plot', default=False)
    parser.add_argument('--inversion-color','-ic', dest='inv_color', help='color for inversions', default='green')
    parser.add_argument('--mum-color','-mc', dest='mum_color', help='color as hex (default: #00A2FF)', default='#00A2FF', type=str)
    parser.add_argument('--alpha','-a', dest='alpha', help='opacity of mums [0-1] (default: 0.8 for blocks, 0.1 for MUMs)', type=float)
    parser.add_argument('--linewidth','-lw', dest='linewidth', help='linewidth of mums (default: 0 for blocks, 0.05 for MUMs)', type=float)
    parser.add_argument('--fout','-o', dest='filename', help='plot fname (default: input_prefix)')
    parser.add_argument('--dims', dest='size', help='fig dimensions (inches) (default: 6.4, 4.8)', default=(6.4,4.8), type=float, nargs=2)
    parser.add_argument('--region', dest='region', help='genomic region to zoom in, global coordinates (assuming contigs are concatenated)', default=None, type=float, nargs=2)
    parser.add_argument('--dpi','-d', dest='dpi', help='dpi', default=500, type=int)
    parser.add_argument('--verbose','-v', dest='verbose', help='verbose mode', action='store_true', default=False)
    parser.add_argument('--no-coll-block','-b', dest='no_coll_block', help='plot only MUMs, not collinear blocks (slower) (default: false)', action='store_true', default=False)
    parser.add_argument('--max-gap-len','-g', dest='max_break', help='maximum break between collinear mums within a collinear block (default: <1px)', default=None, type=int)
    
    # handle multi-fasta inputs
    parser.add_argument('--mode', dest='mode', choices=['normal', 'delineated', 'gapped', 'spaced'],
                       help='Plotting mode for multi-fasta inputs: normal=concat; '
                            'delineated=concat with red contig ticks; '
                            'gapped=column-aligned gaps (equal contig counts); '
                            'spaced=per-genome contig gaps + scale bar (default with --untangle)',
                       default='normal')
    parser.add_argument(
        '--spacer', dest='spacer', type=float, default=0.1,
        help='Inter-contig gap for --mode spaced/gapped, as a fraction of the largest '
             'contig length (default: 0.1). Smaller = tighter; larger = more whitespace.',
    )

    # Contig layout untangling (requires contig-local coords or multilengths to convert)
    parser.add_argument(
        '--untangle', dest='untangle', action='store_true', default=False,
        help='Reorder/orient contigs vs genome 0 (keeps ref fixed; others by median '
             'mapped position / majority strand). Implies --mode spaced. '
             'Requires contig coords or multilengths.',
    )
    empty_ctg = parser.add_mutually_exclusive_group()
    empty_ctg.add_argument(
        '--hide-empty-contigs', dest='hide_empty_contigs', action='store_true',
        default=None,
        help='Drop trailing contigs with no MUM hits (default: on for --mode spaced).',
    )
    empty_ctg.add_argument(
        '--show-empty-contigs', dest='hide_empty_contigs', action='store_false',
        help='Keep trailing contigs even if they have no MUM hits.',
    )

    if args is None:
        args = parser.parse_args()
    else:
        args = parser.parse_args(args)
    
    if args.mumfile:
        # Allow prefix.mums from decompress; strip .mums / .bumbl only
        base, ext = os.path.splitext(args.mumfile)
        if ext in ('.mums', '.bumbl'):
            args.prefix = base
        else:
            args.prefix = args.mumfile
    elif args.prefix:
        if args.prefix.endswith('.mums') or args.prefix.endswith('.bumbl'):
            args.prefix = os.path.splitext(args.prefix)[0]
        if os.path.exists(args.prefix + '.mums'):
            args.mumfile = args.prefix + '.mums'
        elif os.path.exists(args.prefix + '.bumbl'):
            args.mumfile = args.prefix + '.bumbl' 
        else:
            parser.error("No *.mums or *.bumbl file found for prefix %s" % args.prefix)
    else:
        parser.error("Either --mums or --prefix must be provided")
        
    if args.lens is None:
        args.lens = args.prefix + '.lengths'
        if not os.path.exists(args.lens):
            raise FileNotFoundError(f"Lengths file {args.lens} not found, and no lengths file provided")
    
    if not args.filename:
        args.filename = args.prefix
        
    # Set defaults based on no_coll_block flag
    if args.alpha is None:
        args.alpha = 0.05 if args.no_coll_block else 0.8
    if args.linewidth is None:
        args.linewidth = 0.05 if args.no_coll_block else 0

    # Untangle implies spaced display (not red contig ticks)
    if args.untangle and args.mode == 'normal':
        args.mode = 'spaced'

    # Trailing empty-contig trim defaults on for spaced mode
    if args.hide_empty_contigs is None:
        args.hide_empty_contigs = (args.mode == 'spaced')

    # if (args.mode == 'gapped' or args.mode == 'delineated') and args.multilengths is None:
    #     parser.error('multi-lengths file is required for gapped or delineated mode')
    return args


def _nice_scale_bp(visible_span_bp):
    """Pick a round scale length from the visible x-span (Google Maps–style).

    Chooses a 1/2/5 × 10^n length that is ~8–15% of ``visible_span_bp``,
    so zooming (``--region``) or shorter assemblies automatically steps
    through bp → kb → Mb → Gb.
    """
    span = float(visible_span_bp)
    if span <= 0:
        return 1, '1 bp'
    candidates = [
        (1, '1 bp'), (2, '2 bp'), (5, '5 bp'),
        (10, '10 bp'), (20, '20 bp'), (50, '50 bp'),
        (100, '100 bp'), (200, '200 bp'), (500, '500 bp'),
        (1e3, '1 kb'), (2e3, '2 kb'), (5e3, '5 kb'),
        (1e4, '10 kb'), (2e4, '20 kb'), (5e4, '50 kb'),
        (1e5, '100 kb'), (2e5, '200 kb'), (5e5, '500 kb'),
        (1e6, '1 Mb'), (2e6, '2 Mb'), (5e6, '5 Mb'),
        (1e7, '10 Mb'), (2e7, '20 Mb'), (5e7, '50 Mb'),
        (1e8, '100 Mb'), (2e8, '200 Mb'), (5e8, '500 Mb'),
        (1e9, '1 Gb'), (2e9, '2 Gb'), (5e9, '5 Gb'),
    ]
    target = span * 0.12
    # Prefer a bar that is ~5–20% of the visible window
    fits = [c for c in candidates if 0.05 * span <= c[0] <= 0.22 * span]
    pool = fits if fits else candidates
    best = min(pool, key=lambda c: abs(c[0] - target))
    return int(best[0]), best[1]


def add_scale_bar(ax, color='black'):
    """Draw a genomic scale bar in a dedicated strip under the plot.

    Appends a thin axes below ``ax`` via ``make_axes_locatable`` so the bar:
      - never overlaps ribbons,
      - never sits in far-right empty data space,
      - never expands the figure sideways,
      - keeps a correct bp length for the current ``xlim`` (incl. ``--region``).

    Bar is drawn at the left of the strip (under assembly starts). Call after
    ``xlim`` is set; skip ``tight_layout`` afterward (it fights the divider).
    """
    x0, x1 = ax.get_xlim()
    span = abs(x1 - x0)
    if span <= 0:
        return
    length_bp, label = _nice_scale_bp(span)

    divider = make_axes_locatable(ax)
    sax = divider.append_axes('bottom', size='9%', pad=0.12)
    sax.set_xlim(x0, x1)
    sax.set_ylim(0, 1)
    sax.set_xticks([])
    sax.set_yticks([])
    for spine in sax.spines.values():
        spine.set_visible(False)
    sax.patch.set_visible(False)

    left = x0 + 0.01 * span
    right = left + length_bp
    y = 0.62
    tick = 0.28
    sax.plot([left, right], [y, y], color=color, linewidth=1.8,
             solid_capstyle='butt', clip_on=False)
    sax.plot([left, left], [y - tick, y + tick], color=color, linewidth=1.3,
             clip_on=False)
    sax.plot([right, right], [y - tick, y + tick], color=color, linewidth=1.3,
             clip_on=False)
    sax.text((left + right) / 2, 0.08, label, ha='center', va='bottom',
             fontsize=8, color=color, clip_on=False)

def points_to_poly(points):
    starts, ends = tuple(zip(*points))
    points = starts + ends[::-1]
    return points
    
def get_mum_polygons(mums, centering, color='#00A2FF', inv_color='red'):
    polygons = []
    colors = []    
    for (l, starts, strands) in mums:
        inverted = not strands[0]
        points = []
        for idx, (x, strand) in enumerate(zip(starts, strands)):
            if x == -1:
                if len(points) > 2:
                    polygons.append(points_to_poly(points[:-1]))
                    colors.append(color)
                continue
            points.append(((centering[idx] + x, idx), (centering[idx] + x + l, idx)))
            if not inverted and not strand:
                inverted = True
                if len(points) > 2:
                    polygons.append(points_to_poly(points[:-1]))
                    colors.append(color)
                polygons.append(points_to_poly(points[-2:]))
                colors.append(inv_color)
                points = [points[-1]]
            elif inverted and strand:
                inverted = False
                if len(points) > 2:
                    polygons.append(points_to_poly(points[:-1]))
                    colors.append(color)
                polygons.append(points_to_poly(points[-2:]))
                colors.append(inv_color)
                points = [points[-1]]
        if len(points) >= 2:
            polygons.append(points_to_poly(points))
            colors.append(color)
    return polygons, colors

def get_block_polygons(collinear_blocks, mums, centering, color='#00A2FF', inv_color='red'):
    polygons = []
    colors = []
    for _, (l, r) in enumerate(collinear_blocks):
        strands = mums[l].strands
        inverted = not strands[0]
        points = []
        left, right = mums[l].starts.copy(), (mums[r].starts + mums[r].length).copy()
        # After contig flips, identity-ordered endpoints can reverse on a genome
        for idx in range(len(left)):
            if left[idx] != -1 and right[idx] != -1 and left[idx] > right[idx]:
                left[idx], right[idx] = right[idx], left[idx]
        for idx, strand in enumerate(strands):
            if left[idx] == -1:
                if len(points) > 2:
                    polygons.append(points_to_poly(points))
                    colors.append(color)
                points = []
                continue
            points.append(((centering[idx] + left[idx], idx), (centering[idx] + right[idx], idx)))
            if not inverted and not strand:
                inverted = True
                if len(points) > 2:
                    polygons.append(points_to_poly(points[:-1]))
                    colors.append(color)
                polygons.append(points_to_poly(points[-2:]))
                colors.append(inv_color)
                points = [points[-1]]
            elif inverted and strand:
                inverted = False
                if len(points) > 2:
                    polygons.append(points_to_poly(points[:-1]))
                    colors.append(color)
                polygons.append(points_to_poly(points[-2:]))
                colors.append(inv_color)
                points = [points[-1]]
        if len(points) >= 2:
            polygons.append(points_to_poly(points))
            colors.append(color)
    return polygons, colors

def _polygons_xmax(polygons):
    """Rightmost x among block/MUM polygons (display coords)."""
    mx = 0.0
    for poly in polygons:
        for x, _y in poly:
            if x > mx:
                mx = float(x)
    return mx


def mum_dense_xmax(mums, quantile: float = 0.99):
    """Rightmost display-x covering ``quantile`` of MUM midpoints.

    Using the raw max lets a few sparse MUMs on junk contigs stretch xlim and
    leave a blank right margin; a high quantile tracks the dense ribbon mass.
    """
    mids = []
    for g in range(mums.num_seqs):
        present = mums.starts[:, g] >= 0
        if not np.any(present):
            continue
        s = mums.starts[present, g].astype(np.float64)
        mids.append(s + 0.5 * mums.lengths[present])
    if not mids:
        return 0.0
    mids = np.concatenate(mids)
    return float(np.quantile(mids, quantile))


def plot(args, genome_lengths, polygons, colors, centering, dpi=500, size=None, genomes=None, filename=None):
    fig, ax = plt.subplots()
    if size:
        fig.set_size_inches(*size)
    max_length = max(genome_lengths) if genome_lengths else 0
    spacer_bp = int(getattr(args, 'spacer_bp', 0) or 0)

    # Spaced: crop right margin to dense MUM content (not full backbone / junk tails)
    if args.mode == 'spaced' and args.region is None:
        content_xmax = getattr(args, 'plot_xmax', None)
        if content_xmax is None and polygons:
            content_xmax = _polygons_xmax(polygons)
        if content_xmax and content_xmax > 0:
            pad = max(1.0, 0.005 * content_xmax)
            max_length = min(max_length, float(content_xmax) + pad)

    # Plot genome lines based on mode
    if args.mode == 'normal' or args.mode == 'delineated':
        # Just plot simple genome lines
        for idx, g in enumerate(genome_lengths):
            ax.plot([centering[idx] + 0, centering[idx] + g], [idx, idx], 
                    alpha=0.2, linewidth=0.75, c='black')

    elif args.mode == 'spaced':
        # Per-genome contig segments with whitespace gaps (no red ticks).
        # Clip segments to max_length so tails past MUM content are not drawn.
        for idx, ml in enumerate(args.multilengths):
            pos = centering[idx]
            for i, L in enumerate(ml):
                L = int(L)
                if pos >= max_length:
                    break
                if L > 0:
                    x1 = min(pos + L, max_length)
                    if x1 > pos:
                        ax.plot([pos, x1], [idx, idx],
                                alpha=0.35, linewidth=1.0, c='black',
                                solid_capstyle='butt')
                pos += L
                if i < len(ml) - 1:
                    pos += spacer_bp

    elif args.mode == 'gapped':
        # Plot with gaps between subsequences
        offsets = np.array([0] + (args.multilengths.max(axis=0) + args.spacer).cumsum().tolist()[:-1])
        vert_seps = [p - (args.spacer) / 2 for p in offsets] + [args.multilengths.max(axis=0).sum() + args.spacer * (args.multilengths.shape[1] - 1)]
        for p in vert_seps[1:-1]:
            ax.plot([p, p], [0, len(genome_lengths)-1], alpha=0.5, linewidth=1, color='black')
        for idx in range(args.multilengths.shape[0]):
            for i, offset in enumerate(args.multilengths[idx]):
                ax.plot([centering[idx] + offsets[i], centering[idx] + offsets[i] + offset], 
                        [idx, idx], alpha=0.2, linewidth=0.25)
        chr_markers = [vert_seps[idx - 1] + ((vert_seps[idx] - vert_seps[idx-1]) / 2) for idx in range(1, len(vert_seps))]
        ax.set_xticks(chr_markers)
        ax.set_xticklabels(range(1, len(chr_markers) + 1))
    else:
        # Default behavior for non-multifasta inputs, should never be reached
        for idx, g in enumerate(genome_lengths):
            ax.plot([centering[idx] + 0, centering[idx] + g], [idx, idx], 
                   alpha=0.2, linewidth=0.75)

    ax.add_collection(PolyCollection(polygons, linewidths=args.linewidth, alpha=args.alpha, edgecolors=colors, facecolors=colors))
    
    
    if args.mode == 'delineated':
        # Plot lines with delineators for multifasta
        for idx in range(len(args.multilengths)):
            cur_offsets = np.cumsum(args.multilengths[idx])
            last_offset = 0
            for i, offset in enumerate([0] + cur_offsets[:-1]):
                ax.plot([centering[idx] + offset, centering[idx] + offset], 
                        [idx - 0.25, idx + 0.25], alpha=1, linewidth=0.25, color='red')
                last_offset = offset
            
            
    ax.yaxis.set_ticks(list(range(len(genome_lengths))))
    ax.tick_params(axis='y', which='both',length=0)
    if genomes:
        ax.set_yticklabels(genomes, fontsize=8)
    else:
        ax.yaxis.set_ticklabels([])
    if args.mode == 'gapped':
        ax.set_xlabel('chromosome')
    elif args.mode == 'spaced':
        ax.set_xlabel('')
        ax.set_xticks([])
        ax.tick_params(axis='x', which='both', length=0, labelbottom=False)
        ax.spines['bottom'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
    else:
        ax.set_xlabel('genomic position')
    ax.set_ylabel('sequences')
    ax.set_ylim(-0.25, len(genome_lengths) - 1 + 0.25)
    if args.mode == 'gapped':
        ax.set_xlim(0, args.multilengths.max(axis=0).sum() + args.spacer * (args.multilengths.shape[1] - 1))
    elif args.region is not None:
        ax.set_xlim(*args.region)
    else:
        ax.set_xlim(0, max_length)
    ax.invert_yaxis()
    if args.mode == 'spaced':
        # After xlim is final; divider strip must not be followed by tight_layout
        add_scale_bar(ax)
    else:
        fig.set_tight_layout(True)

    if filename:
        if os.path.splitext(filename)[1] not in ['.png', '.pdf', '.svg']:
            filename = filename + '.png'
        if os.path.dirname(filename):
            path = filename
        else:
            path = os.path.join(os.path.dirname(args.mumfile), filename)
        fig.savefig(path, dpi=dpi, bbox_inches='tight', pad_inches=0.25)
    return ax

def offset_mums(args, mums, spacer=None, blocks=None):
    offset = args.multilengths
    if spacer is None:
        spacer = args.spacer
    NUM_SEQS = len(offset)
    ### offset the mums by the contig locations
    offsets = np.cumsum(offset, axis=1)
    ### label each mum with the contig it belongs to
    # breaks = np.hstack((np.array([[0]*NUM_SEQS]).transpose(), offsets))
    contig_idx = np.array([np.searchsorted(offsets[idx], mums.starts[:,idx], side='right') for idx in range(NUM_SEQS)]).transpose()
    # if using blocks, split any that cross contigs
    if blocks is not None:
        new_blocks = []
        for l,r in blocks:
            # if block crosses contig boundaries in any sequence, break it up
            if not np.all(contig_idx[l,:] == contig_idx[r,:]):
                if l == r - 1:
                    continue
                diffs = ~(np.diff(contig_idx[l:r+1,:], axis=0) == 0).all(axis=1)
                old_l = l
                for s in np.where(diffs)[0] + 1:
                    if s == 1:
                        l = old_l + 2
                        continue
                    cur = old_l + s
                    if l < cur - 1 and cur - 1 <= r:
                        new_blocks.append((l, cur - 1))
                        l = cur + 1
                if l < r:
                    new_blocks.append((l, r))
            else:
                new_blocks.append((l, r))
    else:
        new_blocks = None   
    ### get the relative offset of each mum to the start of its contig
    left_start = np.hstack((np.zeros((offsets.shape[0],1), dtype=int), offsets[:,:-1]))
    rel_offsets = mums.starts - left_start[np.arange(NUM_SEQS), contig_idx]
    partial_mask = mums.starts != -1
    new_starts = np.array([0] + (offset.max(axis=0) + spacer).cumsum().tolist()[:-1])[contig_idx] + rel_offsets
    mums.starts[partial_mask] = new_starts[partial_mask]
    return new_blocks


def insert_spaced_spacers(mums, multilengths, spacer, blocks=None):
    """Insert per-genome inter-contig spacers into concat-global MUM starts.

    For genome g, contig k: display_start = concat_start + k * spacer.
    Splits collinear blocks that cross contig boundaries when ``blocks`` given.
    """
    spacer = int(spacer)
    n_seqs = len(multilengths)
    contig_idx = np.zeros((mums.num_mums, n_seqs), dtype=np.int32)
    for g, ml in enumerate(multilengths):
        offsets = np.cumsum(ml)
        starts = mums.starts[:, g]
        present = starts >= 0
        cid = np.zeros(mums.num_mums, dtype=np.int32)
        cid[present] = np.searchsorted(offsets, starts[present], side='right')
        # Clamp past-end to last contig
        cid[present] = np.minimum(cid[present], len(ml) - 1)
        contig_idx[:, g] = cid
        if spacer:
            mums.starts[present, g] = starts[present] + cid[present] * spacer

    if blocks is None:
        return None
    new_blocks = []
    for l, r in blocks:
        if not np.all(contig_idx[l, :] == contig_idx[r, :]):
            if l == r - 1:
                continue
            diffs = ~(np.diff(contig_idx[l:r + 1, :], axis=0) == 0).all(axis=1)
            old_l = l
            cur_l = l
            for s in np.where(diffs)[0] + 1:
                if s == 1:
                    cur_l = old_l + 2
                    continue
                cur = old_l + s
                if cur_l < cur - 1 and cur - 1 <= r:
                    new_blocks.append((cur_l, cur - 1))
                    cur_l = cur + 1
            if cur_l < r:
                new_blocks.append((cur_l, r))
        else:
            new_blocks.append((l, r))
    return new_blocks

def main(args):
    need_multi = (
        args.mode != 'normal'
        or args.untangle
        or args.hide_empty_contigs
    )
    if need_multi:
        try:
            offset = get_sequence_lengths(args.lens, multilengths=True)
            seq_lengths = [sum(o) for o in offset]
            if args.mode == 'gapped' and len(set([len(o) for o in offset])) > 1:
                print('Warning: gapped mode requires the same number of sequences per input FASTA file. Using spaced mode instead.', file=sys.stderr)
                args.mode = 'spaced'
                args.multilengths = offset
            elif args.mode == 'gapped':
                args.multilengths = np.array(offset)
            else:
                args.multilengths = offset
        except ValueError:
            if args.untangle or args.hide_empty_contigs:
                print('Error: --untangle / --hide-empty-contigs require multilengths in %s'
                      % args.lens, file=sys.stderr)
                sys.exit(1)
            print('Warning: Multi-FASTA lengths not available in %s. Treating input FASTAs as a single sequence instead.' % args.lens, file=sys.stderr)
            args.mode = 'normal'
            seq_lengths = get_sequence_lengths(args.lens)
            args.multilengths = None
    else:
        seq_lengths = get_sequence_lengths(args.lens)
        args.multilengths = None

    args.spacer_bp = 0
    if args.mode == 'gapped' and args.multilengths is not None:
        args.spacer = args.spacer * args.multilengths.max(axis=0).max()
        args.spacer_bp = int(args.spacer)
    elif args.mode == 'spaced' and args.multilengths is not None:
        max_contig = max((max(ml) for ml in args.multilengths if len(ml)), default=0)
        args.spacer_bp = int(args.spacer * max_contig) if max_contig else 0
    
    genome_order = None
    if args.filelist:
        cur_order = get_seq_paths(args.lens)
        file_order = [l.split()[0] for l in open(args.filelist, 'r').read().splitlines()]
        try:
            genome_order = np.array([cur_order.index(l) for l in file_order])
        except ValueError:
            print(f'Error: sequence in filelist not found in lengths file. Ensure the filelist paths match the lengths file paths exactly.', file=sys.stderr)
            sys.exit(1)
        seq_lengths = [seq_lengths[i] for i in genome_order]
        if args.mode == 'gapped':
            args.multilengths = args.multilengths[genome_order]
        elif args.multilengths is not None:
            args.multilengths = [args.multilengths[i] for i in genome_order]
    
    if args.labels:
        if args.labels.endswith('.lengths'):
            genome_names = [os.path.splitext(os.path.basename(l))[0] for l in get_seq_paths(args.labels)]
            if genome_order is not None:
                genome_names = [genome_names[i] for i in genome_order]
        else:
            genome_names = open(args.labels, 'r').read().splitlines()
    else:
        genome_names = None
        
    # Use new MUM class
    mums = MUMdata(args.mumfile, lenfilter=args.lenfilter, subsample=args.subsample, verbose=args.verbose)
    if args.verbose:
        print(f'Found {mums.num_mums} MUMs'
              + (' (contig coords)' if mums.has_contig_coords else ''), file=sys.stderr)

    if genome_order is not None:
        mums.starts = mums.starts[:, genome_order]
        mums.strands = mums.strands[:, genome_order]
        if mums.contig_ids is not None:
            mums.contig_ids = mums.contig_ids[:, genome_order]

    # ---- Contig layout / trailing-empty trim ------------------------------
    # Spaced + hide-empty (default) needs contig-local coords to trim, then
    # prepare_contig_plot applies the spacer. Also used for --untangle.
    contig_mode = (
        args.untangle
        or mums.has_contig_coords
        or args.hide_empty_contigs
        or args.mode == 'spaced'
    )
    if contig_mode:
        if args.multilengths is None:
            print('Error: contig coordinates / --untangle / --mode spaced require multilengths',
                  file=sys.stderr)
            sys.exit(1)
        ensure_contig_coords(mums, args.multilengths)
        # Drop partials before layout/blocks
        mums.filter_pmums()
        if len(mums) == 0:
            print('No strict MUMs found after filtering.', file=sys.stderr)
            return

        if args.hide_empty_contigs:
            mums, args.multilengths = trim_trailing_empty_contigs(
                mums, args.multilengths, verbose=args.verbose,
            )
            seq_lengths = [sum(ml) for ml in args.multilengths]

        layout = ContigLayout.identity(args.multilengths)
        if args.untangle:
            layout = untangle_layout(mums, args.multilengths, ref=0, verbose=args.verbose)

        # Blocks before flips; apply = contig-id remap + strand flip only
        if args.max_break is None:
            bp_per_inch = max(seq_lengths) / (args.dpi * args.size[0])
            args.max_break = int(min(bp_per_inch, 100000))
        # Recompute spacer after trim / filelist reorder of multilengths
        if args.mode == 'spaced':
            max_contig = max((max(ml) for ml in args.multilengths if len(ml)), default=0)
            args.spacer_bp = int(args.spacer * max_contig) if max_contig else 0
        mums, blocks, args.multilengths = prepare_contig_plot(
            mums, args.multilengths, layout,
            max_break=args.max_break, verbose=args.verbose,
            spacer=args.spacer_bp if args.mode == 'spaced' else 0,
        )
        mums.blocks = blocks
        if args.mode == 'spaced':
            seq_lengths = [
                gapped_genome_length(ml, args.spacer_bp) for ml in args.multilengths
            ]
        else:
            seq_lengths = [sum(ml) for ml in args.multilengths]
        if args.mode == 'normal':
            args.mode = 'spaced'
        if args.verbose:
            print(f'contig layout applied (ids+strands only); '
                  f'{mums.num_mums} MUMs, {len(blocks)} blocks'
                  + (f'; spacer={args.spacer_bp}bp' if args.mode == 'spaced' else ''),
                  file=sys.stderr)

    max_length = max(seq_lengths)
    
    centering = [0] * len(seq_lengths)
    if args.center:
        centering = [(max_length - g) / 2 for g in seq_lengths]
            
    if args.no_coll_block:
        if args.mode == 'gapped':
            offset_mums(args, mums)
        elif args.mode == 'spaced' and not contig_mode:
            insert_spaced_spacers(mums, args.multilengths, args.spacer_bp)
        polygons, colors = get_mum_polygons(mums, centering, color=args.mum_color, inv_color=args.inv_color)
    else:
        ### filter out pmums for collinear blocks
        if mums.blocks is None:
            if not contig_mode:
                mums.filter_pmums()
            if len(mums) == 0:
                print('No strict MUMs found after filtering. Try turning off collinear blocking with --no-coll-block', file=sys.stderr)
                return
            if args.max_break is None:
                bp_per_inch = max_length / (args.dpi * args.size[0])
                args.max_break = int(min(bp_per_inch, 100000))
            collinear_blocks = find_coll_blocks(mums, max_break=args.max_break, verbose=args.verbose)
            if args.verbose:
                print(f'found {len(collinear_blocks)} collinear blocks', file=sys.stderr)
        else:
            if args.verbose:
                print(f'Using pre-computed collinear blocks: {len(mums.blocks)} blocks', file=sys.stderr)
            collinear_blocks = mums.blocks
        if args.mode == 'gapped':
            collinear_blocks = offset_mums(args, mums, spacer=args.spacer, blocks=collinear_blocks)
        elif args.mode == 'spaced' and not contig_mode:
            collinear_blocks = insert_spaced_spacers(
                mums, args.multilengths, args.spacer_bp, blocks=collinear_blocks,
            )
            
        polygons, colors = get_block_polygons(collinear_blocks, mums, centering, color=args.mum_color, inv_color=args.inv_color)

    # Spaced: xlim follows dense MUM mass, not a few outlier midpoints on junk tails
    if args.mode == 'spaced' and args.region is None:
        dense = mum_dense_xmax(mums, quantile=0.99)
        poly_x = _polygons_xmax(polygons) if polygons else 0.0
        # Prefer dense quantile; never exceed polygon max
        cand = dense if dense > 0 else poly_x
        if poly_x > 0:
            cand = min(cand, poly_x) if cand > 0 else poly_x
        if cand > 0:
            args.plot_xmax = cand
            if args.verbose:
                print(f'spaced xlim from dense MUM content: {cand:.0f} bp '
                      f'(backbone max {max(seq_lengths)})', file=sys.stderr)

    if args.verbose:
        print('Rendering plot...', file=sys.stderr)
    plot(args, seq_lengths, polygons, colors, centering, genomes=genome_names, filename=args.filename, dpi=args.dpi, size=args.size)
    if args.verbose:
        print('Done.', file=sys.stderr)

if __name__ == "__main__":
    args = parse_arguments()
    main(args)
