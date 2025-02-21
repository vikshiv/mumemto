#!/usr/bin/env python3

from matplotlib import pyplot as plt
from matplotlib import cm
from matplotlib.collections import PolyCollection
import os, sys
import numpy as np
import argparse
try:
    from utils import MUMdata, find_coll_blocks, get_sequence_lengths, get_seq_paths
except ImportError:
    from mumemto.utils import MUMdata, find_coll_blocks, get_sequence_lengths, get_seq_paths

def parse_arguments(args=None):    
    parser = argparse.ArgumentParser(description="Plots a synteny plot of MUMs from mumemto")
    # parser.add_argument('--filelist', '-f', dest='filelist', help='path to filelist from mumemto', required=True)
    # parser.add_argument('--mums', '-m', dest='mumfile', help='path to *.mum file from mumemto', required=True)
    # parser.add_argument('--lengths','-l', dest='lens', help='lengths file, first column is seq length in order of filelist', required=True)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--input-prefix', '-i', dest='prefix', help='prefix for filelist, mums, and lengths files')
    group.add_argument('--mums', '-m', dest='mumfile', help='path to *.mum file from mumemto')
    
    parser.add_argument('--lengths','-l', dest='lens', help='lengths file, first column is seq length in order of filelist')
    
    parser.add_argument('--filelist', '-f', dest='filelist', help='if the filelist is provided, then FASTA filenames are used as labels')
    parser.add_argument('--len-filter','-L', dest='lenfilter', help='only plot MUMs longer than threshold', default=0, type=int)
    parser.add_argument('--subsample','-s', dest='subsample', help='subsample every Nth mum', default=1, type=int)
    parser.add_argument('--center','-c', dest='center', action='store_true', help='center plot', default=False)
    parser.add_argument('--inversion-color','-ic', dest='inv_color', help='color for inversions', default='green')
    parser.add_argument('--mum-color','-mc', dest='mum_color', help='color as hex (default: #00A2FF)', default='#00A2FF', type=str)
    parser.add_argument('--alpha','-a', dest='alpha', help='opacity of mums [0-1] (default: 0.8 for blocks, 0.1 for MUMs)', type=float)
    parser.add_argument('--linewidth','-lw', dest='linewidth', help='linewidth of mums (default: 0 for blocks, 0.05 for MUMs)', type=float)
    parser.add_argument('--fout','-o', dest='filename', help='plot fname (default: input_prefix)')
    parser.add_argument('--dims', dest='size', help='fig dimensions (inches) (default: 6.4, 4.8)', default=(6.4,4.8), type=float, nargs=2)
    parser.add_argument('--dpi','-d', dest='dpi', help='dpi', default=500, type=int)
    parser.add_argument('--verbose','-v', dest='verbose', help='verbose mode', action='store_true', default=False)
    parser.add_argument('--no-coll-block','-b', dest='no_coll_block', help='plot only MUMs, not collinear blocks (slower) (default: false)', action='store_true', default=False)
    parser.add_argument('--max-gap-len','-g', dest='max_break', help='maximum break between collinear mums within a collinear block (default: <1px)', default=None, type=int)
    
    # handle multi-fasta inputs
    parser.add_argument('--mode', dest='mode', choices=['normal', 'delineated', 'gapped'],
                       help='Plotting mode for multi-fasta inputs (default: normal)', default='normal')
    parser.add_argument('--spacer', dest='spacer', help='spacer between contigs, expressed as proportion (0 - 1) of largest contig length (default: 0.1)', default=0.1, type=float)

    if args is None:
        args = parser.parse_args()
    else:
        args = parser.parse_args(args)
    
    if args.mumfile:
        args.prefix = os.path.splitext(args.mumfile)[0]
    elif args.prefix:
        if args.prefix.endswith('.mums'):
            args.prefix = args.prefix[:-5]
        args.mumfile = args.prefix + '.mums'
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
        
    # if (args.mode == 'gapped' or args.mode == 'delineated') and args.multilengths is None:
    #     parser.error('multi-lengths file is required for gapped or delineated mode')
    return args

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
        left, right = mums[l].starts, mums[r].starts + mums[r].length
        for idx, strand in enumerate(strands):
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

def plot(args, genome_lengths, polygons, colors, centering, dpi=500, size=None, genomes=None, filename=None):
    fig, ax = plt.subplots()
    max_length = max(genome_lengths)
    # Plot genome lines based on mode
    if args.mode == 'normal':
        # Just plot simple genome lines
        for idx, g in enumerate(genome_lengths):
            ax.plot([centering[idx] + 0, centering[idx] + g], [idx, idx], 
                    alpha=0.2, linewidth=0.75)
            
    elif args.mode == 'delineated':
        # Plot lines with delineators for multifasta
        offsets = args.multilengths.cumsum(axis=1)
        for idx in range(offsets.shape[0]):
            last_offset = 0
            for i, offset in enumerate([0] + offsets[idx][:-1].tolist()):
                ax.plot([centering[idx] + offset, centering[idx] + offset], 
                        [idx - 0.25, idx + 0.25], alpha=0.5, linewidth=0.25, color=cm.tab20((i+1) % 20))
                ax.plot([centering[idx] + last_offset, centering[idx] + offset], 
                        [idx, idx], alpha=0.2, linewidth=0.75, color=cm.tab20(i % 20))
                last_offset = offset
                
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
    
    ax.yaxis.set_ticks(list(range(len(genome_lengths))))
    ax.tick_params(axis='y', which='both',length=0)
    if genomes:
        ax.set_yticklabels(genomes)
    else:
        ax.yaxis.set_ticklabels([])
    if args.mode == 'gapped':
        ax.set_xlabel('chromosome')
    else:
        ax.set_xlabel('genomic position')
    ax.set_ylabel('sequences')
    ax.set_ylim(0, len(genome_lengths)-1)
    if args.mode == 'gapped':
        ax.set_xlim(0, args.multilengths.max(axis=0).sum() + args.spacer * (args.multilengths.shape[1] - 1))
    else:
        ax.set_xlim(0, max_length)
    ax.invert_yaxis()
    fig.set_tight_layout(True)
    # ax.axis('off')
    
    # ax.get_yaxis().set_visible(False)
    # ax.spines['top'].set_visible(False)
    # ax.spines['right'].set_visible(False)
    # ax.spines['left'].set_visible(False)
    # ax.tick_params(axis='y', which='both', left=False, right=False, labelleft=False)
    # ax.grid(False)

    if size:
        fig.set_size_inches(*size)
    if filename:
        if os.path.splitext(filename)[1] not in ['.png', '.pdf', '.svg']:
            filename = filename + '.png'
        if os.path.dirname(filename):
            path = filename
        else:
            path = os.path.join(os.path.dirname(args.mumfile), filename)
        fig.savefig(path, dpi=dpi)
    return ax

def offset_mums(args, mums, spacer=100000):
    offset = args.multilengths
    NUM_SEQS = len(offset)
    ### offset the mums by the contig locations
    offsets = np.cumsum(offset, axis=1)
    ### label each mum with the contig it belongs to
    breaks = np.hstack((np.array([[0]*NUM_SEQS]).transpose(), offsets))
    contig_idx = np.array([np.searchsorted(breaks[idx], mums.starts[:,idx]) - 1 for idx in range(NUM_SEQS)]).transpose()
    ### get the relative offset of each mum to the start of its contig
    left_start = np.hstack((np.zeros((offsets.shape[0],1), dtype=int), offsets[:,:-1]))
    rel_offsets = mums.starts - left_start[np.arange(NUM_SEQS), contig_idx]
    partial_mask = mums.starts != -1
    new_starts = np.array([0] + (offset.max(axis=0) + spacer).cumsum().tolist()[:-1])[contig_idx] + rel_offsets
    mums.starts[partial_mask] = new_starts[partial_mask]

def main(args):
    if args.mode != 'normal':
        try:
            offset = get_sequence_lengths(args.lens, multilengths=True)
            if args.mode == 'gapped' and len(set([len(o) for o in offset])) > 1:
                print('Warning: gapped mode requires the same number of sequences per input FASTA file. Using delineated mode instead.', file=sys.stderr)
                args.mode = 'delineated'
            seq_lengths = offset.sum(axis=1).tolist()
            args.multilengths = offset
        except ValueError:
            print('Warning: Multi-FASTA lengths not available in %s. Treating input FASTAs as a single sequence instead.' % args.lens, file=sys.stderr)
            args.mode = 'normal'
            seq_lengths = get_sequence_lengths(args.lens)
    else:
        seq_lengths = get_sequence_lengths(args.lens)
    
    if args.mode == 'gapped':
        args.spacer = args.spacer * args.multilengths.max(axis=0).max()
    if args.filelist:
        genome_names = get_seq_paths(args.filelist)
        genome_names = [os.path.splitext(os.path.basename(l))[0] for l in genome_names]
    else:
        genome_names = None
    max_length = max(seq_lengths)
    
    # Use new MUM class
    mums = MUMdata(args.mumfile, lenfilter=args.lenfilter, subsample=args.subsample, verbose=args.verbose)
    if args.verbose:
        print(f'Found {mums.num_mums} MUMs', file=sys.stderr)

    centering = [0] * len(seq_lengths)
    if args.center:
        centering = [(max_length - g) / 2 for g in seq_lengths]
            
    if args.no_coll_block:
        if args.mode == 'gapped':
            offset_mums(args, mums)
        polygons, colors = get_mum_polygons(mums, centering, color=args.mum_color, inv_color=args.inv_color)
    else:
        ### filter out pmums for collinear blocks
        if mums.blocks is None:
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
            offset_mums(args, mums, spacer=args.spacer)
            
        polygons, colors = get_block_polygons(collinear_blocks, mums, centering, color=args.mum_color, inv_color=args.inv_color)
    
    if args.verbose:
        print('Rendering plot...', file=sys.stderr)
    plot(args, seq_lengths, polygons, colors, centering, genomes=genome_names, filename=args.filename, dpi=args.dpi, size=args.size)
    if args.verbose:
        print('Done.', file=sys.stderr)

if __name__ == "__main__":
    args = parse_arguments()
    main(args)
