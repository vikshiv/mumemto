from matplotlib import pyplot as plt
from matplotlib.collections import PolyCollection
import os, sys
import argparse
try:
    from utils import MUMdata, find_coll_blocks
except ImportError:
    from mumemto.utils import MUMdata, find_coll_blocks

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
    parser.add_argument('--alpha','-a', dest='alpha', help='opacity of mums [0-1] (default: 0.5)', default=0.5, type=float)
    parser.add_argument('--fout','-o', dest='filename', help='plot fname (default: input_prefix)')
    parser.add_argument('--dims', dest='size', help='fig dimensions (inches) (default: 6.4, 4.8)', default=(6.4,4.8), type=float, nargs=2)
    parser.add_argument('--dpi','-d', dest='dpi', help='dpi', default=500, type=int)
    parser.add_argument('--verbose','-v', dest='verbose', help='verbose mode', action='store_true', default=False)
    parser.add_argument('--no-coll-block','-b', dest='no_coll_block', help='plot only MUMs, not collinear blocks (slower) (default: false)', action='store_true', default=False)
    parser.add_argument('--max-gap-len','-g', dest='max_break', help='maximum break between collinear mums within a collinear block (default: <1px)', default=None, type=int)
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

    if not args.filename:
        args.filename = args.prefix
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

def plot(args, genome_lengths, polygons, colors, centering, alpha=0.5, dpi=500, size=None, genomes=None, filename=None, verbose=False):
    fig, ax = plt.subplots()
    max_length = max(genome_lengths)
    for idx, g in enumerate(genome_lengths):
        ax.plot([centering[idx] + 0,centering[idx] + g], [idx, idx], alpha=0.2, linewidth=0.75)
        
    ax.add_collection(PolyCollection(polygons, linewidths=0, alpha=alpha, edgecolors=colors, facecolors=colors))
    
    ax.yaxis.set_ticks(list(range(len(genome_lengths))))
    ax.tick_params(axis='y', which='both',length=0)
    if genomes:
        ax.set_yticklabels(genomes)
    else:
        ax.yaxis.set_ticklabels([])
    ax.set_xlabel('genomic position')
    ax.set_ylabel('sequences')
    ax.set_ylim(0, len(genome_lengths)-1)
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
        filename = filename + ('' if filename.endswith('.png') else '.png')
        if os.path.dirname(filename):
            path = filename
        else:
            path = os.path.join(os.path.dirname(args.mumfile), filename)
        fig.savefig(path, dpi=dpi)
    return ax

def main(args):
    seq_lengths = [int(l.split()[1]) for l in open(args.lens, 'r').read().splitlines()]
    if args.filelist:
        genome_names = [os.path.splitext(os.path.basename(l.split()[0]))[0] for l in open(args.filelist, 'r').read().splitlines()]
    else:
        genome_names = None
    max_length = max(seq_lengths)
    
    # Use new MUM class
    mums = MUMdata(args.mumfile, seq_lengths=seq_lengths, lenfilter=args.lenfilter, subsample=args.subsample, verbose=args.verbose)
    if args.verbose:
        print(f'Found {mums.num_mums} MUMs', file=sys.stderr)

    centering = [0] * len(seq_lengths)
    if args.center:
        centering = [(max_length - g) / 2 for g in seq_lengths]

    if args.no_coll_block:
        polygons, colors = get_mum_polygons(mums, centering, color=args.mum_color, inv_color=args.inv_color)
    else:
        ### filter out pmums for collinear blocks
        mums.filter_pmums()
        if args.max_break is None:
            bp_per_inch = max_length / (args.dpi * args.size[0])
            args.max_break = min(bp_per_inch, 100000)
        if args.verbose:
            print(f'Finding collinear blocks (max gap = {args.max_break} bp)...', file=sys.stderr, end=' ')
        _, collinear_blocks, _ = find_coll_blocks(mums, max_break=args.max_break, verbose=args.verbose)

        if args.verbose:
            print(f'found {len(collinear_blocks)} collinear blocks', file=sys.stderr)
        polygons, colors = get_block_polygons(collinear_blocks, mums, centering, color=args.mum_color, inv_color=args.inv_color)
    
    if args.verbose:
        print('Rendering plot...', file=sys.stderr)
    plot(args, seq_lengths, polygons, colors, centering, genomes=genome_names, alpha=args.alpha, filename=args.filename, dpi=args.dpi, size=args.size, verbose=args.verbose)
    if args.verbose:
        print('Done.', file=sys.stderr)

if __name__ == "__main__":
    args = parse_arguments()
    main(args)
