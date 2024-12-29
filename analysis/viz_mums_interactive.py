#!/usr/bin/env python3

import plotly.graph_objects as go
import os
import argparse
import numpy as np
import sys
from mumemto.utils import find_coll_blocks, MUMdata

def parse_arguments():    
    parser = argparse.ArgumentParser(description="Plots an interactive synteny plot of MUMs from mumemto")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--input-prefix', '-i', dest='prefix', help='prefix for filelist, mums, and lengths files')
    group.add_argument('--mums', '-m', dest='mumfile', help='path to *.mum file from mumemto')
    
    parser.add_argument('--lengths','-l', dest='lens', help='lengths file, first column is seq length in order of filelist')
    parser.add_argument('--filelist', '-f', dest='filelist', help='if the filelist is provided, then FASTA filenames are used as labels')
    parser.add_argument('--len-filter','-L', dest='lenfilter', help='only plot MUMs longer than threshold', default=0, type=int)
    parser.add_argument('--subsample','-s', dest='subsample', help='subsample every Nth mum', default=1, type=int)
    parser.add_argument('--center','-c', dest='center', action='store_true', help='center plot', default=False)
    parser.add_argument('--inversion-color','-ic', dest='inv_color', help='color for inversions', default='green')
    parser.add_argument('--mum-color','-mc', dest='mum_color', help='color for MUMs', default='rgba(0, 162, 255, 0.5)')
    parser.add_argument('--fout','-o', dest='filename', help='plot fname (default: input_prefix)')
    parser.add_argument('--dims', dest='size', help='fig dimensions (pixels) (default: 1000, 600)', default=(1000, 600), type=int, nargs=2)
    parser.add_argument('--verbose','-v', dest='verbose', help='verbose mode', action='store_true', default=False)
    parser.add_argument('--no-coll-block','-b', dest='no_coll_block', help='plot only MUMs, not collinear blocks (slower) (default: false)', action='store_true', default=False)
    parser.add_argument('--max-gap-len','-g', dest='max_break', help='maximum break between collinear mums within a collinear block (default: 1000)', default=1000, type=int)
    
    args = parser.parse_args()
    if args.mumfile:
        args.prefix = os.path.splitext(args.mumfile)[0]
        lens = args.prefix + '.lengths'
    else:
        args.mumfile = args.prefix + '.mums'
        lens = args.prefix + '.lengths'
        
    if args.lens is None:
        args.lens = lens
    if not args.filename:
        args.filename = args.prefix
    return args

def get_mum_shapes(mums, centering, color='rgba(0, 162, 255, 0.5)', inv_color='green'):
    shapes = []
    for i in range(mums.num_mums):
        mum = mums[i]
        inverted = not mum.strands[0]
        points = []
        for idx, (x, strand) in enumerate(zip(mum.starts, mum.strands)):
            if x == -1:
                continue
            points.append(((centering[idx] + x, idx), (centering[idx] + x + mum.length, idx)))
            if not inverted and not strand:
                inverted = True
                if len(points) > 2:
                    shapes.append(dict(
                        type="path",
                        path=make_polygon_path(points[:-1]),
                        fillcolor=color,
                        line_color=color,
                        opacity=0.5
                    ))
                shapes.append(dict(
                    type="path",
                    path=make_polygon_path(points[-2:]),
                    fillcolor=inv_color,
                    line_color=inv_color,
                    opacity=0.5
                ))
                points = [points[-1]]
            elif inverted and strand:
                inverted = False
                if len(points) > 2:
                    shapes.append(dict(
                        type="path",
                        path=make_polygon_path(points[:-1]),
                        fillcolor=color,
                        line_color=color,
                        opacity=0.5
                    ))
                shapes.append(dict(
                    type="path",
                    path=make_polygon_path(points[-2:]),
                    fillcolor=inv_color,
                    line_color=inv_color,
                    opacity=0.5
                ))
                points = [points[-1]]
        if len(points) >= 2:
            shapes.append(dict(
                type="path",
                path=make_polygon_path(points),
                fillcolor=color,
                line_color=color,
                opacity=0.5
            ))
    return shapes

def get_block_shapes(collinear_blocks, mums, centering, color='rgba(0, 162, 255, 0.5)', inv_color='green'):
    shapes = []
    for l, r in collinear_blocks:
        left_mum = mums[l]
        right_mum = mums[r]
        strands = left_mum.strands
        inverted = not strands[0]
        points = []
        left, right = left_mum.starts, right_mum.starts + right_mum.length
        for idx, strand in enumerate(strands):
            points.append(((centering[idx] + left[idx], idx), (centering[idx] + right[idx], idx)))
            if not inverted and not strand:
                inverted = True
                if len(points) > 2:
                    shapes.append(dict(
                        type="path",
                        path=make_polygon_path(points[:-1]),
                        fillcolor=color,
                        line_color=color,
                        opacity=0.5
                    ))
                shapes.append(dict(
                    type="path",
                    path=make_polygon_path(points[-2:]),
                    fillcolor=inv_color,
                    line_color=inv_color,
                    opacity=0.5
                ))
                points = [points[-1]]
            elif inverted and strand:
                inverted = False
                if len(points) > 2:
                    shapes.append(dict(
                        type="path",
                        path=make_polygon_path(points[:-1]),
                        fillcolor=color,
                        line_color=color,
                        opacity=0.5
                    ))
                shapes.append(dict(
                    type="path",
                    path=make_polygon_path(points[-2:]),
                    fillcolor=inv_color,
                    line_color=inv_color,
                    opacity=0.5
                ))
                points = [points[-1]]
        if len(points) >= 2:
            shapes.append(dict(
                type="path",
                path=make_polygon_path(points),
                fillcolor=color,
                line_color=color,
                opacity=0.5
            ))
    return shapes

def make_polygon_path(points):
    starts, ends = tuple(zip(*points))
    points = starts + ends[::-1]
    path = f"M {points[0][0]},{points[0][1]}"
    for x, y in points[1:]:
        path += f" L{x},{y}"
    path += " Z"
    return path

def plot(genome_lengths, shapes, centering, genomes=None, size=(1000, 600), filename=None):
    max_length = max(genome_lengths)
    
    # Create base lines for genomes
    base_lines = []
    for idx, length in enumerate(genome_lengths):
        base_lines.append(dict(
            type="line",
            x0=centering[idx],
            x1=centering[idx] + length,
            y0=idx,
            y1=idx,
            line=dict(color="gray", width=1)
        ))
    
    # Create the figure
    fig = go.Figure()
    
    # Add shapes
    fig.update_layout(
        shapes=base_lines + shapes,
        showlegend=False,
        xaxis=dict(
            title="genomic position",
            showgrid=False,
            zeroline=False
        ),
        yaxis=dict(
            title="sequences",
            showgrid=False,
            zeroline=False,
            tickmode='array',
            ticktext=genomes if genomes else [f"seq_{i}" for i in range(len(genome_lengths))],
            tickvals=list(range(len(genome_lengths))),
            autorange="reversed"
        ),
        width=size[0],
        height=size[1],
        plot_bgcolor='white'
    )
    
    # Set axis ranges
    fig.update_xaxes(range=[0, max_length])
    fig.update_yaxes(range=[-0.5, len(genome_lengths) - 0.5])
    
    if filename:
        filename = filename + ('' if filename.endswith('.html') else '.html')
        if os.path.dirname(filename):
            path = filename
        else:
            path = os.path.join(os.path.dirname(args.mumfile), filename)
        fig.write_html(path)
    
    return fig

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
        if args.verbose:
            print('Centering plot...', file=sys.stderr)
        centering = [(max_length - g) / 2 for g in seq_lengths]
    
    if args.no_coll_block:
        if args.verbose:
            print('Building synteny plot...', file=sys.stderr)
        shapes = get_mum_shapes(mums, centering, color=args.mum_color, inv_color=args.inv_color)
    else:
        mums.filter_pmums()
        if args.verbose:
            print(f'Finding collinear blocks (max gap = {args.max_break} bp)...', file=sys.stderr, end=' ')
        _, collinear_blocks, _ = find_coll_blocks(mums, max_break=args.max_break, verbose=args.verbose)
        if args.verbose:
            print(f'found {len(collinear_blocks)} collinear blocks', file=sys.stderr)
            print('Building synteny plot...', file=sys.stderr)
        shapes = get_block_shapes(collinear_blocks, mums, centering, color=args.mum_color, inv_color=args.inv_color)
    
    if args.verbose:
        print('Rendering plot...', file=sys.stderr)
    plot(seq_lengths, shapes, centering, genomes=genome_names, size=args.size, filename=args.filename)
    
    if args.verbose:
        print('Done.', file=sys.stderr)

if __name__ == "__main__":
    args = parse_arguments()
    main(args)