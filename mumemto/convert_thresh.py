import numpy as np
import argparse
import os
import sys

from mumemto.utils import MUMdata


def parse_arguments():
    parser = argparse.ArgumentParser(description='Sort thresholds from mumemto <1.3.4')
    parser.add_argument('input_file', help='Path to MUMs file, prefix, or bumbl file')
    parser.add_argument('--output', '-o', help='Path to output converted threshold file', default=None)
    parser.add_argument('--verbose', '-v', action='store_true', help='Print verbose output')
    
    args = parser.parse_args()

    # Determine input file type and set threshold file path
    if args.input_file.endswith('.mums'):
        args.prefix = args.input_file[:-5]
    elif args.input_file.endswith('.bumbl'):
        args.prefix = args.input_file[:-6]
    else:
        args.prefix = args.input_file
        args.input_file += '.mums'
        
    args.thresh_file = args.prefix + '.thresh'
    args.thresh_rev_file = args.prefix + '.thresh_rev'
        
    # Set default output if not provided
    if args.output is None:
        args.output = args.prefix + '_converted'
    
    if not os.path.exists(args.input_file):
        print(f"Error: Input file {args.input_file} does not exist.", file=sys.stderr)
        sys.exit(1)
    # Check if threshold files exist
    if not os.path.exists(args.thresh_file):
        print(f"Error: Threshold file {args.thresh_file} does not exist.", file=sys.stderr)
        sys.exit(1)
    
    if not os.path.exists(args.thresh_rev_file):
        print(f"Error: Reverse threshold file {args.thresh_rev_file} does not exist.", file=sys.stderr)
        sys.exit(1)
    
    return args


def convert_threshold(args):
    """
    Convert thresholds to version >=1.3.4 order.
    """
    if args.verbose:
        print(f"Converting threshold file: {args.thresh_file}", file=sys.stderr)
        print(f"Output file: {args.output}", file=sys.stderr)
    
    with open(args.thresh_file, 'rb') as f:
        thresholds = np.fromfile(f, dtype=np.uint16)
    with open(args.thresh_rev_file, 'rb') as f:
        rev_thresholds = np.fromfile(f, dtype=np.uint16)
    
    new_thresholds = []
    new_rev_thresholds = []

    mums = MUMdata(args.input_file, sort=False)
    lengths = mums.lengths
    starts = [0] + np.cumsum(lengths + 1).tolist()
    order = np.argsort(mums[:, 0].starts)
    for o in order:
        new_thresholds.append(thresholds[starts[o] : starts[o] + lengths[o] + 1])
        new_rev_thresholds.append(rev_thresholds[starts[o] : starts[o] + lengths[o] + 1])
        
    concatenated_thresholds = np.concatenate(new_thresholds)
    concatenated_rev_thresholds = np.concatenate(new_rev_thresholds)
    assert concatenated_thresholds.size == thresholds.size
    assert concatenated_rev_thresholds.size == rev_thresholds.size
    
    with open(args.output + '.thresh', 'wb') as f:
        f.write(concatenated_thresholds.astype(np.uint16).tobytes())
    with open(args.output + '.thresh_rev', 'wb') as f:
        f.write(concatenated_rev_thresholds.astype(np.uint16).tobytes())

    mums[order].write_mums(args.output + '.mums')

def main(args):
    """
    Main function to run the threshold conversion.
    """
    if args.verbose:
        print("Starting threshold conversion...", file=sys.stderr)
    
    convert_threshold(args)
    
    if args.verbose:
        print("Threshold conversion completed.", file=sys.stderr)


if __name__ == "__main__":
    args = parse_arguments()
    main(args)
