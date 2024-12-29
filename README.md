# **mumemto**: finding multi-MUMs and MEMs in pangenomes

<img src="img/polaroid_tattoo.png" alt="logo" width="292" align="left"/>

Mumemto identifies **maximal unique matches (multi-MUMs)** present across a collection of sequences. Multi-MUMs are defined as maximally matching substrings present in each sequence in a collection *exactly once*. Additionally, this tool can identify **multi-MEMs**, maximal exact matches present across sequences, without the uniqueness property. This method is uses the prefix-free parse (PFP) algorithm for suffix array construction on large, repetitive collections of text.

This tool uses PFP to efficiently identify multi-MUM/MEMs. Note that this works best with highly repetitive texts (such as a collection of closely related genomes, likely intra-species such as a pangenome). Mumemto does work for more divergent sequences (inter-species, etc.), however it is less memory efficient and may not scale beyond a few genomes (~100Gbp total).

The base code from this repo was adapted from <a href="https://github.com/maxrossi91/pfp-thresholds">pfp-thresholds</a> repository written by <a href="https://github.com/maxrossi91">Massimiliano Rossi</a> and <a href="https://github.com/oma219/docprofiles">docprofiles</a> repository written by <a href="https://github.com/oma219">Omar Ahmed</a>. 

## Installation

### Conda and pip installation (recommended)
Mumemto is available on `bioconda`, or can be installed with pip:
```
### conda ###
conda install -c bioconda mumemto

### pip ###
git clone https://github.com/vshiv18/mumemto -b dev
pip install .
```

### Docker/Singularity
Mumemto is available on `docker` and `singularity`:
```
### if using docker ###
docker pull vshiv123/mumemto:latest
docker run vshiv123/mumemto:latest mumemto -h

### if using singularity ###
singularity pull docker://vshiv123/mumemto:latest
./mumemto_latest.sif mumemto -h
```

### Compile from scratch
To build from scratch, download the source code and use cmake/make. After running the make command below,
the `mumemto` executable will be found in the `build/` folder. The following are dependencies: cmake, g++, gcc

```sh
git clone https://github.com/vshiv18/mumemto
cd mumemto

mkdir build 
cd build && cmake ..
make install
```

When compiling from scratch, the downstream python scripts will not be in the appropriate $PYTHONPATH. For these scripts, run the relevant python script directly from the `mumemto/` directory (you may need to install the python dependencies separately). 

## Getting started

The basic workflow with `mumemto` is to compute the PFP over a collection of sequences, and identify multi-MUMs while computing the SA/LCP/BWT of the input collection. 

### Find multi-MUMs and MEMs
By default, `mumemto` computes multi-MUMs across a collection, without additional parameters. 
```sh
mumemto -o <output_prefix> [input_fasta [...]]
```

The `mumemto` command takes in a list of fasta files as positional arguments and then generates output files using the output prefix. Alternatively, you can provide a file-list, which specifies a list of fastas and which document/class each file belongs in. Passing in fastas as positional arguments will auto-generate a filelist that defines the order of the sequences in the output.

**Example of file-list file:**
```sh
/path/to/ecoli_1.fna 1
/path/to/salmonella_1.fna 2
/path/to/bacillus_1.fna 3
/path/to/staph_2.fna 4
```

Use the `-h` flag to list additional options and usage: `mumemto -h`.

Mumemto mode options enable the computation of various different classes of exact matches:
<p align="center">
<img src="img/viz_def.png" alt="visual_guide" width="600" align="center"/>
</p>

The multi-MUM properties can be loosened to find different types of matches with three main flags: 
- `-k` determines the minimum number of sequences a match must occur in (e.g. for finding MUMs across smaller subsets)
- `-f` controls the maximum number of occurences in _each_ sequence (e.g. finding duplication regions)
- `-F` controls the total number of occurences in the collection (e.g. filtering out matches that occur frequently due to low complexity)

`-k` is flexible in input format. The user can specify a positive integer, indicating the minimum number of sequences a match should appear in. Passing a negative integer indicates a subset size relative to N, the number of sequences in the collection (i.e. N - k). For instance, to specify a match must appear in at least all sequences _except_ one, we could pass `-k -1`. Similarly, passing negative values to `-F` specifies limits relative to N. Note: when setting `-F` and `-f` together, the max total limit will be the smaller of `F` and `N * f`.

Here are some example use cases:

```
	 # Find all strict multi-MUMs across a collection
     mumemto [OPTIONS] [input_fasta [...]] (equivalently -k 0 -f 1 -F 0)
	 # Find partial multi-MUMs in all sequences but one
     mumemto -k -1 [OPTIONS] [input_fasta [...]]
	 # Find multi-MEMs that appear at most 3 times in each sequence
     mumemto -f 3 [OPTIONS] [input_fasta [...]]
	 # Find all MEMs that appear at most 100 times within a collection
     mumemto -f 0 -k 2 -F 100 [OPTIONS] [input_fasta [...]]
```

**Format of the \*.mums file:**
```sh
[MUM length] [comma-delimited list of offsets within each sequence, in order of filelist] [comma-delimited strand indicators (one of +/-)]
```
If the maximum number of occurences _per_ sequence is set to 1 (indiciating MUMs), a `*.mums` file is generated. This contains each MUM as a separate line, where the first value is the match length, and the second is 
a comma-delimited list of positions where the match begins in each sequence. An empty entry indicates that the MUM was not found in that sequence (only applicable with *-k* flag). The MUMs are sorted in the output file
lexicographically based on the match sequence.

**Format of the \*.mems file:**
```sh
[MEM length] [comma-delimited list of offsets for each occurence] [comma-delimited list of sequence IDs, as defined in the filelist] [comma-delimited strand indicators (one of +/-)]
```
If more than one occurence is allowed per sequence, the output format is in `*.mems` format. This contains each MEM as a separate line with the following fields: (1) the match length, (2)
a comma-delimited list of offsets within a sequence, (3) the corresponding sequence ID for each offset given in (2). Similar to above, MEMs are sorted in the output file
lexicographically based on the match sequence.


## Visualization
<figure>
<img src="img/potato_syn_small.png" alt="potato_synteny"/>
<figcaption> <p align="center">Potato pangenome (assemblies from <a href='https://www.nature.com/articles/s41586-022-04822-x'>[Tang <i>et al.</i>, 2022]</a>)</p></figcaption>
</figure>
Mumemto can visualize multi-MUMs in a synteny-like format, highlighting conservation and genomic structural diversity within a collection of sequences.

After running `mumemto` on a collection of FASTAs, you can generate a visualization using:
```sh
mumemto viz (-i PREFIX | -m MUMFILE)
```
Use `mumemto viz -h` to see options for customizability. As of now, only strict and partial multi-MUMs are supported (rare multi-MEM support coming soon), thus a `*.mums` output is required. An interactive plot (with plotly) can be generated with `mumemto viz --interactive`.


