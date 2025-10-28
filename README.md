![GitHub release (latest by date)](https://img.shields.io/github/v/release/vikshiv/mumemto) ![GitHub](https://img.shields.io/github/license/vikshiv/mumemto?color=green) [![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/mumemto/README.html)
# **mumemto**: finding multi-MUMs and MEMs in pangenomes

<img src="img/polaroid_tattoo.png" alt="logo" width="310" align="left"/>

Mumemto is a tool for analyzing pangenome sequence collections. It identifies **maximal unique/exact matches (multi-MUMs and multi-MEMs)** present across a collection of sequences. Mumemto can **visualize** pangenome synteny, **identify misassemblies**, and provide a unifiying structure to a pangenome.

This method is uses the prefix-free parse (PFP) algorithm for suffix array construction on large, repetitive collections of text. The main workflow of `mumemto` is to compute the PFP over a collection of sequences, and identify multi-MUMs while computing the SA/LCP/BWT of the input collection. Note that this works best with highly repetitive texts (such as a collection of closely related genomes, likely intra-species such as a pangenome).

---
If you use Mumemto, please cite:

    Shivakumar, V. S., & Langmead, B. (2025). Mumemto: efficient maximal matching across pangenomes. Genome Biology, 26(1), 169.

If you use the partition-merging approach (for parallelization, etc.), please cite:

    Shivakumar, V. S., & Langmead, B. (2025). Partitioned Multi-MUM finding for scalable pangenomics. bioRxiv, 2025-05.

Paper available at [https://genomebiology.biomedcentral.com/articles/10.1186/s13059-025-03644-0](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-025-03644-0).

Merging algorithm described in [https://www.biorxiv.org/lookup/doi/10.1101/2025.05.20.654611](https://www.biorxiv.org/lookup/doi/10.1101/2025.05.20.654611).


### For detailed usage, see the [Mumemto wiki](https://github.com/vikshiv/mumemto/wiki)

## Installation

### Conda installation (recommended)
Mumemto is available on `bioconda`. Note conda installation requires python 3.9+. We recommend using a new environment:
```sh
### conda ###
conda create -n mumemto_env python #3.10 or higher
conda activate mumemto_env

conda install -c conda-forge bioconda::mumemto
```

### pip installation
Mumemto can be installed using pip locally. We recommend using a new conda or virtual environment for this.
```sh
git clone https://github.com/vikshiv/mumemto
cd mumemto; pip install .
```

### Docker/Singularity
Mumemto is available on `docker` and `singularity`. 
> [!TIP]
>  You may need to bind a local directory to access files in the container, which may cause issues when globbing input files. Input filelist + docker/singularity bind mount is recommended.
```sh
### if using docker ###
docker pull vshiv123/mumemto:latest
docker run vshiv123/mumemto:latest -h

### if using singularity ###
singularity pull mumemto.sif docker://vshiv123/mumemto:latest
./mumemto.sif -h
 # or any subcommand, e.g. ./mumemto.sif viz
```

### Compile from scratch
To build from scratch, download the source code and use cmake/make. After running the make command below,
the `mumemto` executable will be found in the `build/` folder. The following are dependencies: cmake, g++, gcc

```sh
git clone https://github.com/vikshiv/mumemto
cd mumemto

mkdir build 
cd build && cmake ..
make install
```

Note: For the python scripts, you may need to install dependencies separately. The following are dependencies: matplotlib, numpy, tqdm, plotly (for interactive plots), and numba (for the coverage script).

## Quick start
To visualize the synteny across the FASTA files in a directory `assemblies/` (each sequence is a separate fasta file):
```sh
mumemto assemblies/*.fa -o pangenome
mumemto viz -i pangenome
```

## Getting started

### Find multi-MUMs and MEMs
By default, `mumemto` computes multi-MUMs across a collection, without additional parameters. 
```sh
mumemto -o <output_prefix> [input_fasta [...]]
```

### Command line options
Use the `-h` flag to list additional options and usage: `mumemto -h`.

Mumemto options enable the computation of various different classes of exact matches:
<p align="center">
<img src="img/viz_def.png" alt="visual_guide" width="600" align="center"/>
</p>

The multi-MUM properties can be loosened to find different types of matches with three main flags: 
- `-k` determines the minimum number of sequences a match must occur in (e.g. for finding MUMs across smaller subsets)
- `-f` controls the maximum number of occurences in _each_ sequence (e.g. finding duplication regions)
- `-F` controls the total number of occurences in the collection (e.g. filtering out matches that occur frequently due to low complexity)

>[!TIP] 
> `-k` is flexible in input format. The user can specify a positive integer, indicating the minimum number of sequences a match should appear in. Passing a negative integer indicates a subset size relative to N, the number of sequences in the collection (i.e. N - k). For instance, to specify a match must appear in at least all sequences _except_ one, we could pass `-k -1`. Similarly, passing negative values to `-F` specifies limits relative to N. Note: when setting `-F` and `-f` together, the max total limit will be the smaller of `F` and `N * f`.

Here are some example use cases:

```sh
	 # Find all strict multi-MUMs across a collection
     mumemto [OPTIONS] [input_fasta [...]] (equivalently -k 0 -f 1 -F 0)
	 # Find partial multi-MUMs in all sequences but one
     mumemto -k -1 [OPTIONS] [input_fasta [...]]
	 # Find multi-MEMs that appear at most 3 times in each sequence
     mumemto -f 3 [OPTIONS] [input_fasta [...]]
	 # Find all MEMs that appear at most 100 times within a collection
     mumemto -f 0 -k 2 -F 100 [OPTIONS] [input_fasta [...]]
```

### Multi-MUM merging (***new in v1.3***)

The output multi-MUMs from Mumemto can be merged between runs in v1.3. There are two methods to do this: anchor-based (`-Mn`) and string-based (`-M`). Anchor-based merging requires the *first* sequence in each partition to be the same. String-based merging does not require any overlap between partitions, however is generally slower.

Running Mumemto with `-M` or `-Mn` generates a threshold file, `*.thresh` and `*.thresh_rev` for string-based merging and `*.athresh` for anchor-based merging. To merge partitions, run:
```
mumemto merge p1.mums p2.mums <...> -o <out_prefix>.mums
```
The merge script automatically detects which type of merging is possible and creates an output using `out_prefix` which is identical to if Mumemto was run on the union of the input partitions.

> [!NOTE]
>  Merging is currently limited to strict multi-MUMs. However, partial multi-MUMs for local partitions can be found using string-based merging incrementally.

> [!NOTE]
>  In v1.3.4, the string-based threshold file format changed and is not backwards compatible. To convert string-merging threshold files from v1.3.3 or earlier, use the provided `mumemto/convert_thresh.py`.

> [!TIP] 
> Using either merge mode enables a dynamic updating of multi-MUMs. You can incrementally add assemblies as the pangenome grows and update the global set of multi-MUMs across the collection.

### I/O format
The `mumemto` command takes in a list of fasta files as positional arguments and then generates output files using the output prefix. Alternatively, you can provide a file-list, which specifies a list of fastas (one per line). Passing in fastas as positional arguments will auto-generate a filelist that defines the order of the sequences in the output. 

> [!TIP]
> The output `*.lengths` file can also serve as an input filelist to re-run expts.

**Example of file-list file:**
```sh
/path/to/ecoli_1.fna
/path/to/salmonella_1.fna
/path/to/bacillus_1.fna
/path/to/staph_2.fna
```

**Format of the output \*.mums file:**
```sh
[MUM length] [comma-delimited list of offsets within each sequence, in order of filelist] [comma-delimited strand indicators (one of +/-)]
```
If the maximum number of occurences _per_ sequence is set to 1 (indicating MUMs), a `*.mums` file is generated. This contains each MUM as a separate line, where the first value is the match length, and the second is 
a comma-delimited list of positions where the match begins in each sequence. An empty entry indicates that the MUM was not found in that sequence (only applicable with *-k* flag). The MUMs are sorted in the output file
lexicographically based on the match sequence.

**Format of the output \*.mems file:**
```sh
[MEM length] [comma-delimited list of offsets for each occurence] [comma-delimited list of sequence IDs, as defined in the filelist] [comma-delimited strand indicators (one of +/-)]
```
If more than one occurence is allowed per sequence, the output format is in `*.mems` format. This contains each MEM as a separate line with the following fields: (1) the match length, (2)
a comma-delimited list of offsets within a sequence, (3) the corresponding sequence ID for each offset given in (2). Similar to above, MEMs are sorted in the output file
lexicographically based on the match sequence.

**Bumbl format:**

Mumemto can also output a binary format (`*.bumbl`) for faster I/O with large MUM files. Most Mumemto commands accept either `*.mums` or `*.bumbl` files interchangeably. Mumemto can output a `*.bumbl` file using the ``-b` flag. To convert between formats:
```sh
mumemto convert -m input.mums -o output.bumbl  # text to binary
mumemto convert -b input.bumbl -o output.mums  # binary to text
```

> [!TIP]
> For large pangenomes, using `*.bumbl` files can significantly reduce file sizes and improve loading times for visualization and analysis.

## Visualization
<figure>
<img src="img/potato_syn.png" alt="potato_synteny"/>
<figcaption> <p align="center">Potato pangenome (assemblies from <a href='https://www.nature.com/articles/s41586-024-08476-9'>[Cheng <i>et al.</i>, 2025]</a>)</p></figcaption>
</figure>
Mumemto can visualize multi-MUMs in a synteny-like format, highlighting conservation and genomic structural diversity within a collection of sequences.

After running `mumemto` on a collection of FASTAs, you can generate a visualization using:
```sh
mumemto viz (-i PREFIX | -m MUMFILE)
```
Use `mumemto viz -h` to see options for customizability. As of now, only strict and partial multi-MUMs are supported (rare multi-MEM support coming soon), thus a `*.mums` output is required. 

An interactive plot (with plotly, still experimental) can be generated with `mumemto viz --interactive`.

## Getting Help

If you run into any issues or have any questions, please feel free to reach out to us either (1) through GitHub Issues or (2) reach out to me at vshivak1 [at] jhu.edu

## Acknowledgements

Portions of code from this repo were adapted from <a href="https://github.com/maxrossi91/pfp-thresholds">pfp-thresholds</a>, written by <a href="https://github.com/maxrossi91">Massimiliano Rossi</a> and <a href="https://github.com/oma219/docprofiles">cliffy</a>, written by <a href="https://github.com/oma219">Omar Ahmed</a>. 

## Citing Mumemto and reproducing results

Preprint: [https://doi.org/10.1101/2025.01.05.631388](https://doi.org/10.1101/2025.01.05.631388)

Scripts to reproduce the results in the preprint are found in this repo: https://github.com/vikshiv/mumemto-reproducibility-scripts
