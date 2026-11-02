# Contig-layout untangle eval

Prototype for `--untangle` on bag-of-contigs inputs: reorder/orient non-reference
genomes by median mapped position vs genome 0, plotted with `--mode spaced`.

## Yeast4 eval

Requires `mumento2` at:

`/hps/nobackup/rdf/research/vshiv/mumemto/kmer_mumemto/target/release/mumento2`

Assemblies in `../data/data/` (R64, ASM308665v1, UWOPS83-787, AGK).

```bash
cd test/untangle
./run_yeast4.sh
```

Outputs under `runs/yeast4/`:

- `default.png` — FASTA contig order (`--mode spaced`)
- `untangle.png` — `--untangle` (ref-fixed greedy)

Optional tighter/looser gaps: pass `--spacer 0.05` or `--spacer 0.2` (fraction of largest contig; default `0.1`).

Trailing contigs with no MUM hits are dropped by default in spaced mode (`--show-empty-contigs` to keep them).
