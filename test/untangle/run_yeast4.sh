#!/usr/bin/env bash
# Generate contig-coord MUMs for the 4 yeast assemblies and plot
# default vs --untangle layouts (--mode spaced).
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
HERE="$(cd "$(dirname "$0")" && pwd)"
DATA="$ROOT/test/data/data"
RUNS="$HERE/runs/yeast4"
MUMENTO2="${MUMENTO2:-/hps/nobackup/rdf/research/vshiv/mumemto/kmer_mumemto/target/release/mumento2}"
VIZ="$ROOT/mumemto/viz_mums.py"
MIN_LEN="${MIN_LEN:-20}"
THREADS="${THREADS:-8}"
SPACER="${SPACER:-0.1}"

FASTA=(
  "$DATA/GCA_000146045.2_R64_genomic.fa"
  "$DATA/GCA_003086655.1_ASM308665v1_genomic.fa"
  "$DATA/GCA_902192335.1_UWOPS83-787.3_genomic.fa"
  "$DATA/GCA_949124045.1_AGK.nuclear_genome.ScRAP_genomic.fa"
)

LABELS="$RUNS/labels.txt"
mkdir -p "$RUNS"

cat > "$LABELS" <<'EOF'
R64
ASM308665v1
UWOPS83-787
AGK
EOF

if [[ ! -x "$MUMENTO2" ]]; then
  echo "mumento2 not found/executable: $MUMENTO2" >&2
  exit 1
fi

OUT_ZST="$RUNS/yeast4.tsv.zst"
OUT_MUMS="$RUNS/yeast4.mums"
OUT_LENS="$RUNS/yeast4.lengths"

if [[ ! -f "$OUT_MUMS" ]]; then
  if [[ ! -f "$OUT_ZST" ]]; then
    echo "==> mumento2 --contig-coords (min-len=$MIN_LEN) ..."
    "$MUMENTO2" \
      --output "$OUT_ZST" \
      --contig-coords \
      -l "$MIN_LEN" \
      -j "$THREADS" \
      "${FASTA[@]}"
  else
    echo "==> reusing existing $OUT_ZST"
  fi
  STEM_LENS="${OUT_ZST%.tsv.zst}.lengths"
  if [[ -f "$STEM_LENS" && "$STEM_LENS" != "$OUT_LENS" ]]; then
    mv -f "$STEM_LENS" "$OUT_LENS"
  elif [[ ! -f "$OUT_LENS" ]]; then
    LATEST=$(ls -t "$RUNS"/*.lengths 2>/dev/null | head -1 || true)
    if [[ -n "${LATEST:-}" ]]; then
      cp -f "$LATEST" "$OUT_LENS"
    else
      echo "Could not find .lengths next to mumento2 output" >&2
      ls -la "$RUNS" >&2
      exit 1
    fi
  fi
  echo "==> decompress MUMs to plain text"
  if [[ -f "$OUT_ZST" ]]; then
    zstd -dc "$OUT_ZST" > "$OUT_MUMS"
  else
    echo "Expected zstd MUM output at $OUT_ZST" >&2
    ls -la "$RUNS" >&2
    exit 1
  fi
else
  echo "==> reusing existing $OUT_MUMS"
fi

if [[ ! -f "$OUT_LENS" ]]; then
  echo "Missing lengths: $OUT_LENS" >&2
  exit 1
fi

run_viz() {
  local name="$1"
  local out="$2"
  shift 2
  echo "==> viz $* -> $out"
  python3 "$VIZ" \
    -m "$OUT_MUMS" \
    -l "$OUT_LENS" \
    -n "$LABELS" \
    --mode spaced \
    --spacer "$SPACER" \
    "$@" \
    -o "$out" \
    --dpi 200 \
    --dims 10 6 \
    -v
}

run_viz default "$RUNS/default.png"
run_viz untangle "$RUNS/untangle.png" --untangle

echo "==> done"
ls -la "$RUNS"/*.png 2>/dev/null || true
