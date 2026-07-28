#!/usr/bin/env bash
###############################################################################
# Focused micro-synteny of a reference region across genomes, with one gene's
# synteny ribbons highlighted in red (jcvi.graphics.synteny).
#
# Prereq: pairwise homology files <REF>.<G>.last (blast/last -outfmt 6) for each
# other genome, e.g. produced by run_synteny.sh (jcvi writes ecd.IHx.last), and
# the per-genome <name>.bed / (all in DATADIR).
###############################################################################
set -euo pipefail
SCRIPTDIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# -------- CONFIG --------
REF=ecd
OTHERS="IH1,IH2,IH3,IH4,IH5"
CHROM=A08
START=14935176
END=15537537
TARGET=BraA08t039305E          # gene whose synteny lines are drawn red
HLCOLOR=red
LABELS="ECD04 A08,IH1,IH2,IH3,IH4,IH5"
DATADIR="$(pwd)"               # holds <name>.bed
ALIGNDIR="$(pwd)"              # holds <REF>.<other>.last
OUT=microsynteny_${CHROM}
# ------------------------

mkdir -p "$OUT"; cd "$OUT"
for g in $REF $(echo "$OTHERS"|tr ',' ' '); do ln -sf "$DATADIR/$g.bed" "$g.bed"; done
cat $REF.bed $(echo "$OTHERS"|tr ',' ' '|sed 's/[^ ]*/&.bed/g') > all.bed

python3 "$SCRIPTDIR/build_region_blocks.py" --ref "$REF" --others "$OTHERS" \
  --chrom "$CHROM" --start "$START" --end "$END" --target "$TARGET" \
  --hlcolor "$HLCOLOR" --datadir "$DATADIR" --aligndir "$ALIGNDIR" --out region.blocks

# build layout: stacked tracks, consecutive edges
python3 - "$LABELS" > blocks.layout <<'PY'
import sys
labels=sys.argv[1].split(',')
n=len(labels); ys=[0.93-(0.85*i/(n-1)) for i in range(n)]
print("# x, y, rotation, ha, va, color, ratio, label")
for i,l in enumerate(labels):
    print(f"0.5, {ys[i]:.3f}, 0, left, center, , 1, {l}")
print("# edges")
for i in range(n-1): print(f"e, {i}, {i+1}")
PY

python3 -m jcvi.graphics.synteny region.blocks all.bed blocks.layout \
  --genelabels "$TARGET" --genelabelsize 6 --figsize 11x9 --notex --format pdf
python3 -m jcvi.graphics.synteny region.blocks all.bed blocks.layout \
  --genelabels "$TARGET" --genelabelsize 6 --figsize 11x9 --notex --format png --dpi 220
mv region.png "microsynteny_${CHROM}_${TARGET}.png"
mv region.pdf "microsynteny_${CHROM}_${TARGET}.pdf"
echo "DONE -> $OUT/microsynteny_${CHROM}_${TARGET}.png"
