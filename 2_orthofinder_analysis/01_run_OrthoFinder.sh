#!/usr/bin/env bash
# OrthoFinder on the six NBARC-filtered NLR proteomes.
#
# Six input fastas in ../proteomes/:
#   IH1.fa  IH2.fa  IH3.fa  IH4.fa  IH5.fa  Westar.fa
#
# Recommended parameters for divergent gene families (NLRs):
#   -S diamond_ultra_sens   sensitive Diamond search
#   -M msa                  MSA-based gene trees (more accurate than dendroblast)
#   -A mafft                MAFFT for alignment
#   -T fasttree             FastTree for gene trees
#   -t / -a                 cpu threads (auto-detect below)
#
set -euo pipefail

SCRIPT_DIR="/path/to/new_pan-NLRome/orthofinder/scripts"
BASE_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
INPUT_DIR="$BASE_DIR/proteomes"
OUT_DIR="$BASE_DIR/orthofinder_out"


orthofinder \
    -f "$INPUT_DIR" \
    -o "$OUT_DIR" \
    -S diamond_ultra_sens \
    -M msa \
    -A mafft \
    -T fasttree \
    -I 1.2 \
    -t 20 \
    -a 20 \
    -n pan_nlrome

echo
echo ">>> OrthoFinder finished."
echo ">>> Results in: $OUT_DIR/Results_pan_nlrome"
echo
echo ">>> Next: bash 02_postprocess.sh  (or run 02_postprocess.py directly)"
