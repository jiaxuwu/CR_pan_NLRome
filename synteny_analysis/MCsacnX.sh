#!/usr/bin/env bash

# Global micro-synteny pipeline:  ECD04 (reference) vs IH1..IH5
#
#   Input  (per genome, in $DATADIR):   <name>.bed   <name>.cds
#   Tools  :  jcvi (python -m jcvi...),  LAST (lastdb/lastal),  MCScanX
#   Output :  - jcvi ribbon/karyotype figure (ECD04 + IH1..IH5, adjacent ribbons)
#             - jcvi pairwise dot plots + .anchors/.simple
#             - MCScanX reference-based .collinearity (ECD04 vs each IH)
#             - collinearity_summary.csv

set -euo pipefail

# ----------------------------- CONFIG ---------------------------------------
GENOMES=(ecd IH1 IH2 IH3 IH4 IH5)     # stacking order, top -> bottom
LABELS=(ECD04 IH1 IH2 IH3 IH4 IH5)    # display labels (same length as GENOMES)

# Reference chromosomes shown on the top track (ECD04 pseudo-chromosomes)
REF_SEQS="A01,A02,A03,A04,A05,A06,A07,A08,A09,A10"

ALIGNER=last          # 'last' (uses your LAST install, works on .cds directly)
                      # or 'diamond_blastp' (needs proteins; script will translate)
DBTYPE=nucl           # nucl for .cds ; auto-set to prot if ALIGNER=diamond_blastp
CPUS=8                # CPUs for jcvi / LAST
CSCORE=0.7            # jcvi c-score (0.7 default; raise toward 0.99 for stricter)

TOPN=18               # IH contigs displayed per track in the figure
MINBLOCKS=3           # min syntenic blocks for a contig to be displayed

DATADIR="$(pwd)"      # directory holding <name>.bed and <name>.cds
OUTDIR="synteny_out"
SCRIPTDIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# ----------------------------------------------------------------------------

REF="${GENOMES[0]}"
mkdir -p "$OUTDIR"/{work,mcscanx,figures}
cd "$OUTDIR/work"

echo "==> 0. Staging inputs and sanity checks"
for g in "${GENOMES[@]}"; do
  [[ -s "$DATADIR/$g.bed" ]] || { echo "MISSING $g.bed"; exit 1; }
  [[ -s "$DATADIR/$g.cds" ]] || { echo "MISSING $g.cds"; exit 1; }
  ln -sf "$DATADIR/$g.bed" "$g.bed"
  ln -sf "$DATADIR/$g.cds" "$g.cds"
  echo "    $g: $(cut -f4 "$g.bed"|sort -u|wc -l) bed genes, $(grep -c '^>' "$g.cds") cds seqs"
done

SEQEXT=cds
if [[ "$ALIGNER" == "diamond_blastp" ]]; then
  DBTYPE=prot; SEQEXT=pep
  for g in "${GENOMES[@]}"; do
    [[ -s "$g.pep" ]] || python3 "$SCRIPTDIR/translate_cds.py" "$g.cds" "$g.pep"
  done
fi

# adjacent pairs -> stacked ribbon figure
ADJ=(); for ((i=0;i<${#GENOMES[@]}-1;i++)); do ADJ+=("${GENOMES[$i]}:${GENOMES[$((i+1))]}"); done
# reference pairs -> MCScanX (ECD04 is chromosome-scale => fast & meaningful)
REFP=(); for ((i=1;i<${#GENOMES[@]};i++)); do REFP+=("$REF:${GENOMES[$i]}"); done

run_ortholog () {  # $1=A $2=B  -> A.B.anchors, A.B.simple, A.B.last, A.B.pdf
  local A=$1 B=$2
  if [[ ! -s "$A.$B.anchors" ]]; then
    echo "    jcvi ortholog $A $B  (aligner=$ALIGNER)"
    python3 -m jcvi.compara.catalog ortholog "$A" "$B" \
      --dbtype "$DBTYPE" --align_soft "$ALIGNER" --no_strip_names \
      --cscore "$CSCORE" --notex --cpus "$CPUS"
  fi
  [[ -s "$A.$B.simple" ]] || \
    python3 -m jcvi.compara.synteny simple "$A.$B.anchors" --qbed="$A.bed" --sbed="$B.bed"
}

echo "==> 1. jcvi synteny for adjacent pairs (ribbon figure)"
for p in "${ADJ[@]}";  do run_ortholog "${p%:*}" "${p#*:}"; done
echo "==> 2. jcvi synteny for reference pairs (for MCScanX)"
for p in "${REFP[@]}"; do run_ortholog "${p%:*}" "${p#*:}"; done

echo "==> 3. Build seqids + layout and draw the karyotype ribbon plot"
python3 "$SCRIPTDIR/make_seqids.py" \
  --tracks "$(IFS=,;echo "${GENOMES[*]}")" \
  --refseqs "$REF_SEQS" --reflabels "$REF_SEQS" \
  --labels "$(IFS=,;echo "${LABELS[*]}")" \
  --topn "$TOPN" --minblocks "$MINBLOCKS" --relabel
python3 -m jcvi.graphics.karyotype seqids layout --keep-chrlabels --notex \
  --format pdf -o ../figures/karyotype_ECD04_IH1-5.pdf
python3 -m jcvi.graphics.karyotype seqids layout --keep-chrlabels --notex \
  --format png --dpi 300 -o ../figures/karyotype_ECD04_IH1-5.png

echo "==> 4. MCScanX reference-based collinearity (ECD04 vs each IH)"
declare -A TAG=( [ecd]=EA [IH1]=IA [IH2]=IB [IH3]=IC [IH4]=ID [IH5]=IE )
for ((i=1;i<${#GENOMES[@]};i++)); do
  B="${GENOMES[$i]}"; d="../mcscanx/${REF}_${B}"; mkdir -p "$d"
  python3 "$SCRIPTDIR/prep_mcscanx.py" --tracks "$REF,$B" \
      --tags "${TAG[$REF]},${TAG[$B]}" --pairs "$REF:$B" \
      --beddir . --lastdir . --out "$d/${REF}_${B}"
  ( cd "$d" && MCScanX "${REF}_${B}" >/dev/null 2>&1; find . -name '*.html' -delete )
  echo "    ECD04 vs $B : $(grep -c '## Alignment' "$d/${REF}_${B}.collinearity") blocks"
done

echo "==> 5. Summary table"
python3 "$SCRIPTDIR/synteny_summary.py" --ref "$REF" \
  --others "$(IFS=,;echo "${GENOMES[*]:1}")" --dir ../mcscanx \
  --reftotal "$(cut -f4 $REF.bed|sort -u|wc -l)" \
  --out ../collinearity_summary.csv

echo "==> DONE. See $OUTDIR/figures/ , $OUTDIR/mcscanx/ , $OUTDIR/collinearity_summary.csv"
