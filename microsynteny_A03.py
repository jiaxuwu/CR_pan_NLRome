#!/usr/bin/env python3
"""
Standalone micro-synteny figure - ECD04 A03 locus vs IH1-IH5, target BraA03t012133E.

Draws: ECD04 A03 region full-width in the middle; IH1/IH2 side-by-side on top;
IH3/IH4/IH5 side-by-side on the bottom; each isolate joined independently to ECD04.
Genes are arrows (forward #FFF7C5, reverse #67A2C5, thin black border); the target
gene's ribbons are highlighted #E05454. No text labels.

INPUTS (all in DATADIR, default = current folder):
  ecd.bed, ecd.cds  and  IH1.bed/IH1.cds ... IH5.bed/IH5.cds
    .bed = 6-col BED (chrom,start,end,geneID,score,strand); .cds = nucleotide FASTA
    with headers == BED column 4.
REQUIRES: python with `jcvi`, and LAST (lastdb/lastal) on PATH. Pairwise
  ecd.<IH>.last files are auto-generated with jcvi+LAST if missing (set
  ALIGNER="diamond_blastp" to use DIAMOND on translated proteins instead).

RUN:    python microsynteny_A03.py
OUTPUT: microsynteny_A03.png / .pdf
"""
import os, sys, subprocess, matplotlib
matplotlib.use("Agg")

# ============================== CONFIG ======================================
REF     = "ecd"
OTHERS  = ["IH1", "IH2", "IH3", "IH4", "IH5"]
CHROM, START, END = "A03", 24949048, 25141642
TARGET  = "BraA03t012133E"

# display window per isolate = (contig, start, end) in ACTUAL contig coordinates.
# IH2/IH3/IH4 blocks are reverse-oriented vs ECD04, so values are the reverse-
# converted display coordinates (contig_length - coord). Use end=None for whole contig.
WINDOWS = {
    "IH1": ("contig_4441_polished", 913793, 1097000),
    "IH2": ("contig_1104_polished", 261488, 445488),   # display 1.216-1.4 M reversed
    "IH3": ("contig_4607_polished", 0,      132678),   # display 113.952-246.612 K reversed
    "IH4": ("contig_1572_polished", 62530,  249530),   # display 1.129-1.316 M reversed
    "IH5": ("contig_3416_polished", 0,      None),     # whole contig
}
TOP_COLS, BOTTOM_COLS = [1, 2], [3, 4, 5]     # block columns (0 = reference)

ALIGNER   = "last"          # "last" or "diamond_blastp"
DBTYPE    = "nucl"          # nucl for .cds ; prot for diamond_blastp
FWD_COLOR, REV_COLOR = "#FFF7C5", "#67A2C5"
GENE_EDGE, GENE_LW   = "#000000", 0.5
RIBBON, RIBBON_ALPHA = "#E05454", 0.8         # alpha: 1=solid, 0.2=80% transparent
FIGSIZE, DPI = (12, 8), 220
DATADIR   = "."
OUTPREFIX = "microsynteny_A03"
CANVAS    = 0.65            # jcvi internal constant; leave as-is
# ============================================================================

_CODON = {
 'TTT':'F','TTC':'F','TTA':'L','TTG':'L','CTT':'L','CTC':'L','CTA':'L','CTG':'L',
 'ATT':'I','ATC':'I','ATA':'I','ATG':'M','GTT':'V','GTC':'V','GTA':'V','GTG':'V',
 'TCT':'S','TCC':'S','TCA':'S','TCG':'S','CCT':'P','CCC':'P','CCA':'P','CCG':'P',
 'ACT':'T','ACC':'T','ACA':'T','ACG':'T','GCT':'A','GCC':'A','GCA':'A','GCG':'A',
 'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','CAT':'H','CAC':'H','CAA':'Q','CAG':'Q',
 'AAT':'N','AAC':'N','AAA':'K','AAG':'K','GAT':'D','GAC':'D','GAA':'E','GAG':'E',
 'TGT':'C','TGC':'C','TGA':'*','TGG':'W','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
 'AGT':'S','AGC':'S','AGA':'R','AGG':'R','GGT':'G','GGC':'G','GGA':'G','GGG':'G'}
def _translate(inf, outf):
    def tr(s):
        s=s.upper().replace('U','T')
        p=''.join(_CODON.get(s[i:i+3],'X') for i in range(0,len(s)-2,3))
        return (p[:-1] if p.endswith('*') else p).replace('*','X')
    with open(inf) as f, open(outf,'w') as o:
        hid=None; seq=[]
        def flush():
            if hid and seq:
                p=tr(''.join(seq))
                if p: o.write(">"+hid+"\n"+"\n".join(p[i:i+60] for i in range(0,len(p),60))+"\n")
        for ln in f:
            if ln.startswith('>'): flush(); hid=ln[1:].split()[0].strip(); seq=[]
            else: seq.append(ln.strip())
        flush()

def ensure_last():
    if ALIGNER == "diamond_blastp":
        for g in [REF]+OTHERS:
            pep=os.path.join(DATADIR,f"{g}.pep")
            if not os.path.exists(pep): _translate(os.path.join(DATADIR,f"{g}.cds"), pep)
    for g in OTHERS:
        lf=os.path.join(DATADIR,f"{REF}.{g}.last")
        if os.path.exists(lf) and os.path.getsize(lf)>0: continue
        print(f"[align] {REF} vs {g} ({ALIGNER}) ...")
        subprocess.run([sys.executable,"-m","jcvi.compara.catalog","ortholog",REF,g,
                        "--dbtype",DBTYPE,"--align_soft",ALIGNER,"--no_strip_names",
                        "--notex","--cpus","4"], cwd=DATADIR, check=True)

def build_blocks(path):
    contig={}
    for g in OTHERS:
        for ln in open(os.path.join(DATADIR,f"{g}.bed")):
            c=ln.rstrip("\n").split("\t"); contig[c[3]]=c[0]
    region=[]
    for ln in open(os.path.join(DATADIR,f"{REF}.bed")):
        c=ln.rstrip("\n").split("\t")
        if c[0]==CHROM and START<=int(c[1]) and int(c[2])<=END: region.append((int(c[1]),c[3]))
    region.sort(); region_genes=[g for _,g in region]; rset=set(region_genes)
    genes_on={}
    for g in OTHERS:
        ctg,lo,hi=WINDOWS[g]; hi=10**12 if hi is None else hi; on=[]
        for ln in open(os.path.join(DATADIR,f"{g}.bed")):
            c=ln.rstrip("\n").split("\t")
            if c[0]==ctg and lo<=int(c[1])<=hi: on.append((int(c[1]),c[3]))
        on.sort(); genes_on[g]=[n for _,n in on]
    winset={g:set(genes_on[g]) for g in OTHERS}
    best={g:{} for g in OTHERS}
    for g in OTHERS:
        for ln in open(os.path.join(DATADIR,f"{REF}.{g}.last")):
            c=ln.split("\t"); q=c[0]
            if q not in rset or c[1] not in winset[g]: continue
            bit=float(c[11])
            if q not in best[g] or bit>best[g][q][1]: best[g][q]=(c[1],bit)
    ncol=len(OTHERS); used={g:set() for g in OTHERS}; rows=[]
    for eg in region_genes:
        row=[eg]
        for g in OTHERS:
            v=best[g].get(eg,(".",0))[0]; row.append(v)
            if v!=".": used[g].add(v)
        rows.append(row)
    ctx=[]
    for gi,g in enumerate(OTHERS):
        for gene in genes_on[g]:
            if gene in used[g]: continue
            r=["."]*(ncol+1); r[gi+1]=gene; ctx.append(r)
    with open(path,"w") as o:
        for row in rows:
            line="\t".join(row)
            if row[0]==TARGET: line=f"{RIBBON}*"+line
            o.write(line+"\n")
        for r in ctx: o.write("\t".join(r)+"\n")
    for g in OTHERS: print(f"  {g}: window genes={len(genes_on[g])}, matched={len(used[g])}")

def build_layout(blocks_path, bed_path, out_path):
    from jcvi.formats.bed import Bed
    from jcvi.compara.synteny import BlockFile
    bed=Bed(bed_path); order=bed.order; bf=BlockFile(blocks_path)
    spans={i:bf.get_extent(i,order,debug=False)[-1] for i in range(bf.ncols)}
    maxspan=max(spans.values()); natw=lambda i: spans[i]/maxspan*CANVAS
    def row(cols,y,xlo=0.07,xhi=0.93,gap=0.04):
        n=len(cols); w=((xhi-xlo)-gap*(n-1))/n
        return [(c, xlo+w/2+k*(w+gap), y, w/natw(c)) for k,c in enumerate(cols)]
    placed={0:(0.5,0.5,0.86/natw(0))}
    for c,cx,y,r in row(TOP_COLS,0.9):    placed[c]=(cx,y,r)
    for c,cx,y,r in row(BOTTOM_COLS,0.1): placed[c]=(cx,y,r)
    with open(out_path,"w") as o:
        o.write("# x, y, rotation, ha, va, color, ratio, label\n")
        for i in range(bf.ncols):
            cx,y,r=placed[i]; o.write(f"{cx:.3f}, {y:.3f}, 0, center, center, k, {r:.3f}, .\n")
        o.write("# edges\n")
        for i in range(1,bf.ncols): o.write(f"e, 0, {i}\n")

def render(blocks_path, bed_path, layout_path):
    from matplotlib import rcParams
    rcParams["text.usetex"]=False
    rcParams["font.family"]=["Arial","Liberation Sans","DejaVu Sans"]
    import matplotlib.pyplot as plt
    import jcvi.graphics.glyph as G
    G.OrientationPalette.forward=FWD_COLOR
    G.OrientationPalette.backward=REV_COLOR
    G.OrientationPalette.palette={"+":FWD_COLOR,"-":REV_COLOR}
    import jcvi.compara.synteny as CS, jcvi.graphics.synteny as SY
    _Block=CS.BlockFile
    class BlockFileHex(_Block):
        def __init__(self, filename, defaultcolor="#fb8072", header=False, allow_hex=True):
            super().__init__(filename, defaultcolor=defaultcolor, header=header, allow_hex=True)
    SY.BlockFile=BlockFileHex
    _Glyph=G.Glyph
    def _GlyphBorder(ax,x1,x2,y,height=0.04,gradient=True,fc="gray",ec="gainsboro",lw=0,style="box",**kw):
        return _Glyph(ax,x1,x2,y,height=height,gradient=gradient,fc=fc,ec=GENE_EDGE,lw=GENE_LW,style=style,**kw)
    SY.Glyph=_GlyphBorder
    _Shade=SY.Shade
    def _ShadeAlpha(ax,a,b,ymid_pad=0.0,highlight=False,alpha=RIBBON_ALPHA,**kw):
        if highlight: alpha=RIBBON_ALPHA
        return _Shade(ax,a,b,ymid_pad,highlight=highlight,alpha=alpha,**kw)
    SY.Shade=_ShadeAlpha
    from jcvi.graphics.synteny import Synteny
    fig=plt.figure(figsize=FIGSIZE); root=fig.add_axes((0,0,1,1))
    Synteny(fig,root,blocks_path,bed_path,layout_path,loc_label=True,gene_labels=None,genelabelsize=0)
    for t in list(root.texts): t.remove()
    root.set_xlim(0,1); root.set_ylim(0,1); root.set_axis_off()
    for ext in ("png","pdf"):
        fig.savefig(f"{OUTPREFIX}.{ext}",dpi=DPI,bbox_inches="tight"); print("saved",f"{OUTPREFIX}.{ext}")

def main():
    ensure_last()
    all_bed=os.path.join(DATADIR,"_all.bed")
    with open(all_bed,"w") as o:
        for g in [REF]+OTHERS: o.write(open(os.path.join(DATADIR,f"{g}.bed")).read())
    blocks=os.path.join(DATADIR,f"{OUTPREFIX}.blocks"); layout=os.path.join(DATADIR,f"{OUTPREFIX}.layout")
    build_blocks(blocks); build_layout(blocks,all_bed,layout); render(blocks,all_bed,layout)

if __name__=="__main__":
    main()
