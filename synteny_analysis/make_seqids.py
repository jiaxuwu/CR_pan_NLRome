#!/usr/bin/env python3
"""
Build jcvi `seqids` + `layout` (+ optional relabeled BEDs) for a stacked
karyotype ribbon plot across an ordered list of genomes.

Track 0 (reference) uses a fixed seqid list (e.g. pseudo-chromosomes).
Each next track keeps the top-N sequences by anchored genes to the PREVIOUS
track and orders them to follow that track (fewer ribbon crossings).

With --relabel, chosen sequences are renamed to short display tags and copies
of each BED are written (<track>.kb.bed). Use jcvi karyotype --keep-chrlabels
so those short tags show as clean labels instead of long contig names.
"""
import sys, os, argparse
from collections import defaultdict

def read_bed(path):
    order={}; seqgenes=defaultdict(list); rows=[]
    with open(path) as f:
        for ln in f:
            if not ln.strip() or ln.startswith('#'): continue
            c=ln.rstrip('\n').split('\t')
            seq,start,end,name=c[0],int(c[1]),int(c[2]),c[3]
            seqgenes[seq].append((start,name)); rows.append(c)
    for seq,lst in seqgenes.items():
        lst.sort()
        for i,(s,n) in enumerate(lst): order[n]=(seq,i)
    return order, rows

def read_simple(path):
    with open(path) as f:
        for ln in f:
            if ln.startswith('StartGeneA') or ln.startswith('#') or not ln.strip(): continue
            c=ln.rstrip('\n').split('\t')
            if len(c)<6: continue
            try: sc=int(c[5])
            except: sc=1
            yield c[0],c[2],sc   # geneA(prev), geneB(this), score

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--tracks',required=True)
    ap.add_argument('--refseqs',required=True)
    ap.add_argument('--reflabels',default='',help='display labels for ref seqs (default: same)')
    ap.add_argument('--topn',type=int,default=18)
    ap.add_argument('--minblocks',type=int,default=3)
    ap.add_argument('--beddir',default='.')
    ap.add_argument('--simpledir',default='.')
    ap.add_argument('--out_seqids',default='seqids')
    ap.add_argument('--out_layout',default='layout')
    ap.add_argument('--labels',default='')
    ap.add_argument('--relabel',action='store_true')
    a=ap.parse_args()

    tracks=a.tracks.split(',')
    labels=a.labels.split(',') if a.labels else tracks
    beds={}; rows={}
    for t in tracks:
        o,r=read_bed(os.path.join(a.beddir,f'{t}.bed')); beds[t]=o; rows[t]=r

    chosen={tracks[0]:a.refseqs.split(',')}
    prev_rank={s:i for i,s in enumerate(chosen[tracks[0]])}
    for k in range(1,len(tracks)):
        P,T=tracks[k-1],tracks[k]
        bp,bt=beds[P],beds[T]
        agg=defaultdict(lambda:[0,0.0,0])
        for gA,gB,sc in read_simple(os.path.join(a.simpledir,f'{P}.{T}.simple')):
            if gA not in bp or gB not in bt: continue
            pseq,prank=bp[gA]
            if pseq not in prev_rank: continue
            tseq,_=bt[gB]; pos=prev_rank[pseq]+prank/1e6
            g=agg[tseq]; g[0]+=sc; g[1]+=pos*sc; g[2]+=1
        cand=[(ts,v[0],v[1]/v[0]) for ts,v in agg.items() if v[2]>=a.minblocks]
        cand.sort(key=lambda x:-x[1]); cand=cand[:a.topn]; cand.sort(key=lambda x:x[2])
        chosen[T]=[c[0] for c in cand]
        prev_rank={s:i for i,s in enumerate(chosen[T])}

    # display names
    disp={}
    reflab=a.reflabels.split(',') if a.reflabels else chosen[tracks[0]]
    for i,s in enumerate(chosen[tracks[0]]): disp[(tracks[0],s)]=reflab[i]
    for t in tracks[1:]:
        for i,s in enumerate(chosen[t]): disp[(t,s)]=str(i+1)

    bedfield={t:f'{t}.bed' for t in tracks}
    if a.relabel:
        for t in tracks:
            outb=os.path.join(a.beddir,f'{t}.kb.bed'); bedfield[t]=f'{t}.kb.bed'
            with open(outb,'w') as o:
                for c in rows[t]:
                    seq=c[0]; c=c[:]; c[0]=disp.get((t,seq),seq); o.write('\t'.join(c)+'\n')

    with open(a.out_seqids,'w') as o:
        for t in tracks:
            names=[disp[(t,s)] for s in chosen[t]] if a.relabel else chosen[t]
            o.write(','.join(names)+'\n')

    n=len(tracks); ys=[0.92-(0.86*i/(n-1)) for i in range(n)]
    with open(a.out_layout,'w') as o:
        o.write("# y, xstart, xend, rotation, color, label, va, bed\n")
        for i,t in enumerate(tracks):
            o.write(f"{ys[i]:.3f}, 0.05, 0.95, 0, , {labels[i]}, top, {bedfield[t]}\n")
        o.write("# edges\n")
        for i in range(n-1):
            o.write(f"e, {i}, {i+1}, {tracks[i]}.{tracks[i+1]}.simple\n")
    for t in tracks: print(f"{t}: {len(chosen[t])} seqs")
if __name__=='__main__': main()
