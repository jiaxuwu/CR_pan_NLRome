#!/usr/bin/env python3
"""
Build MCScanX inputs (<out>.gff + <out>.blast) for a GLOBAL multi-genome run.

Each genome gets a unique 2-letter TAG prepended to its sequence names so
MCScanX can separate genomes (it groups by the first 2 characters). Gene names
are left unchanged so they match the BLAST/last hit IDs.

BLAST file = concatenation of the supplied pairwise .last files (blast -outfmt 6);
each hit is also written in reverse orientation so anchors are symmetric.
"""
import sys, os, argparse

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--tracks',required=True)          # ecd,IH1,...
    ap.add_argument('--tags',required=True)            # EA,IA,IB,IC,ID,IE
    ap.add_argument('--pairs',required=True)           # ecd:IH1,IH1:IH2,...
    ap.add_argument('--beddir',default='.')
    ap.add_argument('--lastdir',default='.')
    ap.add_argument('--out',default='global')
    a=ap.parse_args()
    tracks=a.tracks.split(','); tags=a.tags.split(',')
    tag={t:g for t,g in zip(tracks,tags)}

    with open(a.out+'.gff','w') as o:
        for t in tracks:
            with open(os.path.join(a.beddir,f'{t}.bed')) as f:
                for ln in f:
                    if not ln.strip() or ln.startswith('#'): continue
                    c=ln.rstrip('\n').split('\t')
                    seq,start,end,name=c[0],c[1],c[2],c[3]
                    o.write(f"{tag[t]}{seq}\t{name}\t{start}\t{end}\n")

    n=0
    with open(a.out+'.blast','w') as o:
        for pr in a.pairs.split(','):
            A,B=pr.split(':')
            lf=os.path.join(a.lastdir,f'{A}.{B}.last')
            with open(lf) as f:
                for ln in f:
                    c=ln.rstrip('\n').split('\t')
                    if len(c)<12: continue
                    o.write('\t'.join(c)+'\n'); n+=1
                    r=[c[1],c[0],c[2],c[3],c[4],c[5],c[8],c[9],c[6],c[7],c[10],c[11]]
                    o.write('\t'.join(r)+'\n'); n+=1
    print(f"wrote {a.out}.gff and {a.out}.blast ({n} hit lines)")
if __name__=='__main__': main()
