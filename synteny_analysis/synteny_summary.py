#!/usr/bin/env python3
"""Summarise MCScanX .collinearity files into a CSV table."""
import sys,os,re,csv,argparse
ap=argparse.ArgumentParser()
ap.add_argument('--ref',default='ecd')
ap.add_argument('--others',required=True)          # IH1,IH2,...
ap.add_argument('--dir',default='.')
ap.add_argument('--refprefix',default='Bra',help='prefix identifying reference gene IDs')
ap.add_argument('--reftotal',type=int,required=True)
ap.add_argument('--out',default='collinearity_summary.csv')
a=ap.parse_args()
rows=[]
for B in a.others.split(','):
    f=os.path.join(a.dir,f'{a.ref}_{B}',f'{a.ref}_{B}.collinearity')
    blocks=pairs=0; eg=set()
    with open(f) as fh:
        for ln in fh:
            if ln.startswith('## Alignment'): blocks+=1; continue
            c=ln.split('\t')
            if len(c)>=3 and re.match(r'\s*\d+-\s*\d+',c[0]):
                pairs+=1
                for g in (c[1].strip(),c[2].strip()):
                    if g.startswith(a.refprefix): eg.add(g)
    rows.append([f'{a.ref} vs {B}',blocks,pairs,len(eg),round(100*len(eg)/a.reftotal,1)])
with open(a.out,'w',newline='') as o:
    w=csv.writer(o); w.writerow(['comparison','collinear_blocks','anchor_pairs','ref_genes_in_synteny','pct_ref_genes'])
    w.writerows(rows)
print(open(a.out).read())
