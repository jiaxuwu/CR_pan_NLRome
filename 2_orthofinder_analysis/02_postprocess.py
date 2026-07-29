#!/usr/bin/env python3
"""
Pan-NLRome downstream analyses on OrthoFinder output.

"""
from __future__ import annotations
import os, sys, glob, itertools, random, argparse
from collections import Counter, defaultdict

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---------- paths ----------
HERE       = os.path.dirname(os.path.abspath(__file__))
BASE       = os.path.abspath(os.path.join(HERE, ".."))

# CLI: --of-dir <orthofinder Results_* dir>   --results-dir <where to write output>
_pp = argparse.ArgumentParser(add_help=False)
_pp.add_argument("--of-dir",      default=os.path.join(BASE, "orthofinder_out", "Results_pan_nlrome"),
                 help="Path to OrthoFinder Results_<name>/ directory")
_pp.add_argument("--results-dir", default=os.path.join(BASE, "results"),
                 help="Where to write postprocess tables/figures")
_ppargs, _ = _pp.parse_known_args()

RES_DIR    = _ppargs.results_dir
FIG_DIR    = os.path.join(RES_DIR, "figures")
os.makedirs(FIG_DIR, exist_ok=True)

OF_DIR     = _ppargs.of_dir
GC_FILE    = os.path.join(OF_DIR, "Orthogroups", "Orthogroups.GeneCount.tsv")
OG_FILE    = os.path.join(OF_DIR, "Orthogroups", "Orthogroups.tsv")
UA_FILE    = os.path.join(OF_DIR, "Orthogroups", "Orthogroups_UnassignedGenes.tsv")
print(f">>> OrthoFinder dir: {OF_DIR}")
print(f">>> Writing results to: {RES_DIR}")

# NLRtracker annotation files: pass via --tracker-h1 etc. or autodiscover
DEFAULT_TRACKER = {
    # Relative to scripts/ -> IH1-IH5_NLRtracker/<H?>/NLRtracker.tsv  (../../.. = IH1-IH5_NLRtracker)
    "H1":     "../../../IH1/NLRtracker.tsv",
    "H2":     "../../../IH2/NLRtracker.tsv",
    "H3":     "../../../IH3/NLRtracker.tsv",
    "H4":     "../../../IH4/NLRtracker.tsv",
    "H5":     "../../../IH5/NLRtracker.tsv",
    # westar_reannotation is sibling of H1-H5_NLRtracker
    "Westar": "../../../../westar_reannotation/NLRtracker_final/NLRtracker.tsv",
}

# Override DEFAULT_TRACKER from CLI:  --tracker H1=/path/to/NLRtracker.tsv ...
_ap = argparse.ArgumentParser(add_help=False)
_ap.add_argument("--tracker", action="append", default=[])
_args, _ = _ap.parse_known_args()
for kv in _args.tracker:
    if "=" in kv:
        k, v = kv.split("=", 1)
        DEFAULT_TRACKER[k] = v

RESISTANT  = ["H1","H2","H3","H4","H5"]
SUSCEPT    = ["Westar"]
LINES      = RESISTANT + SUSCEPT

# ---------- load orthogroups ----------
print(">>> Loading OrthoFinder output")
gc = pd.read_csv(GC_FILE, sep="\t").set_index("Orthogroup")
# Drop "Total" column if present
if "Total" in gc.columns:
    gc = gc.drop(columns=["Total"])
# Reorder lines
gc = gc[[c for c in LINES if c in gc.columns]]
print(f"    {gc.shape[0]} orthogroups, {gc.shape[1]} lines")

# Also load unassigned genes (singletons) and merge as 1-member OGs
if os.path.exists(UA_FILE):
    ua = pd.read_csv(UA_FILE, sep="\t").set_index("Orthogroup")
    ua_gc = ua.notna().astype(int)
    ua_gc = ua_gc[[c for c in LINES if c in ua_gc.columns]]
    gc = pd.concat([gc, ua_gc])
    print(f"    + {ua_gc.shape[0]} singletons -> {gc.shape[0]} total OGs")

# ---------- presence/absence ----------
print(">>> Building presence/absence matrix")
pa = (gc > 0).astype(int)
pa.to_csv(os.path.join(RES_DIR, "presence_absence.tsv"), sep="\t")
gc.to_csv(os.path.join(RES_DIR, "gene_counts.tsv"), sep="\t")

# ---------- pan-genome category per orthogroup ----------
print(">>> Classifying pan-genome categories")
N = pa.shape[1]
n_lines = pa.sum(axis=1)

def category(n, N):
    """4-category pan scheme (core/soft_core/shell/private).
    For small N (~6 lines) this is more informative than the 5-band Tettelin
    split because 'shell' meaningfully spans n=2..N-2 instead of just n=N/2..N-1.
    """
    if n == N:        return "core"        # present in all lines
    if n == N - 1:    return "soft_core"   # present in all but one (annotation-gap-friendly)
    if n >= 2:        return "shell"       # dispensable/accessory, 2..N-2 lines
    return "private"                       # singleton, only one line

cat = n_lines.apply(lambda n: category(n, N))
cat_df = pd.DataFrame({"n_lines": n_lines, "category": cat})
cat_df = cat_df.join(pa)
cat_df.to_csv(os.path.join(RES_DIR, "pan_categories.tsv"), sep="\t")

cat_counts = cat.value_counts().rename_axis("category").reset_index(name="count")
cat_counts.to_csv(os.path.join(RES_DIR, "pan_categories_counts.tsv"),
                  sep="\t", index=False)
print(cat_counts.to_string(index=False))

# ---------- R vs S differential ----------
print(">>> Resistance-vs-Westar differential analysis")
present_in_R = pa[RESISTANT].sum(axis=1)
present_in_S = pa[SUSCEPT].sum(axis=1)

# Orthogroups in all 5 R but absent from Westar
R_only_all   = pa[(present_in_R == len(RESISTANT)) & (present_in_S == 0)].index
# Orthogroups in any R but absent from Westar
R_any_no_S   = pa[(present_in_R >= 1) & (present_in_S == 0)].index
# Westar-only orthogroups (absent from all R)
S_only       = pa[(present_in_R == 0) & (present_in_S == 1)].index
# Shared (any R + Westar)
shared       = pa[(present_in_R >= 1) & (present_in_S == 1)].index

rs = pd.DataFrame({
    "R_only_all5":           [og in set(R_only_all) for og in pa.index],
    "R_any_absent_in_S":     [og in set(R_any_no_S) for og in pa.index],
    "S_only":                [og in set(S_only)     for og in pa.index],
    "shared_R_and_S":        [og in set(shared)     for og in pa.index],
    "n_resistant_lines":     present_in_R.values,
    "westar_present":        present_in_S.values.astype(bool),
}, index=pa.index)
rs.to_csv(os.path.join(RES_DIR, "R_vs_S.tsv"), sep="\t")

print(f"    OGs present in ALL 5 R-lines but absent from Westar: {len(R_only_all)}")
print(f"    OGs present in >=1 R-line but absent from Westar:    {len(R_any_no_S)}")
print(f"    Westar-only OGs:                                     {len(S_only)}")
print(f"    Shared OGs (>=1 R + Westar):                         {len(shared)}")

# ---------- merge NLRtracker class info ----------
print(">>> Merging NLRtracker class annotations")
tracker_paths = {}
for line, rel in DEFAULT_TRACKER.items():
    p = os.path.normpath(os.path.join(HERE, rel))
    if os.path.exists(p):
        tracker_paths[line] = p

if not tracker_paths:
    print("    (no NLRtracker.tsv files found; skipping)")
    OG_class = pd.DataFrame()
else:
    # Build gene_id -> class
    gene2class: dict[str,str] = {}
    for line, path in tracker_paths.items():
        try:
            t = pd.read_csv(path, sep="\t")
        except Exception as e:
            print(f"    WARN: failed to read {path}: {e}")
            continue
        col_class = None
        for c in t.columns:
            if c.lower() in ("subclass (putative)", "subclass", "class"):
                col_class = c; break
        if col_class is None: continue
        for sid, cls in zip(t["seqname"], t[col_class]):
            sid_eff = str(sid)
            # H1-H5 fastas were renamed HI*->IH*; same here so seqnames align
            if line in ("H1","H2","H3","H4","H5"):
                idx = line.replace("H","")
                sid_eff = sid_eff.replace(f"HI{idx}_", f"IH{idx}_")
            gene2class[sid_eff] = cls

    # Load Orthogroups.tsv to get gene -> OG mapping
    og_table = pd.read_csv(OG_FILE, sep="\t")
    rows = []
    for _, r in og_table.iterrows():
        og = r["Orthogroup"]
        classes = []
        for col in og_table.columns[1:]:
            if pd.isna(r[col]): continue
            for g in str(r[col]).split(", "):
                g = g.strip()
                if g in gene2class:
                    classes.append(gene2class[g])
        if classes:
            c = Counter(classes)
            dominant = c.most_common(1)[0][0]
            rows.append((og, dominant, ",".join(f"{k}:{v}" for k,v in c.most_common())))
    OG_class = pd.DataFrame(rows, columns=["Orthogroup","dominant_class","class_breakdown"]).set_index("Orthogroup")
    OG_class.to_csv(os.path.join(RES_DIR, "OG_with_class.tsv"), sep="\t")
    print(f"    annotated {len(OG_class)} OGs with class")

    # plot per-class pan-category breakdown
    if len(OG_class):
        merged = cat_df.join(OG_class[["dominant_class"]], how="left")
        merged["dominant_class"] = merged["dominant_class"].fillna("Unknown")
        ct = pd.crosstab(merged["dominant_class"], merged["category"])
        ct = ct.reindex(columns=[c for c in order if c in ct.columns])
        ct.plot(kind="bar", stacked=True, figsize=(8,4),
                colormap="viridis")
        plt.ylabel("# orthogroups"); plt.title("NLR-class × pan-category")
        plt.xticks(rotation=20, ha="right"); plt.tight_layout()
        plt.savefig(os.path.join(FIG_DIR, "class_vs_pan_category.png"), dpi=200)
        plt.close()

# ---------- rarefaction ----------
print(">>> Rarefaction (permuted pan/core size)")
N_PERM = 200
rng = random.Random(42)
samples = list(pa.columns)
records = []
for perm in range(N_PERM):
    order_perm = samples[:]; rng.shuffle(order_perm)
    pan_set, core_set = set(), None
    for k, line in enumerate(order_perm, 1):
        present = set(pa.index[pa[line] > 0])
        pan_set |= present
        core_set = present if core_set is None else core_set & present
        records.append((perm, k, "pan",  len(pan_set)))
        records.append((perm, k, "core", len(core_set)))
rare = pd.DataFrame(records, columns=["perm","k","kind","size"])
rare_agg = rare.groupby(["k","kind"])["size"].agg(["mean","std","min","max"]).reset_index()
rare_agg.to_csv(os.path.join(RES_DIR, "rarefaction.tsv"), sep="\t", index=False)

plt.figure(figsize=(6,4))
for kind, col in [("pan","#2c7fb8"),("core","#d9534f")]:
    sub = rare_agg[rare_agg["kind"]==kind]
    plt.plot(sub["k"], sub["mean"], "-o", label=kind, color=col)
    plt.fill_between(sub["k"], sub["mean"]-sub["std"], sub["mean"]+sub["std"], color=col, alpha=0.2)
plt.xlabel("# lines"); plt.ylabel("# orthogroups")
plt.title("Pan-NLRome rarefaction (200 permutations)")
plt.legend(); plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, "rarefaction.png"), dpi=200)
plt.close()

# ---------- summary ----------
print()
print(">>> Wrote results to:", RES_DIR)
for f in sorted(os.listdir(RES_DIR)):
    if os.path.isfile(os.path.join(RES_DIR,f)):
        print("   ", f)
print(">>> Figures in:", FIG_DIR)
