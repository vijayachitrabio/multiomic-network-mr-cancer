#!/usr/bin/env python3
"""
Script 63: Generate complete ARIC + OpenGWAS replication table
Covers all 17 FDR-significant proteins
Status codes: ✓ replicated / — not available / ✗ failed
Date: 2026-05-25
"""
import csv, math, os

PROJ = "/Users/vijayachitramodhukur/Library/Mobile Documents/com~apple~CloudDocs/ECLAI/MultiOmic_Network_MR_Project"

# ── Load discovery results (all 17 proteins) ────────────────────────────────
discovery = {}
with open(f"{PROJ}/results/tables/STable2_17_FDR_hits_complete.csv") as f:
    for r in csv.DictReader(f):
        discovery[(r['protein'], r['cancer'])] = r

# ── Load ARIC EA-cis results (primary ARIC analysis) ────────────────────────
aric = {}
with open(f"{PROJ}/results/replication/aric_replication_mr_results_combined.csv") as f:
    for r in csv.DictReader(f):
        if r['aric_analysis'] == 'ea_cis':
            key = (r['exposure'], r['outcome'])
            if key not in aric:
                aric[key] = r

# ── Load ARIC coverage (to know why proteins are missing) ───────────────────
aric_cov = {}
with open(f"{PROJ}/results/replication/aric_replication_coverage_combined.csv") as f:
    for r in csv.DictReader(f):
        if r['aric_analysis'] == 'ea_cis':
            aric_cov[(r['protein'], r['outcome'])] = r['status']

# ── Load OpenGWAS results ────────────────────────────────────────────────────
# Primary method per protein: IVW if multi-SNP, Wald ratio if single
opengwas = {}
with open(f"{PROJ}/results/opengwas/opengwas_5protein_replication_mr_results.csv") as f:
    for r in csv.DictReader(f):
        prot = r['protein']
        outcome = r['outcome']
        method = r['method']
        # Keep IVW preferentially
        key = (prot, outcome)
        if key not in opengwas or method == 'Inverse variance weighted':
            opengwas[key] = r

# ── All 17 proteins with their discovery cancer ──────────────────────────────
proteins_17 = [
    # Tier 1
    ("EFNA1",    "Breast_GCST90018757",   "T1"),
    ("TNFRSF6B", "Breast_GCST90018757",   "T1"),
    ("ATRAID",   "Breast_GCST90018757",   "T1"),
    ("FGF5",     "Breast_GCST90018757",   "T1"),
    ("UMOD",     "Breast_GCST90018757",   "T1"),
    ("ABO",      "Endometrial_GCST006464","T1"),
    # Tier 2a
    ("SNX15",    "Breast_GCST90018757",   "T2a"),
    ("PM20D1",   "Breast_GCST90018757",   "T2a"),
    # Tier 2b
    ("TSPAN8",   "Breast_GCST90018757",   "T2b"),
    ("APOE",     "Breast_GCST90018757",   "T2b"),
    # Tier 2c
    ("TSPAN8",   "Breast_GCST90018757",   "T2c"),  # placeholder — confirm with manuscript
    ("IL34",     "Breast_GCST90018757",   "T2c"),
    ("APOB",     "Breast_GCST90018757",   "T2c"),
    ("FGFR4",    "Breast_GCST90018757",   "T2c"),
    ("ITIH3",    "Breast_GCST90018757",   "T2c"),
    ("SNX13",    "Breast_GCST90018757",   "T2c"),
    ("TNFRSF11A","Breast_GCST90018757",   "T2c"),
    # Additional proteins with ARIC hits
    ("KLB",      "Breast_GCST90018757",   "T2c"),
    ("INHBB",    "Breast_GCST90018757",   "T2c"),
    ("SWAP70",   "Breast_GCST90018757",   "T2c"),
]

# Deduplicate and use discovery data to get the actual 17
actual_17 = []
seen = set()
for r in discovery.values():
    key = (r['protein'], r['cancer'])
    if key not in seen:
        seen.add(key)
        actual_17.append(r)

def fmt_or(b, se):
    """Return OR [95%CI] string"""
    try:
        b, se = float(b), float(se)
        OR = math.exp(b)
        lo = math.exp(b - 1.96*se)
        hi = math.exp(b + 1.96*se)
        return f"{OR:.3f} [{lo:.3f}-{hi:.3f}]"
    except:
        return "—"

def fmt_p(p):
    try:
        p = float(p)
        if p < 0.001:
            return f"{p:.2e}"
        return f"{p:.4f}"
    except:
        return "—"

def direction(b_disc, b_rep):
    """↑↑ concordant, ↑↓ discordant"""
    try:
        return "↑" if float(b_rep) > 0 else "↓"
    except:
        return "—"

def disc_direction(b):
    try:
        return "↑" if float(b) > 0 else "↓"
    except:
        return "—"

def status_icon(p, b_disc, b_rep, threshold=0.05):
    """✓ / ✗ / — """
    try:
        p = float(p)
        concordant = sign(float(b_disc)) == sign(float(b_rep))
        if p < threshold and concordant:
            return "✓"
        elif p < threshold and not concordant:
            return "✗ (opposite direction)"
        else:
            return "✗ (p≥0.05)"
    except:
        return "—"

def sign(x): return 1 if x >= 0 else -1

# ── Build output rows ─────────────────────────────────────────────────────────
rows = []

for disc in actual_17:
    prot   = disc['protein']
    cancer = disc['cancer']
    tier   = "T1" if disc.get('protein','') in ["EFNA1","TNFRSF6B","ATRAID","FGF5","UMOD","ABO"] else "T2"
    b_disc = disc.get('beta', disc.get('b', ''))
    se_disc = disc.get('se', '')
    p_disc = disc.get('pvalue', '')

    # ── ARIC ──────────────────────────────────────────────────────────────────
    aric_key = (prot, cancer)
    aric_row = aric.get(aric_key)
    aric_cov_status = aric_cov.get(aric_key, 'unknown')

    if aric_row:
        aric_or_ci = fmt_or(aric_row['b'], aric_row['se'])
        aric_p     = fmt_p(aric_row['pval'])
        aric_dir   = direction(b_disc, aric_row['b'])
        aric_nsnp  = aric_row.get('nsnp', '?')
        aric_method = aric_row.get('method', '?')
        aric_status = status_icon(aric_row['pval'], b_disc, aric_row['b'])
    elif aric_cov_status == 'no_aric_cis_instrument':
        aric_or_ci = "—"
        aric_p     = "—"
        aric_dir   = "—"
        aric_nsnp  = "0"
        aric_method = "—"
        aric_status = "— (no instrument)"
    else:
        aric_or_ci = "—"
        aric_p     = "—"
        aric_dir   = "—"
        aric_nsnp  = "—"
        aric_method = "—"
        aric_status = f"— ({aric_cov_status})"

    # ── OpenGWAS (INTERVAL SomaScan) ─────────────────────────────────────────
    # Map cancer label
    cancer_short = "Breast" if "Breast" in cancer else ("Endometrial" if "Endometrial" in cancer else "Ovarian")
    og_key = (prot, cancer_short)
    og_row = opengwas.get(og_key)

    if og_row:
        og_or_ci  = fmt_or(og_row['b'], og_row['se'])
        og_p      = fmt_p(og_row['pval'])
        og_dir    = direction(b_disc, og_row['b'])
        og_nsnp   = og_row.get('nsnp', '?')
        og_method = og_row.get('method', '?')
        og_status = status_icon(og_row['pval'], b_disc, og_row['b'])
    else:
        og_or_ci  = "—"
        og_p      = "—"
        og_dir    = "—"
        og_nsnp   = "0"
        og_method = "—"
        og_status = "— (no instrument)"

    rows.append({
        'protein':            prot,
        'cancer':             cancer_short,
        'tier':               tier,
        'discovery_OR_95CI':  fmt_or(b_disc, se_disc),
        'discovery_p':        fmt_p(p_disc),
        'discovery_dir':      disc_direction(b_disc),
        # ARIC
        'ARIC_nSNP':          aric_nsnp,
        'ARIC_method':        aric_method,
        'ARIC_OR_95CI':       aric_or_ci,
        'ARIC_p':             aric_p,
        'ARIC_direction':     aric_dir,
        'ARIC_status':        aric_status,
        # OpenGWAS
        'OpenGWAS_nSNP':      og_nsnp,
        'OpenGWAS_method':    og_method,
        'OpenGWAS_OR_95CI':   og_or_ci,
        'OpenGWAS_p':         og_p,
        'OpenGWAS_direction': og_dir,
        'OpenGWAS_status':    og_status,
    })

# ── Save CSV ──────────────────────────────────────────────────────────────────
out_path = f"{PROJ}/results/tables/STable_Replication_ARIC_OpenGWAS_2026-05-25.csv"
fields = list(rows[0].keys())
with open(out_path, 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=fields)
    w.writeheader()
    w.writerows(rows)
print(f"Saved: {out_path}\n")

# ── Print formatted table ─────────────────────────────────────────────────────
print(f"{'Protein':<12} {'Cancer':<12} {'Tier':<4}  {'Discovery OR [95%CI]':<22} {'p':>8}  │  "
      f"{'ARIC nSNP':<9} {'ARIC OR [95%CI]':<22} {'p_ARIC':>9}  {'ARIC':>20}  │  "
      f"{'OG nSNP':<8} {'OpenGWAS OR [95%CI]':<22} {'p_OG':>9}  {'OpenGWAS':>25}")
print("─"*185)

for r in rows:
    print(f"{r['protein']:<12} {r['cancer']:<12} {r['tier']:<4}  "
          f"{r['discovery_OR_95CI']:<22} {r['discovery_p']:>8}  │  "
          f"{r['ARIC_nSNP']:<9} {r['ARIC_OR_95CI']:<22} {r['ARIC_p']:>9}  {r['ARIC_status']:>20}  │  "
          f"{r['OpenGWAS_nSNP']:<8} {r['OpenGWAS_OR_95CI']:<22} {r['OpenGWAS_p']:>9}  {r['OpenGWAS_status']:>25}")

# ── Summary counts ────────────────────────────────────────────────────────────
aric_rep   = sum(1 for r in rows if r['ARIC_status'].startswith('✓'))
aric_avail = sum(1 for r in rows if r['ARIC_status'] != '— (no instrument)')
og_rep     = sum(1 for r in rows if r['OpenGWAS_status'].startswith('✓'))
og_avail   = sum(1 for r in rows if r['OpenGWAS_status'] != '— (no instrument)')

print(f"\nSummary:")
print(f"  ARIC:     {aric_rep}/{aric_avail} replicated (of {len(rows)} proteins with instruments available)")
print(f"  OpenGWAS: {og_rep}/{og_avail} replicated (of {len(rows)} proteins with instruments available)")
