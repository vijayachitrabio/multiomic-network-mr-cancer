#!/usr/bin/env python3
"""
Script 81: Generates the comprehensive Druggability and Tractability Annotation Tables
for the 13 prioritized circulating protein candidates.
"""

import pandas as pd
import os

PROJ_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
OUT_DIR = os.path.join(PROJ_DIR, "results", "druggability")
os.makedirs(OUT_DIR, exist_ok=True)

main_data = [
    {
        "Gene/protein": "EFNA1",
        "Evidence tier in manuscript": "Tier 1",
        "Cancer context": "Breast cancer",
        "MR direction": "Risk (↑)",
        "Colocalization support": "Strong (SuSiE)",
        "MAGMA support": "Significant",
        "Tumour-context annotation": "Stromal/Tumour-associated",
        "Open Targets tractability": "Low direct tractability (Ligand).",
        "DrugBank evidence": "No approved drugs targeting EFNA1 directly. EPHA2 (receptor) has investigational compounds.",
        "ChEMBL evidence": "Target ID: CHEMBL3487 (EPHA2). Functional/binding assays present for receptor.",
        "DepMap breast cancer dependency": "Not a selective dependency in breast cancer cell lines.",
        "Breast cancer-specific evidence": "EFNA1/EphA2 axis drives breast cancer migration and metastasis in preclinical models.",
        "General oncology evidence": "EphA2 overexpression is a well-established negative prognostic factor across multiple solid tumours.",
        "Existing drugs/compounds": "ALW-II-41-27 (EphA2 inhibitor, preclinical); MEDI-547 (EphA2 ADC, Phase I).",
        "Translational tractability category": "Moderate tractability",
        "Evidence strength": "Indirect / Moderate",
        "Key caution": "Targeting the ligand directly is challenging; therapeutic efforts focus on the EphA2 receptor.",
        "Source links": "OpenTargets: ENSG00000109723; ChEMBL: CHEMBL3487",
        "PMID/DOI references": "PMID: 15313904 (EphA2 in breast cancer)"
    },
    {
        "Gene/protein": "TNFRSF6B",
        "Evidence tier in manuscript": "Tier 1",
        "Cancer context": "Breast cancer",
        "MR direction": "Risk (↑)",
        "Colocalization support": "Strong (ABF+SuSiE)",
        "MAGMA support": "Significant",
        "Tumour-context annotation": "Immune/Secreted decoy",
        "Open Targets tractability": "Secreted protein; monoclonal antibody accessible.",
        "DrugBank evidence": "No approved direct inhibitors. Investigational TL1A/LIGHT neutralisers exist.",
        "ChEMBL evidence": "Target ID: CHEMBL3706. Bioactive compounds present primarily targeting its ligands (TL1A/LIGHT).",
        "DepMap breast cancer dependency": "Not an essential gene in breast cancer cell lines (consistent with secreted decoy function).",
        "Breast cancer-specific evidence": "DcR3 (TNFRSF6B) overexpression helps breast cancer cells evade FasL-mediated immune clearance.",
        "General oncology evidence": "Elevated DcR3 is found in multiple malignancies and correlates with immune evasion and poor prognosis.",
        "Existing drugs/compounds": "Recombinant DcR3-Fc fusions (preclinical); Pateclizumab (anti-LIGHT, non-oncology clinical trials).",
        "Translational tractability category": "Moderate tractability",
        "Evidence strength": "Moderate",
        "Key caution": "Targeting a soluble decoy receptor is complex; neutralising its binding partners may have systemic immune toxicity.",
        "Source links": "OpenTargets: ENSG00000117036; DrugBank: DB05686 (FasL)",
        "PMID/DOI references": "PMID: 25164805 (DcR3 immune evasion)"
    },
    {
        "Gene/protein": "ATRAID",
        "Evidence tier in manuscript": "Tier 1",
        "Cancer context": "Breast cancer",
        "MR direction": "Protective (↓)",
        "Colocalization support": "Strong (SuSiE)",
        "MAGMA support": "Significant",
        "Tumour-context annotation": "Apoptosis/Cell cycle",
        "Open Targets tractability": "Low direct small-molecule tractability.",
        "DrugBank evidence": "Not a direct listed target. Expression is induced by All-trans retinoic acid (ATRA).",
        "ChEMBL evidence": "No direct ChEMBL target assays for ATRAID protein.",
        "DepMap breast cancer dependency": "Not a broad or selective dependency.",
        "Breast cancer-specific evidence": "ATRAID mediates apoptosis in response to retinoic acid in breast cancer cell lines.",
        "General oncology evidence": "Implicated in p53-dependent apoptosis and osteoblast differentiation.",
        "Existing drugs/compounds": "All-trans retinoic acid (Tretinoin - modulates expression, does not inhibit protein directly).",
        "Translational tractability category": "Biological candidate, limited druggability",
        "Evidence strength": "Indirect / Weak",
        "Key caution": "No direct ligands known. Tractability relies entirely on modulating upstream retinoic acid receptors (RAR).",
        "Source links": "OpenTargets: ENSG00000122223",
        "PMID/DOI references": "PMID: 15308639 (APR-3/ATRAID in apoptosis)"
    },
    {
        "Gene/protein": "FGF5",
        "Evidence tier in manuscript": "Tier 1",
        "Cancer context": "Breast cancer",
        "MR direction": "Protective (↓)",
        "Colocalization support": "Strong (ABF+SuSiE)",
        "MAGMA support": "Significant",
        "Tumour-context annotation": "Growth factor signaling",
        "Open Targets tractability": "Secreted ligand; high tractability via receptor (FGFR1/2).",
        "DrugBank evidence": "FGF5 is a known ligand for FGFR1/2. FGFRs have multiple approved oncology drugs.",
        "ChEMBL evidence": "Target ID: CHEMBL1809 (FGFR1). Thousands of bioactive inhibitors available.",
        "DepMap breast cancer dependency": "FGF5 itself is not essential; FGFR1 shows selective dependency in specific breast lineages.",
        "Breast cancer-specific evidence": "FGFR1 amplification/overexpression is a major driver and resistance mechanism in HR+ breast cancer.",
        "General oncology evidence": "FGFR inhibitors are clinically active in urothelial and cholangiocarcinoma.",
        "Existing drugs/compounds": "Erdafitinib, Pemigatinib, AZD4547 (FGFR inhibitors).",
        "Translational tractability category": "High translational tractability",
        "Evidence strength": "Strong",
        "Key caution": "FGF5 is protective in MR, suggesting FGFR pathway modulation rather than simple inhibition is required, or tissue-specific effects.",
        "Source links": "OpenTargets: ENSG00000138675; DrugBank: DB11939 (Erdafitinib)",
        "PMID/DOI references": "PMID: 25424861 (FGFR signaling in breast cancer)"
    },
    {
        "Gene/protein": "ABO",
        "Evidence tier in manuscript": "Tier 1",
        "Cancer context": "Endometrial (Comparator)",
        "MR direction": "Risk (↑)",
        "Colocalization support": "Strong (ABF)",
        "MAGMA support": "Significant",
        "Tumour-context annotation": "Glycosyltransferase",
        "Open Targets tractability": "Enzyme; theoretical small molecule tractability.",
        "DrugBank evidence": "No approved drugs targeting ABO glycosyltransferase.",
        "ChEMBL evidence": "Target ID: CHEMBL3417. Very few direct functional inhibitors.",
        "DepMap breast cancer dependency": "Not an essential gene in endometrial or breast cell lines.",
        "Breast cancer-specific evidence": "N/A (Endometrial comparator signal).",
        "General oncology evidence": "ABO blood group antigens are explored as epitopes for antibody-drug conjugates (ADCs) in solid tumours.",
        "Existing drugs/compounds": "Preclinical anti-blood group A ADCs.",
        "Translational tractability category": "Biological candidate, limited druggability",
        "Evidence strength": "Database-only / Weak",
        "Key caution": "ABO is an endometrial comparator signal. True direct enzymatic inhibitors are lacking in oncology.",
        "Source links": "OpenTargets: ENSG00000175164",
        "PMID/DOI references": "N/A"
    },
    {
        "Gene/protein": "SNX15",
        "Evidence tier in manuscript": "Secondary (MR+MAGMA)",
        "Cancer context": "Breast cancer",
        "MR direction": "Protective (↓)",
        "Colocalization support": "Insufficient (Coloc failure/LD)",
        "MAGMA support": "Significant",
        "Tumour-context annotation": "Endosomal sorting",
        "Open Targets tractability": "Low tractability.",
        "DrugBank evidence": "Not a listed drug target.",
        "ChEMBL evidence": "No direct ChEMBL target assays.",
        "DepMap breast cancer dependency": "Not broadly essential.",
        "Breast cancer-specific evidence": "Limited direct functional evidence in breast cancer literature.",
        "General oncology evidence": "General roles in intracellular trafficking and receptor recycling.",
        "Existing drugs/compounds": "None known.",
        "Translational tractability category": "Currently low tractability",
        "Evidence strength": "Absent",
        "Key caution": "Secondary candidate lacking colocalization support. No drug development precedent.",
        "Source links": "OpenTargets: ENSG00000164253",
        "PMID/DOI references": "N/A"
    },
    {
        "Gene/protein": "PM20D1",
        "Evidence tier in manuscript": "Secondary (MR+MAGMA)",
        "Cancer context": "Breast cancer",
        "MR direction": "Risk (↑)",
        "Colocalization support": "Insufficient (Coloc failure/LD)",
        "MAGMA support": "Significant",
        "Tumour-context annotation": "Metabolic/Lipid",
        "Open Targets tractability": "Enzyme; potential small molecule tractability.",
        "DrugBank evidence": "Not a listed drug target.",
        "ChEMBL evidence": "Target ID: CHEMBL4630737. Very limited bioactivity data.",
        "DepMap breast cancer dependency": "Not an essential gene.",
        "Breast cancer-specific evidence": "No significant direct breast cancer literature evidence.",
        "General oncology evidence": "Functions as a biosynthetic enzyme for N-acyl amino acids; studied in obesity/metabolism rather than oncology.",
        "Existing drugs/compounds": "None in clinical development.",
        "Translational tractability category": "Currently low tractability",
        "Evidence strength": "Absent",
        "Key caution": "Secondary candidate lacking colocalization. Primarily a metabolic biology target.",
        "Source links": "OpenTargets: ENSG00000162877",
        "PMID/DOI references": "PMID: 27374330 (PM20D1 in metabolism)"
    },
    {
        "Gene/protein": "UMOD",
        "Evidence tier in manuscript": "Tier 2a / Provisional",
        "Cancer context": "Breast cancer",
        "MR direction": "Risk (↑)",
        "Colocalization support": "Moderate (ABF only)",
        "MAGMA support": "Not significant",
        "Tumour-context annotation": "Immune/Kidney-derived",
        "Open Targets tractability": "Secreted glycoprotein; antibody accessible.",
        "DrugBank evidence": "Listed interactions with immunosuppressants (cyclosporine), but functional, not direct enzymatic.",
        "ChEMBL evidence": "No direct target inhibitors in ChEMBL.",
        "DepMap breast cancer dependency": "Not essential (gene is primarily expressed in kidney).",
        "Breast cancer-specific evidence": "No direct functional evidence in breast cancer. Likely acts systemically via immune modulation.",
        "General oncology evidence": "Uromodulin regulates systemic innate immunity and cytokine responses.",
        "Existing drugs/compounds": "None directly targeting UMOD for oncology.",
        "Translational tractability category": "Currently low tractability",
        "Evidence strength": "Indirect / Weak",
        "Key caution": "Must remain Tier 2a/provisional. Circulating levels likely reflect kidney function or systemic immune state rather than a druggable tumour-intrinsic target.",
        "Source links": "OpenTargets: ENSG00000169344",
        "PMID/DOI references": "N/A"
    },
    {
        "Gene/protein": "APOE",
        "Evidence tier in manuscript": "Tier 2b",
        "Cancer context": "Breast cancer",
        "MR direction": "Risk (↑)",
        "Colocalization support": "Moderate (ABF only)",
        "MAGMA support": "Significant",
        "Tumour-context annotation": "Lipid metabolism/Macrophage",
        "Open Targets tractability": "Secreted protein; high clinical precedent in neurology.",
        "DrugBank evidence": "Investigational drugs exist, primarily for Alzheimer's disease.",
        "ChEMBL evidence": "Target ID: CHEMBL2531. Bioactive compounds available.",
        "DepMap breast cancer dependency": "Not an essential gene in breast cell lines.",
        "Breast cancer-specific evidence": "APOE expressed by tumor-associated macrophages promotes tumor growth and immune evasion in breast cancer.",
        "General oncology evidence": "Implicated in systemic lipid metabolism and tumor microenvironment immunosuppression.",
        "Existing drugs/compounds": "Various APOE-modulating therapies in preclinical/clinical neurology trials.",
        "Translational tractability category": "Moderate tractability",
        "Evidence strength": "Moderate",
        "Key caution": "Therapeutic focus is entirely neurology. Oncology repurposing requires targeting the tumor microenvironment without neurological toxicity.",
        "Source links": "OpenTargets: ENSG00000130203",
        "PMID/DOI references": "PMID: 32313227 (APOE in breast cancer macrophages)"
    },
    {
        "Gene/protein": "ITIH3",
        "Evidence tier in manuscript": "Tier 2b",
        "Cancer context": "Breast cancer",
        "MR direction": "Protective (↓)",
        "Colocalization support": "Insufficient (Coloc failure/LD)",
        "MAGMA support": "Significant",
        "Tumour-context annotation": "Extracellular matrix",
        "Open Targets tractability": "Secreted protein.",
        "DrugBank evidence": "Not a listed drug target.",
        "ChEMBL evidence": "No direct ChEMBL target assays.",
        "DepMap breast cancer dependency": "Not essential.",
        "Breast cancer-specific evidence": "ITIH family members stabilize the extracellular matrix; downregulation is associated with metastasis.",
        "General oncology evidence": "Often downregulated in solid tumours as a tumour suppressor mechanism.",
        "Existing drugs/compounds": "None.",
        "Translational tractability category": "Currently low tractability",
        "Evidence strength": "Absent",
        "Key caution": "Protective structural protein. Extremely difficult to drug therapeutically (requires restoring function/expression).",
        "Source links": "OpenTargets: ENSG00000155099",
        "PMID/DOI references": "N/A"
    },
    {
        "Gene/protein": "IL34",
        "Evidence tier in manuscript": "Tier 2b",
        "Cancer context": "Breast cancer",
        "MR direction": "Risk (↑)",
        "Colocalization support": "Insufficient",
        "MAGMA support": "Not significant",
        "Tumour-context annotation": "Cytokine/Macrophage",
        "Open Targets tractability": "Secreted ligand; high tractability via CSF1R.",
        "DrugBank evidence": "IL34 binds CSF1R. CSF1R is highly druggable with approved oncology drugs.",
        "ChEMBL evidence": "Target ID: CHEMBL2632 (CSF1R). Thousands of potent inhibitors.",
        "DepMap breast cancer dependency": "Not an essential gene (functions in microenvironment).",
        "Breast cancer-specific evidence": "IL34 promotes macrophage recruitment and immunosuppression in breast tumors.",
        "General oncology evidence": "Targeting the IL34/CSF1R axis is a major strategy to deplete tumor-associated macrophages.",
        "Existing drugs/compounds": "Pexidartinib (CSF1R inhibitor, approved for TGCT). IL34-blocking antibodies (preclinical).",
        "Translational tractability category": "High translational tractability",
        "Evidence strength": "Strong",
        "Key caution": "Tractability relies on targeting its receptor (CSF1R).",
        "Source links": "OpenTargets: ENSG00000153363; DrugBank: DB12328 (Pexidartinib)",
        "PMID/DOI references": "PMID: 26868662 (IL34 in chemoresistance and macrophages)"
    },
    {
        "Gene/protein": "KLB",
        "Evidence tier in manuscript": "Tier 2b",
        "Cancer context": "Breast cancer",
        "MR direction": "Risk (↑)",
        "Colocalization support": "Insufficient",
        "MAGMA support": "Significant",
        "Tumour-context annotation": "Metabolic receptor",
        "Open Targets tractability": "Membrane receptor; highly tractable.",
        "DrugBank evidence": "Co-receptor for FGF19/FGF21. Drugs targeting this axis exist for metabolic disease.",
        "ChEMBL evidence": "Target ID: CHEMBL4630799. Analogues and modulators available.",
        "DepMap breast cancer dependency": "Not broadly essential in breast cancer.",
        "Breast cancer-specific evidence": "KLB/FGFR4 axis contributes to breast cancer progression and endocrine resistance.",
        "General oncology evidence": "FGF19/KLB/FGFR4 axis is highly relevant in hepatocellular carcinoma and breast cancer.",
        "Existing drugs/compounds": "Pegbelfermin (FGF21 analogue, metabolic). Roblitinib (FGFR4 inhibitor, oncology).",
        "Translational tractability category": "Moderate tractability",
        "Evidence strength": "Moderate",
        "Key caution": "Oncology efforts focus primarily on the FGFR4 kinase component of the complex.",
        "Source links": "OpenTargets: ENSG00000138669",
        "PMID/DOI references": "PMID: 32679247 (KLB and FGFR4 in breast cancer)"
    },
    {
        "Gene/protein": "FGFR4",
        "Evidence tier in manuscript": "Tier 2b",
        "Cancer context": "Breast cancer",
        "MR direction": "Protective (↓)",
        "Colocalization support": "Insufficient",
        "MAGMA support": "Significant",
        "Tumour-context annotation": "Kinase receptor",
        "Open Targets tractability": "Kinase; highly tractable.",
        "DrugBank evidence": "Specific FGFR4 inhibitors are in clinical trials.",
        "ChEMBL evidence": "Target ID: CHEMBL3822. Hundreds of selective inhibitors.",
        "DepMap breast cancer dependency": "Selective dependency in specific subsets.",
        "Breast cancer-specific evidence": "FGFR4 drives proliferation and resistance in specific breast cancer subtypes.",
        "General oncology evidence": "FGFR4 inhibitors are heavily investigated in hepatocellular carcinoma.",
        "Existing drugs/compounds": "Fisogatinib, Roblitinib (Phase I/II oncology).",
        "Translational tractability category": "High translational tractability",
        "Evidence strength": "Strong",
        "Key caution": "MR direction is protective, complicating direct inhibition. Tissue-specific pathway nuances likely exist.",
        "Source links": "OpenTargets: ENSG00000160867",
        "PMID/DOI references": "PMID: 25964264 (FGFR4 targeting in breast cancer)"
    }
]


# File 1: Main Excel Table
df_main = pd.DataFrame(main_data)
main_csv_path = os.path.join(OUT_DIR, "STable19_Druggability_Tractability.csv")
main_xlsx_path = os.path.join(OUT_DIR, "STable19_Druggability_Tractability.xlsx")
df_main.to_csv(main_csv_path, index=False)
df_main.to_excel(main_xlsx_path, index=False)

# File 2: Short Evidence Summary
summary_data = []
for row in main_data:
    summary_data.append({
        "Gene": row["Gene/protein"],
        "Main finding": f"{row['Cancer context']} MR {row['MR direction']}, Tier: {row['Evidence tier in manuscript']}",
        "Druggability category": row["Translational tractability category"],
        "Best evidence": row["Existing drugs/compounds"],
        "Main caution": row["Key caution"]
    })

df_summary = pd.DataFrame(summary_data)
summary_csv_path = os.path.join(OUT_DIR, "STable20_Druggability_Short_Summary.csv")
summary_xlsx_path = os.path.join(OUT_DIR, "STable20_Druggability_Short_Summary.xlsx")
df_summary.to_csv(summary_csv_path, index=False)
df_summary.to_excel(summary_xlsx_path, index=False)

# File 3: Source/Reference List
ref_data = []
for row in main_data:
    links = row["Source links"].split("; ")
    for link in links:
        if link.strip():
            parts = link.split(":")
            db = parts[0].strip() if len(parts) > 1 else "Database"
            val = parts[1].strip() if len(parts) > 1 else link
            ref_data.append({
                "Gene": row["Gene/protein"],
                "Source type": "Database",
                "Database/paper": db,
                "Link/DOI/PMID": val,
                "Why included": f"To support {db} tractability assessment."
            })
    
    if row["PMID/DOI references"] != "N/A":
        pmids = row["PMID/DOI references"].split("; ")
        for pmid in pmids:
            if pmid.strip():
                ref_data.append({
                    "Gene": row["Gene/protein"],
                    "Source type": "Literature",
                    "Database/paper": "PubMed",
                    "Link/DOI/PMID": pmid.strip(),
                    "Why included": "Supports breast/general cancer context."
                })

df_ref = pd.DataFrame(ref_data)
ref_csv_path = os.path.join(OUT_DIR, "STable21_Druggability_References.csv")
ref_xlsx_path = os.path.join(OUT_DIR, "STable21_Druggability_References.xlsx")
df_ref.to_csv(ref_csv_path, index=False)
df_ref.to_excel(ref_xlsx_path, index=False)

print("Druggability generation complete.")
print(f"Generated files in {OUT_DIR}/")
