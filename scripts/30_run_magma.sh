#!/usr/bin/env bash

# Script 30: Run MAGMA gene-level analysis for breast and endometrial cancer

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
INPUT_DIR="$PROJECT_DIR/results/pathway/magma/inputs"
OUT_DIR="$PROJECT_DIR/results/pathway/magma"
ANNOT_DIR="$OUT_DIR/annot"
RESULTS_DIR="$OUT_DIR/results"

MAGMA_BIN="$PROJECT_DIR/../uterine_fibroids/magma_inputs/magma"
GENE_LOC="$PROJECT_DIR/../uterine_fibroids/magma_inputs/NCBI37.3.gene.loc"
LD_REF="$PROJECT_DIR/../uterine_fibroids/mr_analysis_2026_04_16/reference_data/g1000_eur"

mkdir -p "$ANNOT_DIR" "$RESULTS_DIR"

check_file() {
  if [[ ! -f "$1" ]]; then
    echo "ERROR: missing required file: $1" >&2
    exit 1
  fi
}

check_file "$MAGMA_BIN"
check_file "$GENE_LOC"
check_file "$LD_REF.bed"
check_file "$LD_REF.bim"
check_file "$LD_REF.fam"
check_file "$INPUT_DIR/breast.snploc"
check_file "$INPUT_DIR/breast.pval"
check_file "$INPUT_DIR/endometrial.snploc"
check_file "$INPUT_DIR/endometrial.pval"

chmod +x "$MAGMA_BIN"

echo "============================================================"
echo "  MAGMA gene-level analysis — breast & endometrial cancer"
echo "  $(date)"
echo "============================================================"

for trait in breast endometrial; do
  echo ""
  echo "Annotating $trait ..."
  "$MAGMA_BIN" \
    --annotate \
    --snp-loc "$INPUT_DIR/${trait}.snploc" \
    --gene-loc "$GENE_LOC" \
    --out "$ANNOT_DIR/${trait}"
done

for trait in breast endometrial; do
  N="$(awk 'NR==2{print $3}' "$INPUT_DIR/${trait}.pval")"
  echo ""
  echo "Running gene analysis for $trait (N=$N) ..."
  "$MAGMA_BIN" \
    --bfile "$LD_REF" \
    --pval "$INPUT_DIR/${trait}.pval" use=SNP,P ncol=N \
    --gene-annot "$ANNOT_DIR/${trait}.genes.annot" \
    --out "$RESULTS_DIR/${trait}_genes"
done

echo ""
echo "Summarising MAGMA outputs ..."
Rscript "$SCRIPT_DIR/31_summarise_magma_results.R"

echo ""
echo "MAGMA complete. Results: $RESULTS_DIR"
