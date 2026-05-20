#!/usr/bin/env bash
## Script 49: Download deCODE UKB Olink pQTL files for priority proteins
## Usage: bash scripts/49_decode_download.sh <TOKEN>
## Example: bash scripts/49_decode_download.sh fc9c4647-1f4f-4a51-a67a-5281462a06e9

set -euo pipefail

TOKEN="${1:?Usage: $0 <TOKEN>}"
BASE_URL="https://download.decode.is"
PROJ="/Users/vijayachitramodhukur/Library/Mobile Documents/com~apple~CloudDocs/ECLAI/MultiOmic_Network_MR_Project"
OUT_DIR="${PROJ}/data/decode_pqtl"
mkdir -p "${OUT_DIR}"

## Priority proteins (Tier 1 first, then Tier 2)
PROTEINS=("EFNA1" "ATRAID" "TNFRSF6B" "SNX15" "UMOD" "IL34" "PM20D1" "CGREF1" "ITIH3" "SWAP70" "INHBB" "APOE")

echo "=== Step 1: Listing all files in deCODE proteomics folder ==="
LIST_JSON="${OUT_DIR}/decode_file_list.json"

curl -fsSL "${BASE_URL}/s3/folder?token=${TOKEN}" -o "${LIST_JSON}"
echo "File list saved to: ${LIST_JSON}"
echo "Total files: $(python3 -c "import json; d=json.load(open('${LIST_JSON}')); print(len(d))" 2>/dev/null || echo "run python3 to count")"

echo ""
echo "=== Step 2: Finding European UKB files for priority proteins ==="
MATCH_LIST="${OUT_DIR}/decode_targets.txt"

python3 - <<'PYEOF'
import json, sys, os

json_path = os.path.expanduser(
    "/Users/vijayachitramodhukur/Library/Mobile Documents/com~apple~CloudDocs/ECLAI/MultiOmic_Network_MR_Project/data/decode_pqtl/decode_file_list.json"
)
out_path = os.path.expanduser(
    "/Users/vijayachitramodhukur/Library/Mobile Documents/com~apple~CloudDocs/ECLAI/MultiOmic_Network_MR_Project/data/decode_pqtl/decode_targets.txt"
)

proteins = ["EFNA1","ATRAID","TNFRSF6B","SNX15","UMOD","IL34","PM20D1","CGREF1","ITIH3","SWAP70","INHBB","APOE"]

with open(json_path) as f:
    data = json.load(f)

# Handle both list-of-strings and list-of-dicts
if isinstance(data, list) and len(data) > 0 and isinstance(data[0], dict):
    files = [item.get("name","") or item.get("key","") or item.get("fileName","") for item in data]
elif isinstance(data, list):
    files = data
elif isinstance(data, dict):
    files = data.get("files", data.get("contents", list(data.keys())))
else:
    files = []

print(f"Total files found: {len(files)}")
targets = []
for prot in proteins:
    # Match European-only: contains protein name, starts with GBR_UKB_OLINK (not Africa/SAsia)
    matches = [f for f in files if
               f"_{prot}_" in f and
               f.startswith("GBR_UKB_OLINK") and
               "_Africa_" not in f and
               "_SAsia_" not in f and
               f.endswith(".txt.gz")]
    if matches:
        chosen = matches[0]
        targets.append(chosen)
        print(f"  ✓ {prot:12s} → {chosen}")
    else:
        print(f"  ✗ {prot:12s} → NOT FOUND (check name in portal)")

with open(out_path, "w") as f:
    f.write("\n".join(targets) + "\n")

print(f"\n{len(targets)}/{len(proteins)} proteins matched → saved to decode_targets.txt")
PYEOF

echo ""
echo "=== Step 3: Downloading matched files ==="
echo "(~910 MB each — downloads run sequentially with progress)"
echo ""

while IFS= read -r filename; do
    [[ -z "$filename" ]] && continue
    dest="${OUT_DIR}/${filename}"
    md5dest="${OUT_DIR}/${filename%.gz}.md5sum"

    if [[ -f "$dest" ]]; then
        echo "  [SKIP] Already exists: $filename"
        continue
    fi

    echo "  ↓ Downloading: $filename"
    curl -L --progress-bar \
         "${BASE_URL}/s3/download?token=${TOKEN}&file=${filename}" \
         -o "${dest}"

    # Download md5 checksum
    md5file="${filename%.txt.gz}.txt.md5sum"
    curl -fsSL \
         "${BASE_URL}/s3/download?token=${TOKEN}&file=${md5file}" \
         -o "${OUT_DIR}/${md5file}" 2>/dev/null || true

    echo "  ✓ Done: $filename"
    echo ""
done < "${OUT_DIR}/decode_targets.txt"

echo "=== All downloads complete ==="
echo "Files saved to: ${OUT_DIR}"
ls -lh "${OUT_DIR}"/*.txt.gz 2>/dev/null || echo "(no .gz files yet)"
