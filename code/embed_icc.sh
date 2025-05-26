#!/usr/bin/env bash
# run like: ./embed_icc.sh /path/to/images icc/ProPhoto.icc

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
INPUT_DIR="$1"
ICC_PATH="$SCRIPT_DIR/$2"

# derive suffix from profile name, e.g. "ProPhoto_D50.icc" → "D50"
profile_file="$(basename "$ICC_PATH")"          # ProPhoto_D50.icc
profile_base="${profile_file%.*}"               # ProPhoto_D50
suffix="${profile_base#*_}"                     # D50

cd "$INPUT_DIR"

# 1) make the output subfolder
mkdir -p prophoto

for f in *.png; do
  # skip files we've already created
  if [[ "$f" == icc_*.png ]]; then
    continue
  fi

  # strip extension from original file
  base="${f%.*}"

  # 2) build output filename inside `prophoto/`
  out="prophoto/icc_${base}_${suffix}.png"

  magick "$f" -profile "$ICC_PATH" "$out"
done
