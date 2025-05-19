#!/bin/bash
# run like: ./embed_icc.sh /path/to/images ../data/icc/ProPhoto.icc

INPUT_DIR="$1"
ICC_PATH="$2"
#ICC_PATH="../data/icc/ProPhoto.icc"

cd "$INPUT_DIR" || exit 1

for f in *.png; do
  magick "$f" -profile "$ICC_PATH" "icc_$f"
done
