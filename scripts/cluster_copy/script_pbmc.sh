#!/bin/bash
# Copy filtered_feature_bc_matrix.h5 for selected samples into BIO559R/<sample>/

SRC="$HOME/groups/fslg_dnasc/nobackup/archive/RIKEN_OUTPUT/CELLRANGER_OUTPUT"
DEST="$SRC/BIO559R"    # change this if you want BIO559R somewhere else

for s in 10x_A116 10x_A117 10x_A118 10x_A170 10x_A171 10x_A172; do
    f="$SRC/${s}_cellranger/outs/filtered_feature_bc_matrix.h5"
    if [[ -f "$f" ]]; then
        mkdir -p "$DEST/$s"
        cp -v "$f" "$DEST/$s/"
    else
        echo "WARNING: missing $f" >&2
    fi
done
