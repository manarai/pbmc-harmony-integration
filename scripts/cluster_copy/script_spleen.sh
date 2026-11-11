#!/bin/bash
# Copy filtered_feature_bc_matrix.h5 for selected samples into BIO559R/<sample>/

SRC="$HOME/groups/fslg_dnasc/nobackup/archive/RIKEN_OUTPUT/CELLRANGER_OUTPUT"
DEST="$SRC/BIO559R"    # change this if you want BIO559R somewhere else

for s in A024 A025 A026 A027 A028 A029; do
    f="$SRC/${s}_cellranger/outs/filtered_feature_bc_matrix.h5"
    if [[ -f "$f" ]]; then
        mkdir -p "$DEST/$s"
        cp -v "$f" "$DEST/$s/"
    else
        echo "WARNING: missing $f" >&2
    fi
done
