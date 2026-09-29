#!/bin/bash
set -euo pipefail
: "${QSTARS_CMFGEN_ROOT:?Absolute extracted CMFGEN source directory required}"
# Use sequential archive updates: the upstream makefiles share libraries.
export MACOSX_DEPLOYMENT_TARGET=14.0
make -C "$QSTARS_CMFGEN_ROOT" d_tools d_blas d_dust d_lpack d_nm_subs d_stark d_unix d_subs d_plane d_newsubs d_pgplt
make -C "$QSTARS_CMFGEN_ROOT/lte_hydro" OBJ_FILES
make -C "$QSTARS_CMFGEN_ROOT" d_new_main
