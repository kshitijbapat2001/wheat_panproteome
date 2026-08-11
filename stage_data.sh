#!/bin/bash
# =====================================================================================
# stage_data.sh  --  kept for backward compatibility.
#
# All staging (conda environment, proteomes, UniProt references, external database
# cache) now happens in setup_env.sh, which is the single command to run on a login
# node before submitting. This wrapper just calls it.
# =====================================================================================
set -euo pipefail
DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
echo "[stage_data] delegating to setup_env.sh ..."
exec bash "$DIR/setup_env.sh" "$@"
