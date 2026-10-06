#!/bin/bash
set -euo pipefail
DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
echo "[stage_data] delegating to setup_env.sh ..."
exec bash "$DIR/setup_env.sh" "$@"
