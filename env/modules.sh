#!/usr/bin/env bash
# Source this to load cluster modules based on config values.
# Usage: source env/modules.sh config/config.yaml
# (Falls back to no-op if module command is unavailable.)
set -uo pipefail
CFG="${1:-config/config.yaml}"

# minimal yaml getter (key under compute:) — for quick use; real scripts use lib/yaml.sh
_get() { awk -v k="$1" '$1==k":"{ $1=""; sub(/^ /,""); gsub(/"/,""); print; exit }' "$CFG"; }

if command -v module >/dev/null 2>&1; then
  module load "$(_get plink1_module)" 2>/dev/null || true
  module load "$(_get plink2_module)" 2>/dev/null || true
  module load "$(_get r_module)"      2>/dev/null || true
else
  echo "[env] 'module' not found — assuming tools are on PATH." >&2
fi
