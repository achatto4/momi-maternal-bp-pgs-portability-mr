#!/usr/bin/env bash
# Shared helpers for all stages. Source this at the top of each runner.
set -euo pipefail

log()  { printf '[%s] %s\n' "$(date '+%F %T')" "$*" >&2; }
die()  { printf '[%s] ERROR: %s\n' "$(date '+%F %T')" "$*" >&2; exit 1; }
need() { command -v "$1" >/dev/null 2>&1 || die "required tool not found: $1"; }

# run + echo a command, abort on failure
run() { log "RUN: $*"; "$@"; }
