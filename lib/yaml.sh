#!/usr/bin/env bash
# Tiny YAML reader for flat scalar keys. For nested/list config the R/Python
# stages parse YAML properly; this is only for simple bash-side lookups.
# Usage: yaml_get <file> <dotted.key>   (best-effort, top-level + one nest)
yaml_get() {
  local f="$1" key="$2"
  python3 - "$f" "$key" <<'PY'
import sys, yaml
f, key = sys.argv[1], sys.argv[2]
d = yaml.safe_load(open(f))
for part in key.split('.'):
    d = d[int(part)] if part.isdigit() else d[part]
print(d)
PY
}
