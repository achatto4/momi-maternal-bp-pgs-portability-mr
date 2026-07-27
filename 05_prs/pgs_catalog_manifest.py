#!/usr/bin/env python3
# ============================================================
# pgs_catalog_manifest.py
# Enumerate ALL PGS Catalog polygenic scores for SBP and DBP, with their training
# ancestry, variant count, and GRCh38 harmonized download URL. Run on JHPCE
# (needs internet). Output: a manifest TSV you can filter to an ancestry-specific
# set for the PRS-transferability study.
#
# Usage:
#   python3 pgs_catalog_manifest.py --out pgs_bp_manifest.tsv
#   (then inspect; pick one score per ancestry x trait, or all, for scoring)
# ============================================================
import json, sys, time, urllib.request, urllib.error, argparse

REST = "https://www.pgscatalog.org/rest"
TRAITS = {"SBP": "EFO_0006335", "DBP": "EFO_0006336"}
FTP = "https://ftp.ebi.ac.uk/pub/databases/spot/pgs/scores/{id}/ScoringFiles/Harmonized/{id}_hmPOS_GRCh38.txt.gz"

def get(url, tries=4):
    for i in range(tries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "momi-pipeline"})
            with urllib.request.urlopen(req, timeout=60) as r:
                return json.load(r)
        except (urllib.error.URLError, TimeoutError) as e:
            sys.stderr.write(f"  retry {i+1}/{tries} ({e}) {url}\n"); time.sleep(3*(i+1))
    raise SystemExit(f"FAILED: {url}")

def dominant_ancestry(score):
    """Best-effort: top ancestry of the GWAS (source) stage, else development stage."""
    ad = score.get("ancestry_distribution") or {}
    for stage in ("gwas", "dev", "eval"):
        dist = (ad.get(stage) or {}).get("dist") or {}
        if dist:
            top = max(dist.items(), key=lambda kv: kv[1])
            return f"{top[0]}({top[1]}%)", ";".join(f"{k}:{v}" for k, v in sorted(dist.items()))
    return "NA", "NA"

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="pgs_bp_manifest.tsv")
    args = ap.parse_args()
    rows = []
    for trait, efo in TRAITS.items():
        url = f"{REST}/score/search?trait_id={efo}"
        n = 0
        while url:
            page = get(url)
            for s in page.get("results", []):
                pid = s["id"]
                anc, dist = dominant_ancestry(s)
                rows.append([pid, trait, efo, anc, dist,
                             str(s.get("variants_number", "NA")),
                             (s.get("publication") or {}).get("firstauthor", "NA"),
                             (s.get("publication") or {}).get("PMID", "NA"),
                             FTP.format(id=pid)])
                n += 1
            url = page.get("next")
            time.sleep(0.3)
        sys.stderr.write(f"{trait} ({efo}): {n} scores\n")
    rows.sort(key=lambda r: (r[1], r[3]))
    with open(args.out, "w") as fh:
        fh.write("PGS_id\ttrait\tEFO\tgwas_ancestry\tancestry_dist\tn_variants\tfirst_author\tPMID\tGRCh38_url\n")
        for r in rows:
            fh.write("\t".join(map(str, r)) + "\n")
    sys.stderr.write(f"\nwrote {len(rows)} rows -> {args.out}\n")
    # quick ancestry x trait summary
    from collections import Counter
    c = Counter((r[1], r[3].split("(")[0]) for r in rows)
    sys.stderr.write("=== scores per (trait, ancestry) ===\n")
    for (t, a), k in sorted(c.items()):
        sys.stderr.write(f"  {t:4s} {a:6s} : {k}\n")

if __name__ == "__main__":
    main()
