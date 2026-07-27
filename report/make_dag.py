#!/usr/bin/env python3
"""
make_dag.py -- Figure 2, the assumed causal model, drawn with Graphviz (dot) for a clean,
publication-grade directed acyclic graph.

  python3 report/make_dag.py --out DIR      # writes F1B_dag.pdf (and .png)

Design: the instrument (polygenic score) acts only on maternal blood pressure; blood pressure
affects the perinatal outcome both directly and through hypertensive disorders of pregnancy
(a mediator, retained); measured confounders bias the observational association only;
population structure is the single path that would invalidate the instrument.
"""
import argparse, os, shutil
ap = argparse.ArgumentParser()
ap.add_argument("--out", required=True)
A = ap.parse_args()
os.makedirs(A.out, exist_ok=True)

try:
    import graphviz
    graphviz.Digraph().pipe(format="pdf")   # verifies the `dot` binary is present
except Exception as e:
    print(f"  F1B_dag SKIPPED (graphviz/dot not available: {e}); keeping existing PDF.\n"
          f"    install with:  brew install graphviz && pip3 install graphviz")
    raise SystemExit(0)

# Okabe-Ito, matching momi_figtheme
INSTR, EXPO, MED, OUT = "#D6E4F5", "#FFFFFF", "#FBE6C9", "#CBEBDA"
INSTR_E, MED_E, OUT_E = "#0072B2", "#E69F00", "#009E73"
LAT = "#EFEFEF"; LAT_E = "#9A9A9A"

g = graphviz.Digraph("causal_model", format="pdf")
g.attr(rankdir="LR", nodesep="0.5", ranksep="0.9", splines="spline",
       fontname="Helvetica", bgcolor="white", pad="0.25")
g.attr("node", shape="box", style="rounded,filled", fontname="Helvetica",
       fontsize="12", margin="0.18,0.11", penwidth="1.6", color="#333333")
g.attr("edge", color="#333333", penwidth="1.5", arrowsize="0.8", fontname="Helvetica",
       fontsize="10", fontcolor="#555555")

# --- nodes ---
g.node("Z",  "Blood-pressure\npolygenic score", fillcolor=INSTR, color=INSTR_E)
g.node("BP", "Maternal\nblood pressure",         fillcolor=EXPO,  color="#333333")
g.node("PE", "Hypertensive disorders\nof pregnancy", fillcolor=MED, color=MED_E)
g.node("Y",  "Perinatal outcome\n(BWT, LBW, SGA, PTB)", fillcolor=OUT, color=OUT_E)
g.node("U",  "Confounders\n(age, gravidity,\nBMI, socioeconomic)",
       fillcolor=LAT, color=LAT_E, style="rounded,filled,dashed")
g.node("PS", "Population\nstructure", fillcolor=LAT, color=LAT_E,
       style="rounded,filled,dashed")

# --- edges ---
g.edge("Z", "BP", label="instrument")
g.edge("BP", "Y", label="direct effect", penwidth="2.2")   # main causal path
g.edge("BP", "PE")                                          # mediated path
g.edge("PE", "Y")
g.edge("U", "BP", style="dashed", color=LAT_E)             # confounding (observational only)
g.edge("U", "Y",  style="dashed", color=LAT_E)
g.edge("PS", "Z", style="dashed", color=LAT_E)             # would violate the IV assumption
g.edge("PS", "Y", style="dashed", color=LAT_E)

# keep the main chain on one horizontal rank; PE above it, U above, PS below
with g.subgraph() as s:
    s.attr(rank="same"); s.node("Z"); s.node("BP"); s.node("Y")

import tempfile
tmp = tempfile.mkdtemp(); base = os.path.join(tmp, "F1B_dag")
g.render(base, cleanup=True)                       # pdf
g.format = "png"; g.attr(dpi="300"); g.render(base, cleanup=True)
for ext in ("pdf", "png"):
    shutil.copy2(base + "." + ext, os.path.join(A.out, "F1B_dag." + ext))
print("wrote", os.path.join(A.out, "F1B_dag") + ".pdf / .png")
