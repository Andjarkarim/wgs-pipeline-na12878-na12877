import csv
import os
import re

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

BLAU, AKZENT, TEXT, GITTER, GRAU = "#2f6db0", "#d9742a", "#222222", "#e3e3e3", "#8a8a8a"
plt.rcParams.update({"font.size": 10, "axes.edgecolor": "#999999", "text.color": TEXT,
                     "axes.labelcolor": TEXT, "xtick.color": TEXT, "ytick.color": TEXT})

STRENG = re.compile(r"^(Pathogenic|Likely_pathogenic)")
PROTEIN = re.compile(r"missense|stop_gained|frameshift|splice_acceptor|splice_donor|start_lost|stop_lost|inframe")
PROBEN = [("Mutter (NA12878, ERR194147)", "ERR194147"), ("Vater (NA12877, ERR194146)", "ERR194146")]
STIL = {"belastbar": (BLAU, "o"), "unsicher": (AKZENT, "s"), "artefaktverdacht": (GRAU, "X")}


def qualitaet(gt, ab, dp, gq, ref="", alt=""):
    try:
        ab, dp, gq = float(ab), int(dp), int(gq)
    except ValueError:
        return "unsicher"
    if gt in ("1/1", "1|1"):
        if ab < 0.5:
            return "artefaktverdacht"
        return "belastbar" if (ab >= 0.85 and dp >= 15 and gq >= 30) else "unsicher"
    if ab < 0.25:
        indel = any(len(a) != len(ref) for a in alt.split(","))
        return "unsicher" if indel else "artefaktverdacht"
    if 0.35 <= ab <= 0.75 and dp >= 15 and gq >= 30:
        return "belastbar"
    return "unsicher"


def lade(sample):
    zeilen, gesehen = [], set()
    with open(f"varianten/{sample}.pathogene_tabelle.tsv", encoding="utf-8") as f:
        for z in csv.DictReader(f, delimiter="\t"):
            kons = z["KONSEQUENZ"]
            if not (STRENG.match(z["CLINVAR"]) and PROTEIN.search(kons) and "NMD" not in kons):
                continue
            k = (z["CHROM"], z["POS"], z["REF"], z["ALT"], z["GEN"])
            if k in gesehen:
                continue
            gesehen.add(k)
            try:
                ab = float(z["AB"])
            except ValueError:
                continue
            q = qualitaet(z["GENOTYP"], z["AB"], z["DP"], z["GQ"], z["REF"], z["ALT"])
            zeilen.append((z["GEN"], f'{z["CHROM"]}:{z["POS"]}', ab, q))
    return sorted(zeilen, key=lambda r: r[2])


daten = [(titel, lade(s)) for titel, s in PROBEN]
n_max = max(len(z) for _, z in daten)
fig, axes = plt.subplots(1, 2, figsize=(11.5, 0.42 * n_max + 2.6), sharex=True)

for ax, (titel, zeilen) in zip(axes, daten):
    ys = list(range(len(zeilen)))
    ax.axvspan(0.35, 0.75, color=BLAU, alpha=0.09, lw=0)
    ax.axvspan(0.25, 0.35, color=AKZENT, alpha=0.12, lw=0)
    ax.axvline(0.5, color="#777777", lw=1, ls=(0, (4, 3)))
    for y, (gen, pos, ab, q) in zip(ys, zeilen):
        farbe, marker = STIL[q]
        ax.scatter(ab, y, s=55, color=farbe, marker=marker, zorder=3)
        ax.text(ab + 0.025, y, f"{ab:.2f}", va="center", fontsize=8.5, color="#555555")
    ax.set_yticks(ys)
    ax.set_yticklabels([f"{g}  {p}" for g, p, _, _ in zeilen], fontsize=9)
    ax.set_xlim(0, 1)
    ax.set_ylim(-0.7, n_max - 0.3)
    ax.set_xlabel("Allelbalance (Anteil Reads mit alternativem Allel)")
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.xaxis.grid(True, color=GITTER, lw=0.8)
    ax.set_axisbelow(True)
    n_v = sum(1 for r in zeilen if r[3] == "artefaktverdacht")
    ax.set_title(f"{titel}: {len(zeilen)} Kandidaten, {n_v} verworfen", loc="left",
                 fontsize=9.5, color="#555555")

griff = [Line2D([], [], marker=m, color="none", markerfacecolor=c, markeredgecolor=c, markersize=8, label=l)
         for l, (c, m) in (("belastbar (AB 0,35–0,75, DP ≥ 15, GQ ≥ 30)", STIL["belastbar"]),
                           ("unsicher (auch Indels mit AB < 0,25)", STIL["unsicher"]),
                           ("Artefaktverdacht (SNV, AB < 0,25, verworfen)", STIL["artefaktverdacht"]))]
fig.suptitle("Allelbalance der ClinVar-Kandidaten vor dem Qualitätsfilter", x=0.01, ha="left",
             fontsize=13, fontweight="bold")
fig.legend(handles=griff, loc="lower center", ncol=3, frameon=False, fontsize=9)
fig.tight_layout(rect=(0, 0.07, 1, 0.95))
os.makedirs("abbildungen", exist_ok=True)
fig.savefig("abbildungen/11_allelbalance_eltern.png", dpi=200, facecolor="white")
print("Geschrieben: abbildungen/11_allelbalance_eltern.png")
for titel, zeilen in daten:
    print(titel)
    for g, p, ab, q in zeilen:
        print(f"  {g:9s} {p:14s} AB={ab:.2f}  {q}")
