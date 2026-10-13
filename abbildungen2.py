import gzip, json, os, subprocess, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

import os
SAMPLE = os.environ.get("SAMPLE", "ERR194147")
NAME = os.environ.get("NAME", "NA12878")
OUT = f"abbildungen/{SAMPLE}"
os.makedirs(OUT, exist_ok=True)
RAW = f"varianten/{SAMPLE}.vcf.gz"
KAND = f"varianten/{SAMPLE}.kandidaten_streng.tsv"
BAM = f"alignment/{SAMPLE}.dedup.bam"
CLINVAR = "referenz/clinvar.vcf.gz"
FAI = "referenz/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz.fai"
DICHTE_JSON = f"{OUT}/dichte.json"
KAND_JSON = f"{OUT}/kandidaten_werte.json"
os.makedirs("abbildungen", exist_ok=True)

BIN = 1_000_000
CHROMS = [str(i) for i in range(1, 23)] + ["X", "Y"]
BLAU_HELL, AKZENT, TEXT, GRAU = "#dbe7f5", "#d9742a", "#222222", "#e6e6e6"
STERNE = {"practice_guideline": 4, "reviewed_by_expert_panel": 3,
          "criteria_provided,_multiple_submitters,_no_conflicts": 2,
          "criteria_provided,_single_submitter": 1,
          "criteria_provided,_conflicting_classifications": 1,
          "criteria_provided,_conflicting_interpretations": 1}


# ---------- Daten sammeln ----------
def zaehle_dichte():
    laenge = {}
    for line in open(FAI):
        p = line.split("\t")
        if p[0] in CHROMS:
            laenge[p[0]] = int(p[1])
    z = {c: [0] * (laenge[c] // BIN + 1) for c in CHROMS}
    with gzip.open(RAW, "rt") as f:
        for line in f:
            if line[0] == "#":
                continue
            c, pos = line.split("\t", 2)[:2]
            if c in z:
                z[c][int(pos) // BIN] += 1
    return z


def lese_kandidaten():
    rows = []
    with open(KAND) as f:
        next(f)
        for line in f:
            c = line.rstrip("\n").split("\t")
            rows.append(dict(chrom=c[0], pos=int(c[1]), ref=c[2], alt=c[3], gen=c[4]))
    return rows


def zahl_oder_none(info, k):
    try:
        return float(info[k])
    except (KeyError, ValueError):
        return None


def vcf_werte(rows):
    keys = {(r["chrom"], str(r["pos"])): r for r in rows}
    gefunden = 0
    with gzip.open(RAW, "rt") as f:
        for line in f:
            if line[0] == "#":
                continue
            c = line.rstrip("\n").split("\t")
            r = keys.get((c[0], c[1]))
            if r is None or c[3] != r["ref"] or c[4] != r["alt"]:
                continue
            info = dict(x.split("=", 1) for x in c[7].split(";") if "=" in x)
            fd = dict(zip(c[8].split(":"), c[9].split(":")))
            ad = [int(x) for x in fd["AD"].split(",")]
            r.update(qual=float(c[5]), gt=fd["GT"], alt_anteil=ad[1] / (ad[0] + ad[1]),
                     QD=zahl_oder_none(info, "QD"), FS=zahl_oder_none(info, "FS"),
                     SOR=zahl_oder_none(info, "SOR"), MQ=zahl_oder_none(info, "MQ"))
            gefunden += 1
            if gefunden == len(rows):
                break


def krankheit_kurz(clndn):
    teile = [t.replace("_", " ") for t in clndn.split("|") if t and t not in ("not_provided", "not_specified")]
    text = "; ".join(teile[:2]) if teile else "keine Angabe"
    return text if len(text) <= 52 else text[:49] + "…"


def clinvar_werte(rows):
    keys = {(r["chrom"], str(r["pos"]), r["ref"], r["alt"]): r for r in rows}
    with gzip.open(CLINVAR, "rt") as f:
        for line in f:
            if line[0] == "#":
                continue
            c = line.rstrip("\n").split("\t")
            r = keys.get((c[0], c[1], c[3], c[4]))
            if r is None:
                continue
            info = dict(x.split("=", 1) for x in c[7].split(";") if "=" in x)
            r["krankheit"] = krankheit_kurz(info.get("CLNDN", ""))
            r["sterne"] = STERNE.get(info.get("CLNREVSTAT", ""), 0)


def mapq_werte(r):
    reg = f'{r["chrom"]}:{r["pos"]}-{r["pos"]}'
    out = subprocess.run(["samtools", "view", BAM, reg], capture_output=True, text=True, check=True).stdout
    q = [int(l.split("\t")[4]) for l in out.splitlines()]
    r["reads"] = len(q)
    r["mapq_mittel"] = sum(q) / len(q) if q else None
    r["mapq_schlecht"] = 100 * sum(x < 20 for x in q) / len(q) if q else None


def bewerte(r):
    snp = len(r["ref"]) == 1 and len(r["alt"]) == 1
    s = {}

    def pruefe(k, fehler, gilt=True):
        w = r.get(k)
        s[k] = "na" if (not gilt or w is None) else ("fail" if fehler(w) else "ok")

    pruefe("qual", lambda w: w < 30)
    pruefe("QD", lambda w: w < 2)
    pruefe("SOR", lambda w: w > 3, snp)
    pruefe("FS", lambda w: w > (60 if snp else 200))
    pruefe("MQ", lambda w: w < 40, snp)
    pruefe("alt_anteil", lambda w: not 0.3 <= w <= 0.7, snp)
    pruefe("mapq_mittel", lambda w: w < 40)
    pruefe("mapq_schlecht", lambda w: w > 10)
    r["status"] = s
    r["fehler"] = sum(v == "fail" for v in s.values())


# ---------- Abbildungen ----------
plt.rcParams.update({"font.size": 10, "text.color": TEXT})


def plot_dichte(z, rows):
    n = max(len(v) for v in z.values())
    M = np.full((len(CHROMS), n), np.nan)
    for i, c in enumerate(CHROMS):
        M[i, : len(z[c])] = z[c]
    cmap = plt.get_cmap("Blues").copy()
    cmap.set_bad("white")
    fig, ax = plt.subplots(figsize=(13, 7))
    im = ax.imshow(np.ma.masked_invalid(M), aspect="auto", cmap=cmap, vmin=0,
                   vmax=np.nanpercentile(M, 99), interpolation="nearest")
    ax.set_yticks(range(len(CHROMS)))
    ax.set_yticklabels(CHROMS)
    ax.set_xlim(-3, n - 0.5)
    ax.set_xticks(range(0, n, 25))
    ax.set_xticklabels([str(x) for x in range(0, n, 25)])
    ax.set_xlabel("Position im Chromosom (Mb)")
    for s in ax.spines.values():
        s.set_visible(False)
    fig.colorbar(im, ax=ax, pad=0.01, fraction=0.025, label="Varianten pro Mb")
    ax.text(len(z["Y"]) + 2, CHROMS.index("Y"), "Y: Artefakte (weibliche Probe)", va="center",
            fontsize=8, color="#555555")
    markiert = {}
    for r in rows:
        k = (CHROMS.index(r["chrom"]), r["pos"] // BIN)
        markiert.setdefault(k, []).append(r)
    for (i, x), rs in markiert.items():
        ok = all(r["fehler"] == 0 for r in rs)
        ax.scatter([x], [i], s=120, zorder=3, linewidths=2,
                   facecolor=AKZENT if ok else "white", edgecolor="white" if ok else AKZENT)
        namen = " / ".join(sorted({r["gen"] for r in rs}))
        ax.text(x + 2.5, i, namen, va="center", fontsize=8, zorder=4,
                bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.85))
    ax.legend(handles=[
        Line2D([], [], marker="o", ls="", color=AKZENT, markeredgecolor="white", markersize=9,
               label="Kandidat ohne Auffälligkeit im Qualitätscheck"),
        Line2D([], [], marker="o", ls="", markerfacecolor="white", markeredgecolor=AKZENT,
               markeredgewidth=2, markersize=9, label="Kandidat mit ≥ 1 Auffälligkeit")],
        loc="lower right", frameon=True, facecolor="white", edgecolor="#cccccc", fontsize=9)
    fig.suptitle("Variantendichte im Genom mit den Kandidaten", x=0.01, ha="left", fontsize=13,
                 fontweight="bold")
    ax.set_title(f"{SAMPLE}: Varianten pro 1-Mb-Fenster, dunkler = mehr Varianten", loc="left",
                 fontsize=9.5, color="#555555")
    fig.tight_layout()
    fig.savefig(f"{OUT}/5_dichte_kandidaten.png", dpi=200, facecolor="white")
    plt.close(fig)


def plot_tabelle(rows):
    spalten = [("QUAL", "qual", "{:.0f}"), ("QD", "QD", "{:.1f}"), ("SOR", "SOR", "{:.1f}"),
               ("FS", "FS", "{:.1f}"), ("MQ", "MQ", "{:.0f}"), ("Alt-\nAnteil", "alt_anteil", "{:.2f}"),
               ("MAPQ\nØ", "mapq_mittel", "{:.0f}"), ("Reads\nMAPQ<20 (%)", "mapq_schlecht", "{:.0f}")]
    zeilen = sorted(rows, key=lambda r: (r["fehler"], CHROMS.index(r["chrom"]), r["pos"]))
    b_gen, b_krank, b_stern, b_m, b_f = 3.3, 6.8, 1.3, 1.25, 1.3
    x_krank = b_gen
    x_stern = x_krank + b_krank
    x_m0 = x_stern + b_stern
    x_f = x_m0 + b_m * len(spalten)
    total = x_f + b_f
    n = len(zeilen)
    fig, ax = plt.subplots(figsize=(16, 0.55 * n + 3.2))
    ax.set_xlim(0, total)
    ax.set_ylim(n + 0.2, -1.4)
    ax.axis("off")
    kopf = dict(va="center", ha="center", fontsize=9, fontweight="bold")
    ax.text(0.1, -0.55, "Gen · Position", va="center", ha="left", fontsize=9, fontweight="bold")
    ax.text(x_krank + 0.1, -0.55, "Erkrankung laut ClinVar", va="center", ha="left", fontsize=9,
            fontweight="bold")
    ax.text(x_stern + b_stern / 2, -0.55, "ClinVar-\nEvidenz", **kopf)
    for j, (k, _, _) in enumerate(spalten):
        ax.text(x_m0 + b_m * (j + 0.5), -0.55, k, **kopf)
    ax.text(x_f + b_f / 2, -0.55, "Auffällig-\nkeiten", **kopf)
    farbe = {"ok": BLAU_HELL, "fail": AKZENT, "na": GRAU}
    for i, r in enumerate(zeilen):
        ax.text(0.1, i + 0.5, f'{r["gen"]}  {r["chrom"]}:{r["pos"]:,}'.replace(",", "."), va="center",
                ha="left", fontsize=9.5, fontweight="bold")
        ax.text(x_krank + 0.1, i + 0.5, r.get("krankheit", "keine Angabe"), va="center", ha="left",
                fontsize=9)
        st = r.get("sterne", 0)
        ax.text(x_stern + b_stern / 2, i + 0.5, "★" * st if st else "–", va="center", ha="center",
                fontsize=9)
        for j, (_, k, fmt) in enumerate(spalten):
            s = r["status"][k]
            x = x_m0 + b_m * j
            ax.add_patch(plt.Rectangle((x + 0.05, i + 0.07), b_m - 0.1, 0.86, fc=farbe[s], ec="none"))
            w = r.get(k)
            ax.text(x + b_m / 2, i + 0.5, fmt.format(w) if w is not None else "–", va="center",
                    ha="center", fontsize=9, color="white" if s == "fail" else TEXT,
                    fontweight="bold" if s == "fail" else "normal")
        f = r["fehler"]
        ax.add_patch(plt.Rectangle((x_f + 0.05, i + 0.07), b_f - 0.1, 0.86,
                                   fc=AKZENT if f else BLAU_HELL, ec="none"))
        ax.text(x_f + b_f / 2, i + 0.5, str(f), va="center", ha="center", fontsize=10,
                color="white" if f else TEXT, fontweight="bold" if f else "normal")
    ax.legend(handles=[Patch(fc=BLAU_HELL, label="Schwelle eingehalten"),
                       Patch(fc=AKZENT, label="Schwelle verletzt"),
                       Patch(fc=GRAU, label="nicht anwendbar / kein Wert")],
              loc="upper left", bbox_to_anchor=(0, -0.01), ncol=3, frameon=False, fontsize=9)
    fig.suptitle("Die Kandidaten im Qualitätscheck", x=0.01, ha="left", fontsize=13, fontweight="bold")
    fig.text(0.01, 0.015,
             "Schwellen QUAL, QD, SOR, FS, MQ: GATK-Empfehlungen für Keimbahnvarianten (SNP: QD<2, QUAL<30, SOR>3, FS>60, MQ<40;\n"
             "Indel: QD<2, QUAL<30, FS>200). Alt-Anteil (SNP, erwartet 0,3–0,7) und MAPQ (Ø<40, >10 % der Reads mit MAPQ<20) sind\n"
             "Faustregeln. Bei Indels ist der Alt-Anteil durch HaplotypeCaller oft unterschätzt und wird nicht bewertet.",
             ha="left", va="bottom", fontsize=8, color="#555555")
    fig.subplots_adjust(left=0.01, right=0.99, top=0.91, bottom=0.14)
    fig.savefig(f"{OUT}/6_kandidaten_qualitaet.png", dpi=200, facecolor="white")
    plt.close(fig)


# ---------- Ablauf ----------
neu = "--neu" in sys.argv
if os.path.exists(KAND_JSON) and not neu:
    rows = json.load(open(KAND_JSON))
else:
    print("Lese Kandidaten, VCF, ClinVar und BAM ...", flush=True)
    rows = lese_kandidaten()
    vcf_werte(rows)
    clinvar_werte(rows)
    fehlt = [r["gen"] for r in rows if "qual" not in r]
    if fehlt:
        sys.exit(f"Kein VCF-Eintrag gefunden für: {fehlt}")
    for r in rows:
        mapq_werte(r)
    json.dump(rows, open(KAND_JSON, "w"), indent=1, ensure_ascii=False)
for r in rows:
    bewerte(r)

if os.path.exists(DICHTE_JSON) and not neu:
    z = json.load(open(DICHTE_JSON))
else:
    print("Zähle Variantendichte (dauert wenige Minuten) ...", flush=True)
    z = zaehle_dichte()
    json.dump(z, open(DICHTE_JSON, "w"))

plot_dichte(z, rows)
plot_tabelle(rows)
for r in sorted(rows, key=lambda r: (r["fehler"], r["gen"])):
    print(f'{r["gen"]:9s} {r["chrom"]}:{r["pos"]}  Auffälligkeiten={r["fehler"]}  {r.get("krankheit", "keine Angabe")}')
print("Fertig: abbildungen/5_dichte_kandidaten.png, abbildungen/6_kandidaten_qualitaet.png")
