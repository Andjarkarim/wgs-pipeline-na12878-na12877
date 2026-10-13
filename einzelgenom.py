import gzip, json, os, subprocess, sys
from collections import Counter
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
FAI = "referenz/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz.fai"
KAND_JSON = f"{OUT}/kandidaten_werte.json"
ZAHLEN = f"{OUT}/einzelgenom.json"
os.makedirs("abbildungen", exist_ok=True)

BIN = 1_000_000
CHROMS = [str(i) for i in range(1, 23)] + ["X", "Y"]
BLAU, BLAU_DUNKEL, BLAU_HELL, AKZENT, TEXT, GRAU = "#2f6db0", "#1b4676", "#9bbbe0", "#d9742a", "#222222", "#7a7a7a"
KOMPL = {"A": "T", "C": "G", "G": "C", "T": "A"}
TRANSITION = {"AG", "GA", "CT", "TC"}
KLASSEN = ["C>A", "C>G", "C>T", "T>A", "T>C", "T>G"]
plt.rcParams.update({"font.size": 10, "text.color": TEXT, "axes.labelcolor": TEXT,
                     "xtick.color": TEXT, "ytick.color": TEXT})


def zahl(n):
    return f"{int(n):,}".replace(",", ".")


# ---------- Daten sammeln ----------
def lese_laengen():
    laenge = {}
    for line in open(FAI):
        p = line.split("\t")
        if p[0] in CHROMS:
            laenge[p[0]] = int(p[1])
    return laenge


def zaehle_vcf(laenge):
    nb = {c: laenge[c] // BIN + 1 for c in CHROMS}
    bins = {c: {k: [0] * nb[c] for k in ("snv", "indel", "het", "hom")} for c in CHROMS}
    typ, subst, indel_len = Counter(), Counter(), Counter()
    ti = tv = het = hom = 0
    with gzip.open(RAW, "rt") as f:
        for line in f:
            if line[0] == "#":
                continue
            c = line.split("\t", 10)
            chrom = c[0]
            if chrom not in bins:
                continue
            ref, alt = c[3], c[4]
            b = int(c[1]) // BIN
            al = c[9].split(":", 1)[0].strip().replace("|", "/").split("/")
            ist_het = ist_hom = False
            if len(al) == 2 and "." not in al:
                ist_het, ist_hom = al[0] != al[1], al[0] == al[1]
            snvartig = False
            if "," in alt:
                typ["Multiallelisch"] += 1
                snvartig = len(ref) == 1 and all(len(a) == 1 for a in alt.split(","))
            elif len(ref) == 1 and len(alt) == 1:
                typ["SNV"] += 1
                snvartig = True
                r, a = ref.upper(), alt.upper()
                if r in KOMPL and a in KOMPL:
                    if (r + a) in TRANSITION:
                        ti += 1
                    else:
                        tv += 1
                    if r in "AG":
                        r, a = KOMPL[r], KOMPL[a]
                    subst[f"{r}>{a}"] += 1
            elif len(alt) > len(ref):
                typ["Insertion"] += 1
                indel_len[min(len(alt) - len(ref), 15)] += 1
            elif len(alt) < len(ref):
                typ["Deletion"] += 1
                indel_len[max(len(alt) - len(ref), -15)] += 1
            else:
                typ["MNP / komplex"] += 1
            bins[chrom]["snv" if snvartig else "indel"][b] += 1
            if ist_het:
                bins[chrom]["het"][b] += 1
                het += 1
            elif ist_hom:
                bins[chrom]["hom"][b] += 1
                hom += 1
    return dict(bins=bins, typ=dict(typ), subst=dict(subst), ti=ti, tv=tv, het=het, hom=hom,
                indel_len={str(k): v for k, v in indel_len.items()})


def anzahl_auffaelligkeiten(r):
    snp = len(r["ref"]) == 1 and len(r["alt"]) == 1
    w = lambda k: r.get(k)
    tests = [(w("qual"), lambda x: x < 30, True), (w("QD"), lambda x: x < 2, True),
             (w("SOR"), lambda x: x > 3, snp), (w("FS"), lambda x: x > (60 if snp else 200), True),
             (w("MQ"), lambda x: x < 40, snp), (w("alt_anteil"), lambda x: not 0.3 <= x <= 0.7, snp),
             (w("mapq_mittel"), lambda x: x < 40, True), (w("mapq_schlecht"), lambda x: x > 10, True)]
    return sum(1 for wert, f, gilt in tests if gilt and wert is not None and f(wert))


def lese_kandidaten():
    kand = []
    fehler = {}
    if os.path.exists(KAND_JSON):
        for r in json.load(open(KAND_JSON)):
            r["fehler"] = anzahl_auffaelligkeiten(r)
            fehler[(r["chrom"], r["pos"])] = r
    with open(KAND) as f:
        next(f)
        for line in f:
            c = line.rstrip("\n").split("\t")
            kand.append(dict(chrom=c[0], pos=int(c[1]), gen=c[4], r=fehler.get((c[0], int(c[1])))))
    return kand


def coverage():
    out = subprocess.run(["samtools", "idxstats", BAM], capture_output=True, text=True, check=True).stdout
    daten = {}
    for line in out.splitlines():
        p = line.split("\t")
        if p[0] in CHROMS:
            daten[p[0]] = (int(p[1]), int(p[2]))
    sam = subprocess.run(["samtools", "view", BAM, "1:100000000-100005000"], capture_output=True,
                         text=True, check=True).stdout.splitlines()
    rl = float(np.mean([len(l.split("\t")[9]) for l in sam[:2000]]))
    return {c: m * rl / l for c, (l, m) in daten.items()}, rl


# ---------- Abbildung 7: Circos ----------
def plot_circos(d, laenge, kand):
    from pycirclize import Circos
    circos = Circos({c: laenge[c] for c in CHROMS}, space=1.5)
    alle = np.concatenate([d["bins"][c]["snv"] for c in CHROMS])
    vmax = float(np.percentile(alle, 99.5))
    for sector in circos.sectors:
        c = sector.name
        b = d["bins"][c]
        nb = len(b["snv"])
        x = [min((i + 0.5) * BIN, sector.size) for i in range(nb)]
        außen = sector.add_track((92, 100))
        außen.axis(fc="#eeeeee", ec="#999999", lw=0.5)
        außen.text(c, size=8)
        außen.xticks_by_interval(50_000_000, outer=True, tick_length=1.2, show_label=False,
                                 line_kws=dict(lw=0.6, color="#777777"))
        t1 = sector.add_track((72, 90))
        t1.axis(fc="white", ec="#dddddd", lw=0.4)
        t1.fill_between(x, np.minimum(b["snv"], vmax), vmin=0, vmax=vmax, color=BLAU, lw=0)
        t2 = sector.add_track((52, 70))
        t2.axis(fc="white", ec="#dddddd", lw=0.4)
        xs = [x[i] for i in range(nb) if b["het"][i] + b["hom"][i] >= 10]
        ys = [b["het"][i] / (b["het"][i] + b["hom"][i]) for i in range(nb) if b["het"][i] + b["hom"][i] >= 10]
        if xs:
            t2.bar(xs, ys, width=BIN * 0.9, vmin=0, vmax=1, color=GRAU, lw=0)
        t3 = sector.add_track((40, 50))
        t3.axis(fc="none", ec="#eeeeee", lw=0.4)
        gruppen = {}
        for k in kand:
            if k["chrom"] == c:
                gruppen.setdefault(k["pos"] // BIN, []).append(k)
        for _, ks in gruppen.items():
            pos = ks[0]["pos"]
            ok = all(k["r"] is None or k["r"].get("fehler", 0) == 0 for k in ks)
            t3.scatter([pos], [0.5], vmin=0, vmax=1, s=55, zorder=5, lw=1.6,
                       facecolor=AKZENT if ok else "white", edgecolor=AKZENT)
            außen.annotate(pos, " / ".join(sorted({k["gen"] for k in ks})), label_size=7,
                           min_r=101, max_r=108)
    circos.text(NAME + "\n" + SAMPLE, size=11)
    fig = circos.plotfig(dpi=200, figsize=(10, 10))
    fig.legend(handles=[
        Patch(fc=BLAU, label="Außen: SNV-Dichte pro Mb"),
        Patch(fc=GRAU, label="Mitte: Anteil heterozygoter Varianten pro Mb"),
        Line2D([], [], marker="o", ls="", color=AKZENT, markeredgecolor="white", markersize=8,
               label="Kandidat ohne Auffälligkeit"),
        Line2D([], [], marker="o", ls="", markerfacecolor="white", markeredgecolor=AKZENT,
               markeredgewidth=1.6, markersize=8, label="Kandidat mit Auffälligkeit")],
        loc="lower center", ncol=2, frameon=False, fontsize=9)
    fig.suptitle("Zirkuläre Genomkarte (Chromosom 1–22, X, Y)", x=0.02, y=0.97, ha="left", fontsize=13,
                 fontweight="bold")
    fig.savefig(f"{OUT}/7_circos.png", dpi=200, facecolor="white")
    plt.close(fig)


# ---------- Abbildung 8: Variantenstatistik ----------
def plot_statistik(d, laenge):
    fig, axs = plt.subplots(2, 2, figsize=(13, 9))

    ax = axs[0, 0]
    namen = [n for n in ["SNV", "Insertion", "Deletion", "MNP / komplex", "Multiallelisch"] if d["typ"].get(n)]
    werte = [d["typ"][n] for n in namen]
    ys = list(range(len(namen)))[::-1]
    ax.barh(ys, werte, color=BLAU, height=0.6)
    ax.set_yticks(ys)
    ax.set_yticklabels(namen)
    for y, w in zip(ys, werte):
        ax.text(w + max(werte) * 0.01, y, zahl(w), va="center", fontsize=9)
    ax.set_xlim(0, max(werte) * 1.18)
    ax.set_title("a  Variantentypen", loc="left", fontweight="bold", fontsize=11)
    ax.set_xlabel("Anzahl Varianten")
    ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: zahl(v)))

    ax = axs[0, 1]
    sub = [d["subst"].get(k, 0) for k in KLASSEN]
    gesamt = sum(sub)
    anteile = [100 * s / gesamt for s in sub]
    farben = [BLAU_DUNKEL if k in ("C>T", "T>C") else BLAU_HELL for k in KLASSEN]
    ax.bar(range(6), anteile, color=farben, width=0.65)
    ax.set_xticks(range(6))
    ax.set_xticklabels(KLASSEN)
    for i, a in enumerate(anteile):
        ax.text(i, a + 0.6, f"{a:.1f} %".replace(".", ","), ha="center", fontsize=9)
    ax.set_ylim(0, max(anteile) * 1.18)
    ax.set_ylabel("Anteil an allen SNVs (%)")
    ax.set_title("b  Substitutionsspektrum", loc="left", fontweight="bold", fontsize=11)
    ax.legend(handles=[Patch(fc=BLAU_DUNKEL, label="Transition"), Patch(fc=BLAU_HELL, label="Transversion")],
              frameon=False, loc="upper right", fontsize=9)
    titv = d["ti"] / d["tv"] if d["tv"] else float("nan")
    ax.text(0.02, 0.96, f"Ti/Tv = {titv:.2f}".replace(".", ","), transform=ax.transAxes, va="top",
            fontsize=11, fontweight="bold")

    ax = axs[1, 0]
    il = {int(k): v for k, v in d["indel_len"].items()}
    xs = [k for k in range(-15, 16) if k != 0]
    ax.bar(xs, [il.get(k, 0) for k in xs], color=BLAU, width=0.8)
    ax.set_xticks([-15, -10, -5, -1, 1, 5, 10, 15])
    ax.set_xticklabels(["≤−15", "−10", "−5", "−1", "+1", "+5", "+10", "≥15"])
    ax.set_xlabel("Längenänderung in bp  (negativ = Deletion, positiv = Insertion)")
    ax.set_ylabel("Anzahl Indels")
    ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: zahl(v)))
    ax.set_title("c  Indel-Längen", loc="left", fontweight="bold", fontsize=11)

    ax = axs[1, 1]
    hf = []
    for c in CHROMS:
        h, o = sum(d["bins"][c]["het"]), sum(d["bins"][c]["hom"])
        hf.append(100 * h / (h + o) if h + o else 0)
    ax.bar(range(len(CHROMS)), hf, color=BLAU, width=0.7)
    ax.set_xticks(range(len(CHROMS)))
    ax.set_xticklabels(CHROMS, fontsize=8)
    ax.set_ylabel("Anteil heterozygoter Varianten (%)")
    ax.set_ylim(0, 100)
    ax.set_title("d  Heterozygotie pro Chromosom", loc="left", fontweight="bold", fontsize=11)
    ax.text(0.99, 0.97, f"Het/Hom gesamt = {d['het'] / d['hom']:.2f}".replace(".", ","),
            transform=ax.transAxes, ha="right", va="top", fontsize=10, fontweight="bold",
            bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="#cccccc"))

    for a in axs.ravel():
        for s in ("top", "right"):
            a.spines[s].set_visible(False)
        a.grid(True, axis="y" if a is not axs[0, 0] else "x", color="#e6e6e6", linewidth=0.8)
        a.set_axisbelow(True)
    fig.suptitle(f"Variantenstatistik von {SAMPLE} ({NAME})", x=0.01, ha="left", fontsize=13,
                 fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(f"{OUT}/8_variantenstatistik.png", dpi=200, facecolor="white")
    plt.close(fig)


# ---------- Abbildung 9: Coverage pro Chromosom ----------
def plot_coverage(tiefe, rl):
    aut = [tiefe[str(i)] for i in range(1, 23)]
    med = float(np.median(aut))
    rel = [tiefe[c] / med for c in CHROMS]
    fig, ax = plt.subplots(figsize=(11, 4.8))
    ax.bar(range(len(CHROMS)), rel, color=BLAU, width=0.7)
    ax.axhline(1, color=GRAU, lw=1, ls="--")
    ax.set_xticks(range(len(CHROMS)))
    ax.set_xticklabels(CHROMS)
    for i, r in enumerate(rel):
        ax.text(i, r + 0.02, f"{r:.2f}".replace(".", ","), ha="center", fontsize=8)
    ax.set_ylim(0, max(rel) * 1.15)
    ax.set_ylabel("Tiefe relativ zum Median der Autosomen")
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.grid(True, axis="y", color="#e6e6e6", linewidth=0.8)
    ax.set_axisbelow(True)
    fig.suptitle("Sequenziertiefe pro Chromosom", x=0.01, ha="left", fontsize=13, fontweight="bold")
    ax.set_title(f"Geschätzt aus gemappten Reads (Leselänge {rl:.0f} bp, Duplikate inklusive): Autosomen ca. {med:.0f}-fach.\n"
                 "Erwartung weiblich: X ≈ 1, Y ≈ 0 · männlich: X ≈ 0,5, Y ≈ 0,5",
                 loc="left", fontsize=9, color="#555555")
    fig.tight_layout()
    fig.savefig(f"{OUT}/9_coverage_chromosomen.png", dpi=200, facecolor="white")
    plt.close(fig)
    return med, dict(zip(CHROMS, rel))


# ---------- Ablauf ----------
laenge = lese_laengen()
if os.path.exists(ZAHLEN) and "--neu" not in sys.argv:
    d = json.load(open(ZAHLEN))
else:
    print("Zähle Varianten (dauert einige Minuten) ...", flush=True)
    d = zaehle_vcf(laenge)
    json.dump(d, open(ZAHLEN, "w"))
kand = lese_kandidaten()
plot_circos(d, laenge, kand)
plot_statistik(d, laenge)
tiefe, rl = coverage()
med, rel = plot_coverage(tiefe, rl)
print("Variantentypen:", d["typ"])
print(f"Ti/Tv = {d['ti'] / d['tv']:.2f}   Het/Hom = {d['het'] / d['hom']:.2f}")
print(f"Tiefe Autosomen (Median) ≈ {med:.0f}x, Leselänge ≈ {rl:.0f} bp")
print("Relative Tiefe X = {:.2f}, Y = {:.3f}".format(rel["X"], rel["Y"]))
print("Fertig: abbildungen/7_circos.png, 8_variantenstatistik.png, 9_coverage_chromosomen.png")
