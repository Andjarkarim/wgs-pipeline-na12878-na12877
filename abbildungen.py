import gzip, json, os, re, sys
from collections import Counter

import os
SAMPLE = os.environ.get("SAMPLE", "ERR194147")
NAME = os.environ.get("NAME", "NA12878")
OUT = f"abbildungen/{SAMPLE}"
os.makedirs(OUT, exist_ok=True)
RAW = f"varianten/{SAMPLE}.vcf.gz"
VEP = f"varianten/{SAMPLE}.vep.vcf.gz"
ZAHLEN = f"{OUT}/zahlen.json"
os.makedirs("abbildungen", exist_ok=True)

CHROMS = [str(i) for i in range(1, 23)] + ["X", "Y", "MT"]

# Schweregrad-Gruppen, schwerste zuerst (vereinfacht gegenüber der VEP-Rangfolge)
GRUPPEN = [
    ("Funktionsverlust (Stop, Frameshift, Spleißstelle)",
     {"transcript_ablation", "splice_acceptor_variant", "splice_donor_variant",
      "stop_gained", "frameshift_variant", "stop_lost", "start_lost"}),
    ("Missense / Inframe",
     {"missense_variant", "inframe_insertion", "inframe_deletion", "protein_altering_variant"}),
    ("Synonym",
     {"synonymous_variant", "start_retained_variant", "stop_retained_variant"}),
    ("UTR",
     {"5_prime_UTR_variant", "3_prime_UTR_variant"}),
    ("Intron / Spleißregion",
     {"intron_variant", "splice_region_variant", "splice_donor_5th_base_variant",
      "splice_donor_region_variant", "splice_polypyrimidine_tract_variant"}),
    ("Nicht-kodierende RNA",
     {"non_coding_transcript_exon_variant", "non_coding_transcript_variant",
      "mature_miRNA_variant"}),
    ("Up-/Downstream",
     {"upstream_gene_variant", "downstream_gene_variant"}),
    ("Intergenisch",
     {"intergenic_variant", "regulatory_region_variant", "TF_binding_site_variant"}),
]
TERM2RANK = {t: i for i, (_, ts) in enumerate(GRUPPEN) for t in ts}

STRENG = re.compile(r"^(Pathogenic|Likely_pathogenic)")
PROTEIN = re.compile(r"missense|stop_gained|frameshift|splice_acceptor|splice_donor|start_lost|stop_lost|inframe")


def zaehle_roh():
    gesamt, chrom = 0, Counter()
    with gzip.open(RAW, "rt") as f:
        for line in f:
            if line[0] == "#":
                continue
            gesamt += 1
            name = line.split("\t", 1)[0]
            chrom[name if name in CHROMS else "Scaffolds"] += 1
    return gesamt, dict(chrom)


def zaehle_vep():
    d = dict(n_vep=0, mit_clinvar=0, pathogen_inkl=0, streng=0, streng_protein=0)
    konsequenz, clin, clin_fein = Counter(), Counter(), Counter()
    i_cons = i_clin = None
    with gzip.open(VEP, "rt") as f:
        for line in f:
            if line.startswith("##INFO=<ID=CSQ"):
                felder = line.split("Format: ")[1].rstrip('">\n').split("|")
                i_cons, i_clin = felder.index("Consequence"), felder.index("ClinVar_CLNSIG")
                continue
            if line[0] == "#":
                continue
            d["n_vep"] += 1
            info = line.split("\t", 8)[7]
            k = info.find("CSQ=")
            if k < 0:
                continue
            bester, cs_var, ist_protein = 99, "", False
            for e in info[k + 4:].split(";")[0].split(","):
                p = e.split("|")
                cons = p[i_cons]
                for t in cons.split("&"):
                    r = TERM2RANK.get(t)
                    if r is not None and r < bester:
                        bester = r
                cs = p[i_clin] if len(p) > i_clin else ""
                if cs and not cs_var:
                    cs_var = cs
                if STRENG.match(cs) and PROTEIN.search(cons) and "NMD" not in cons:
                    ist_protein = True
            konsequenz[GRUPPEN[bester][0] if bester < 99 else "Andere"] += 1
            if cs_var:
                d["mit_clinvar"] += 1
                tok = cs_var.split("&")[0]
                clin_fein[tok] += 1
                if tok.startswith("Conflicting"):
                    g = "Widersprüchlich"
                elif STRENG.match(tok):
                    g = "Pathogen / wahrscheinlich pathogen"
                    d["streng"] += 1
                elif "enign" in tok:
                    g = "Benign"
                elif tok.startswith("Uncertain"):
                    g = "Unklar (VUS)"
                else:
                    g = "Andere"
                clin[g] += 1
                if "athogenic" in cs_var:
                    d["pathogen_inkl"] += 1
            if ist_protein:
                d["streng_protein"] += 1
    d.update(konsequenz=dict(konsequenz), clinvar=dict(clin), clinvar_fein=dict(clin_fein))
    return d


if os.path.exists(ZAHLEN) and "--neu" not in sys.argv:
    daten = json.load(open(ZAHLEN))
else:
    print("Zähle Rohdatei ...", flush=True)
    gesamt, chrom = zaehle_roh()
    print("Zähle annotierte Datei (dauert einige Minuten) ...", flush=True)
    daten = zaehle_vep()
    daten.update(gesamt=gesamt, chrom=chrom)
    json.dump(daten, open(ZAHLEN, "w"), indent=1, ensure_ascii=False)

print(json.dumps({k: v for k, v in daten.items() if k not in ("clinvar_fein",)}, indent=1, ensure_ascii=False))

# ---------------- Abbildungen ----------------
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BLAU, AKZENT, TEXT, GITTER = "#2f6db0", "#d9742a", "#222222", "#e3e3e3"
plt.rcParams.update({"font.size": 10, "axes.edgecolor": "#999999", "text.color": TEXT,
                     "axes.labelcolor": TEXT, "xtick.color": TEXT, "ytick.color": TEXT})


def zahl(n):
    return f"{n:,}".replace(",", ".")


def stil(ax, achse):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    getattr(ax, f"{achse}axis").grid(True, color=GITTER, linewidth=0.8)
    ax.set_axisbelow(True)


def horizontal(namen, werte, farben, titel, datei, untertitel):
    fig, ax = plt.subplots(figsize=(9.5, 0.6 * len(namen) + 1.8))
    ys = list(range(len(namen)))[::-1]
    ax.barh(ys, werte, color=farben, height=0.6)
    ax.set_yticks(ys)
    ax.set_yticklabels(namen)
    ax.set_xscale("log")
    ax.set_xlim(1, max(werte) * 12)
    for y, w in zip(ys, werte):
        ax.text(w * 1.2, y, zahl(w), va="center", fontsize=10)
    ax.set_xlabel("Anzahl Varianten (logarithmische Achse)")
    stil(ax, "x")
    fig.suptitle(titel, x=0.01, ha="left", fontsize=13, fontweight="bold")
    ax.set_title(untertitel, loc="left", fontsize=9.5, color="#555555")
    fig.tight_layout()
    fig.savefig(f"{OUT}/{datei}", dpi=200, facecolor="white")
    plt.close(fig)


# 1 Trichter
st = [("Alle Varianten (HaplotypeCaller)", daten["gesamt"]),
      ("Von VEP annotiert (Chromosomen)", daten["n_vep"]),
      ("Mit ClinVar-Eintrag", daten["mit_clinvar"]),
      ("ClinVar „pathogen“ inkl. widersprüchlich", daten["pathogen_inkl"]),
      ("Pathogen / wahrscheinlich pathogen", daten["streng"]),
      ("… und proteinverändernd", daten["streng_protein"]),
      ("Nach Qualitätsfilter (AB, DP, GQ)", sum(1 for _ in open(f"varianten/{SAMPLE}.kandidaten_streng.tsv")) - 1)]
st = [(n, w) for n, w in st if w > 0]
horizontal([n for n, _ in st], [w for _, w in st],
           [BLAU] * (len(st) - 1) + [AKZENT], "Vom Genom zu den Kandidaten",
           "1_trichter.png", f"{SAMPLE}: Anzahl Varianten pro Filterstufe")

# 2 Konsequenzen
reihenfolge = [n for n, _ in GRUPPEN] + ["Andere"]
ks = [(n, daten["konsequenz"][n]) for n in reihenfolge if daten["konsequenz"].get(n)]
horizontal([n for n, _ in ks], [w for _, w in ks], [BLAU] * len(ks),
           "Funktionelle Folge der Varianten (VEP)",
           "2_konsequenzen.png", "Pro Variante zählt die schwerste Folge über alle Transkripte")

# 3 ClinVar
rf = ["Benign", "Unklar (VUS)", "Widersprüchlich", "Pathogen / wahrscheinlich pathogen", "Andere"]
cs = [(n, daten["clinvar"][n]) for n in rf if daten["clinvar"].get(n)]
horizontal([n for n, _ in cs], [w for _, w in cs],
           [AKZENT if n.startswith("Pathogen") else BLAU for n, _ in cs],
           "ClinVar-Einstufung der Varianten",
           "3_clinvar.png", "Nur Varianten mit ClinVar-Eintrag")

# 4 Varianten pro Chromosom
namen = CHROMS + ["Scaffolds"]
werte = [daten["chrom"].get(n, 0) for n in namen]
fig, ax = plt.subplots(figsize=(11, 4.6))
ax.bar(range(len(namen)), werte, color=BLAU, width=0.7)
ax.set_xticks(range(len(namen)))
ax.set_xticklabels(namen, rotation=0, fontsize=9)
ax.set_ylabel("Anzahl Varianten")
ax.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: zahl(int(v))))
stil(ax, "y")
fig.suptitle("Varianten pro Chromosom", x=0.01, ha="left", fontsize=13, fontweight="bold")
fig.tight_layout()
fig.savefig(f"{OUT}/4_chromosomen.png", dpi=200, facecolor="white")
plt.close(fig)
print("Fertig: abbildungen/")
