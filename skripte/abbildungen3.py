import json, re
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import os
SAMPLE = os.environ.get("SAMPLE", "ERR194147")
NAME = os.environ.get("NAME", "NA12878")
OUT = f"abbildungen/{SAMPLE}"
os.makedirs(OUT, exist_ok=True)
BLAU, BLAU_HELL, AKZENT, TEXT, GRAU = "#2f6db0", "#9bbbe0", "#d9742a", "#222222", "#7a7a7a"
plt.rcParams.update({"font.size": 10, "text.color": TEXT, "axes.labelcolor": TEXT,
                     "xtick.color": TEXT, "ytick.color": TEXT})


def de(x, nk=1):
    return f"{x:,.{nk}f}".replace(",", "X").replace(".", ",").replace("X", ".")


fp = json.load(open(f"logs/{SAMPLE}.fastp.json"))
fs = {}
for line in open(f"qc/{SAMPLE}.flagstat.txt"):
    m = re.match(r"(\d+) \+ \d+ (.+)", line.strip())
    if m:
        fs[re.sub(r"\s*\(.*\)", "", m.group(2))] = int(m.group(1))

roh = fp["summary"]["before_filtering"]["total_reads"]
nach = fp["summary"]["after_filtering"]["total_reads"]
gemappt = fs["primary mapped"]
gepaart = fs["properly paired"]
dup = fs["primary duplicates"]
q30_v = fp["summary"]["before_filtering"]["q30_rate"] * 100
q30_n = fp["summary"]["after_filtering"]["q30_rate"] * 100

fig, (a, b) = plt.subplots(1, 2, figsize=(13, 4.8), gridspec_kw=dict(width_ratios=[1.35, 1]))

# a: Anteil der Rohdaten-Reads pro Schritt (Punktdiagramm, Achse beginnt bewusst bei 94 %)
namen = ["Rohdaten", "Nach fastp", "Primär gemappt", "Korrekt gepaart"]
werte = [roh, nach, gemappt, gepaart]
ys = list(range(len(namen)))[::-1]
anteile = [100 * w / roh for w in werte]
a.hlines(ys, 94, anteile, color=BLAU_HELL, lw=2)
a.scatter(anteile, ys, s=90, color=BLAU, zorder=3)
a.set_yticks(ys)
a.set_yticklabels(namen)
for y, w, p in zip(ys, werte, anteile):
    a.text(p + 0.25, y, f"{de(w / 1e9, 3)} Mrd. Reads  ·  {de(p)} %", va="center", fontsize=9)
a.set_xlim(94, 103.5)
a.set_xticks([94, 96, 98, 100])
a.set_xlabel("Anteil der Rohdaten-Reads (%)  –  Achse beginnt bei 94 %")
a.set_title("a  Verbleibende Reads pro Schritt", loc="left", fontweight="bold", fontsize=11)

# b: Verluste und Markierungen
n2 = ["Niedrige Qualität\n(verworfen)", "Zu viele N\n(verworfen)", "Adapter gekürzt\n(behalten)", "Duplikate\n(markiert, behalten)"]
w2 = [fp["filtering_result"]["low_quality_reads"], fp["filtering_result"]["too_many_N_reads"],
      fp["adapter_cutting"]["adapter_trimmed_reads"], dup]
f2 = [AKZENT, AKZENT, BLAU_HELL, BLAU_HELL]
ys2 = list(range(len(n2)))[::-1]
b.barh(ys2, [w / 1e6 for w in w2], color=f2, height=0.6)
b.set_yticks(ys2)
b.set_yticklabels(n2, fontsize=9)
for y, w in zip(ys2, w2):
    b.text(w / 1e6 + 0.8, y, f"{de(w / 1e6)} Mio.  ({de(100 * w / roh, 2)} % der Rohdaten)", va="center", fontsize=9)
b.set_xlim(0, 75)
b.set_xlabel("Reads (Millionen)")
b.set_title("b  Entfernt oder markiert", loc="left", fontweight="bold", fontsize=11)
a.text(0.99, 0.22, f"Q30 (Basenqualität): {de(q30_v, 2)} % → {de(q30_n, 2)} %", transform=a.transAxes, ha="right",
       fontsize=10, fontweight="bold", bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="#cccccc"))

for ax in (a, b):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.grid(True, axis="x", color="#e6e6e6", linewidth=0.8)
    ax.set_axisbelow(True)
fig.suptitle(f"Sequenzierung und Mapping von {SAMPLE} ({NAME})", x=0.01, ha="left", fontsize=13, fontweight="bold")
fig.tight_layout(rect=(0, 0, 1, 0.94))
fig.savefig(f"{OUT}/10_sequenzierung_mapping.png", dpi=200, facecolor="white")
print("Fertig: abbildungen/10_sequenzierung_mapping.png")
