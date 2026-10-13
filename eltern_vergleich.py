import csv

def lade(pfad):
    with open(pfad, encoding="utf-8") as f:
        return list(csv.DictReader(f, delimiter="\t"))

mutter = lade("varianten/ERR194147.kandidaten_streng.tsv")
vater = lade("varianten/ERR194146.kandidaten_streng.tsv")

def schluessel(z):
    return (z["CHROM"], z["POS"], z["REF"], z["ALT"])

def nach_gen(zeilen):
    d = {}
    for z in zeilen:
        d.setdefault(z["GEN"], []).append(z)
    return d

gm, gv = nach_gen(mutter), nach_gen(vater)
gemeinsam = sorted(set(gm) & set(gv))
nur_m = sorted(set(gm) - set(gv))
nur_v = sorted(set(gv) - set(gm))

ausgabe = []
def zeile(kat, gen, z, wer):
    hinweis = "X/Y: haploid beim Vater, Genotyp prüfen" if z["CHROM"] in ("X", "Y") else ""
    ausgabe.append([kat, gen, wer, z["CHROM"], z["POS"], z["REF"], z["ALT"],
                    z["GENOTYP"], z["KONSEQUENZ"], z["CLINVAR"], hinweis])

print(f"Kandidaten: Mutter {len(mutter)}, Vater {len(vater)}")
print(f"\nGene bei BEIDEN Eltern ({len(gemeinsam)}): mögliches 25 %-Risiko, falls rezessiv")
for g in gemeinsam:
    km = {schluessel(z) for z in gm[g]}
    kv = {schluessel(z) for z in gv[g]}
    art = "gleiche Variante" if km & kv else "verschiedene Varianten (compound-heterozygot möglich)"
    print(f"  {g}: {art}")
    for z in gm[g]: zeile("beide:" + art, g, z, "Mutter")
    for z in gv[g]: zeile("beide:" + art, g, z, "Vater")

print(f"\nNur Mutter ({len(nur_m)}): {', '.join(nur_m)}")
for g in nur_m:
    for z in gm[g]: zeile("nur_Mutter", g, z, "Mutter")
print(f"\nNur Vater ({len(nur_v)}): {', '.join(nur_v)}")
for g in nur_v:
    for z in gv[g]: zeile("nur_Vater", g, z, "Vater")

with open("varianten/eltern_vergleich.tsv", "w", encoding="utf-8", newline="") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(["KATEGORIE", "GEN", "ELTERNTEIL", "CHROM", "POS", "REF", "ALT",
                "GENOTYP", "KONSEQUENZ", "CLINVAR", "HINWEIS"])
    w.writerows(ausgabe)
print("\nGeschrieben: varianten/eltern_vergleich.tsv")
