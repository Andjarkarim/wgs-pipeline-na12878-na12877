# WGS-Variantenpipeline mit Snakemake: NA12878 und NA12877

Pipeline von FASTQ-Rohdaten zu annotierten ClinVar-Kandidaten für zwei Gesamtgenome des Platinum-Genomes-Pedigrees (NA12878, Mutter, ERR194147; NA12877, Vater, ERR194146) mit einem Vergleich der Eltern.

**Wichtig:** Alle Ergebnisse sind ClinVar-Annotationen auf öffentlichen Beispieldaten. Sie sind keine Diagnosen und keine Risikoabschätzung für reale Personen.

## Ablauf

FASTQ -> FastQC -> fastp (Trimmen) -> BWA-MEM (GRCh38) -> samtools sort -> GATK MarkDuplicates -> GATK HaplotypeCaller -> VEP (Release 116) mit ClinVar -> filter_vep -> Kandidatentabelle -> strenger Filter (Pathogenic / Likely pathogenic) -> Qualitätsklassen (Allelbalance, Tiefe, Genotypqualität) -> Elternvergleich.

## Ergebnisse

| | Mutter (NA12878) | Vater (NA12877) |
|---|---|---|
| Varianten (ungefiltert) | 4.893.226 | 4.889.622 |
| mit ClinVar-Eintrag | 44.476 | 43.890 |
| streng (P/LP) | 11 | 11 |
| nach Qualitätsfilter | 7 | 6 |
| davon belastbar | 3 (TPP1, PAH, PLG) | 4 (POLG, CCDC40, CFTR, CLCN1) |

Im Elternvergleich gibt es kein Gen mit einer belastbaren Variante bei beiden Eltern. SERPINA1 erscheint bei beiden mit derselben Variante, aber niedriger Allelbalance und ohne Bestätigung in den Roh-Alignments (vermutlich ein Artefakt). Details stehen in `bericht.tex` und `dokumentation.tex`.

## Aufruf

Mutter (Standardprobe) und Vater:

```
snakemake -n
snakemake --cores 4
snakemake --cores 4 --config sample=ERR194146 name=NA12877
```

Qualitätszusammenfassung, Abbildungen und Elternvergleich laufen außerhalb von Snakemake:

```
SAMPLE=ERR194146 NAME=NA12877 python3 abbildungen.py
SAMPLE=ERR194146 NAME=NA12877 python3 abbildungen2.py --neu
SAMPLE=ERR194146 NAME=NA12877 python3 einzelgenom.py
python3 abbildungen3.py
python3 eltern_vergleich.py
python3 abbildung_paar.py
```

## Software

Snakemake 6.1.2, fastp, BWA 0.7.x, samtools, GATK 4.6.x, bcftools, Ensembl VEP 116 (eigene conda-Umgebung), Python 3.10 mit matplotlib und pycirclize. Die Umgebung steht in `environment.yml`.

## Einschränkungen

- Keine Validierung gegen ein Truth Set (Platinum Genomes oder GIAB). Sensitivität und Präzision sind nicht quantifiziert.
- Das VCF ist ungefiltert (kein VQSR, kein Hard Filter, kein BQSR). Die Qualitätsklassen gelten nur für die ClinVar-Kandidaten, ihre Schwellen sind nicht kalibriert.
- Treffer in Regionen mit Paralogen, Pseudogenen oder hohem GC-Gehalt (zum Beispiel PRSS1, CNN2, PI4KA) sind anfällig für Mapping-Artefakte.
- Etwa 0,6 % der gerufenen Varianten fehlen in der VEP-Ausgabe; die Ursache ist nicht vollständig geklärt.
- ClinVar ändert sich laufend (hier Stand 2026-10-04). Krankheitsnamen und Prüfsterne sind wörtlich übernommen.

## Daten

Rohdaten: European Nucleotide Archive, ERR194147 und ERR194146. Sie liegen nicht im Repository. Referenz: Ensembl GRCh38 primary assembly (Release 111).

## Autor

Karim Andjar, Master Bioinformatik und Systembiologie, Justus-Liebig-Universität Gießen.

## Entstehung

Die Pipeline, die Auswertung und die Texte habe ich mit Unterstützung von Claude (Anthropic) entwickelt. Claude hat beim Schreiben der Skripte, beim Debuggen und bei der Formulierung der Texte geholfen. Die Auswahl der Daten, die Ausführung aller Schritte auf meinem Rechner, die Prüfung der Ergebnisse und die inhaltliche Verantwortung liegen bei mir.

## Lizenz

MIT
