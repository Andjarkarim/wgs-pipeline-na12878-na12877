#!/usr/bin/env bash
# Validierung der Variantenbestimmung gegen die Platinum-Genomes-Referenzcalls (2017-1.0, GRCh38).
# Aufruf im Projektordner:  conda activate validierung && bash skripte/validierung.sh
# Die Ausgabeordner validierung/vcfeval_* dürfen vorher nicht existieren (rtg bricht sonst ab).
set -euo pipefail

B=https://platinum-genomes.s3.amazonaws.com/2017-1.0/hg38/small_variants
mkdir -p validierung/truth validierung/ergebnisse
cd validierung/truth

# 1. Truth Set und Konfidenzregionen laden (nur, was noch fehlt)
for f in ConfidentRegions.bed.gz ConfidentRegions.bed.gz.tbi \
         NA12877/NA12877.vcf.gz NA12877/NA12877.vcf.gz.tbi \
         NA12878/NA12878.vcf.gz NA12878/NA12878.vcf.gz.tbi; do
  [ -s "$(basename "$f")" ] || curl -fL --retry 3 -O "$B/$f"
done

# 2. Chromosomennamen angleichen: chr1 -> 1 (Schreibweise der Referenz); das Truth Set kennt nur 1-22 und X
for c in $(seq 1 22) X; do echo "chr$c $c"; done > rename_chr.txt
R=$(for c in $(seq 1 22) X; do printf 'chr%s,' "$c"; done | sed 's/,$//')
for s in NA12877 NA12878; do
  bcftools view -r "$R" "$s.vcf.gz" -Ou | bcftools annotate --rename-chrs rename_chr.txt -Oz -o "$s.ens.vcf.gz"
  bcftools index -t "$s.ens.vcf.gz"
done
gzip -dc ConfidentRegions.bed.gz | sed 's/^chr//' > ConfidentRegions.ens.bed
cd ..

# 3. Referenz ins rtg-Format umwandeln (einmalig)
[ -d GRCh38.sdf ] || rtg format -o GRCh38.sdf ../referenz/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

# 4. Vergleich: NA12878 = Mutter (ERR194147), NA12877 = Vater (ERR194146)
rtg vcfeval -b truth/NA12878.ens.vcf.gz -c ../varianten/ERR194147.vcf.gz -t GRCh38.sdf \
  --bed-regions truth/ConfidentRegions.ens.bed -o vcfeval_NA12878 --vcf-score-field QUAL --threads 4
rtg vcfeval -b truth/NA12877.ens.vcf.gz -c ../varianten/ERR194146.vcf.gz -t GRCh38.sdf \
  --bed-regions truth/ConfidentRegions.ens.bed -o vcfeval_NA12877 --vcf-score-field QUAL --threads 4

# 5. Ergebnisse zusammenfassen (die letzte Zeile der ROC-Dateien entspricht den ungefilterten Calls)
{
  echo "Spalten: Schwelle, TP (Truth Set), FP, TP (Calls), FN, Praezision, Sensitivitaet, F-Mass"
  for s in NA12878 NA12877; do
    echo; echo "== $s =="
    cat vcfeval_$s/summary.txt
    echo "-- SNVs (letzte Zeile von snp_roc) --"
    gzip -dc vcfeval_$s/snp_roc.tsv.gz | tail -1
    echo "-- Indels u.a. (letzte Zeile von non_snp_roc) --"
    gzip -dc vcfeval_$s/non_snp_roc.tsv.gz | tail -1
  done
} > ergebnisse/zusammenfassung.txt
cat ergebnisse/zusammenfassung.txt
