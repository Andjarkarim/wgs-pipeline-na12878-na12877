import json, sys

fastp_json, dedup_txt, flagstat_txt = sys.argv[1:4]

d = json.load(open(fastp_json))
s = d["summary"]
print("== fastp (Rohdaten -> getrimmt) ==")
for k, name in (("before_filtering", "vorher"), ("after_filtering", "nachher")):
    x = s[k]
    print(f"{name}: Reads {x['total_reads']:,}  Basen {x['total_bases']:,}  "
          f"Q30 {x['q30_rate'] * 100:.2f} %  GC {x['gc_content'] * 100:.1f} %  "
          f"Leselänge Ø {x['read1_mean_length']} bp")
print("Filterergebnis:", d.get("filtering_result"))
print("Adapter-Trimming:", d.get("adapter_cutting"))
print("Duplikationsrate (fastp-Schätzung):", d.get("duplication"))

print("\n== Duplikate (GATK MarkDuplicates) ==")
zeilen = open(dedup_txt).read().splitlines()
i = next(k for k, l in enumerate(zeilen) if l.startswith("## METRICS CLASS"))
m = dict(zip(zeilen[i + 1].split("\t"), zeilen[i + 2].split("\t")))
for k in ("READ_PAIRS_EXAMINED", "UNPAIRED_READS_EXAMINED", "READ_PAIR_DUPLICATES",
          "PERCENT_DUPLICATION", "ESTIMATED_LIBRARY_SIZE"):
    print(f"{k}: {m.get(k)}")

print("\n== Mapping (samtools flagstat, dedup.bam) ==")
print("".join(open(flagstat_txt).readlines()[:12]))
