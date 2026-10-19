import sys
from collections import Counter
fields = None
seen = set()
zaehler = Counter()
print("CHROM\tPOS\tREF\tALT\tGEN\tKONSEQUENZ\tCLINVAR\tGENOTYP\tAD\tDP\tGQ\tAB")
for line in open(sys.argv[1]):
    if line.startswith("##INFO=<ID=CSQ"):
        fields = line.split("Format: ")[1].rstrip('">\n').split("|")
    elif line.startswith("#"):
        continue
    else:
        c = line.rstrip("\n").split("\t")
        info = dict(x.split("=", 1) for x in c[7].split(";") if "=" in x)
        gt, ad, dp, gq, ab = "", ".", ".", ".", "NA"
        if len(c) > 9:
            fmt = dict(zip(c[8].split(":"), c[9].split(":")))
            gt = fmt.get("GT", "")
            ad = fmt.get("AD", ".")
            dp = fmt.get("DP", ".")
            gq = fmt.get("GQ", ".")
            try:
                n = [int(x) for x in ad.split(",")]
                if sum(n) > 0:
                    ab = f"{sum(n[1:]) / sum(n):.2f}"
            except ValueError:
                pass
        for e in info.get("CSQ", "").split(","):
            d = dict(zip(fields, e.split("|")))
            cs = d.get("ClinVar_CLNSIG", "")
            if "athogenic" in cs:
                key = (c[0], c[1], c[3], c[4], d.get("SYMBOL", ""), d.get("Consequence", ""), cs, gt, ad, dp, gq, ab)
                if key not in seen:
                    seen.add(key)
                    zaehler[cs] += 1
                    print("\t".join(key))
print("\nZusammenfassung:", file=sys.stderr)
for k, v in zaehler.most_common():
    print(f"  {v:4d}  {k}", file=sys.stderr)
