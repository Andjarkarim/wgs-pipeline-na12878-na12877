import re
import sys

STRENG = re.compile(r"^(Pathogenic|Likely_pathogenic)")
PROTEIN = re.compile(r"missense|stop_gained|frameshift|splice_acceptor|splice_donor|start_lost|stop_lost|inframe")

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

gesehen = set()
verworfen = []
with open(sys.argv[1]) as f:
    print(f.readline().rstrip("\n") + "\tQUALITAET")
    for line in f:
        c = line.rstrip("\n").split("\t")
        chrom, pos, ref, alt, gen, kons, clin, gt, ad, dp, gq, ab = c[:12]
        if STRENG.match(clin) and PROTEIN.search(kons) and "NMD" not in kons:
            schluessel = (chrom, pos, ref, alt, gen)
            if schluessel in gesehen:
                continue
            gesehen.add(schluessel)
            q = qualitaet(gt, ab, dp, gq, ref, alt)
            if q == "artefaktverdacht":
                verworfen.append(f"  {gen:9s} {chrom}:{pos} AB={ab} DP={dp} GQ={gq}")
            else:
                print("\t".join(c[:12]) + "\t" + q)
print(f"\nVerworfen (Artefaktverdacht, AB < 0.25): {len(verworfen)}", file=sys.stderr)
for v in verworfen:
    print(v, file=sys.stderr)
