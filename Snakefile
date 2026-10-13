SAMPLE = config.get("sample", "ERR194147")
NAME = config.get("name", "NA12878")
REF = "referenz/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
THREADS = 4

rule all:
    input:
        f"varianten/{SAMPLE}.kandidaten_streng.tsv"

rule fastqc:
    input:
        r1 = f"rohdaten/{SAMPLE}_1.fastq.gz",
        r2 = f"rohdaten/{SAMPLE}_2.fastq.gz"
    output:
        html1 = f"qc/{SAMPLE}_1_fastqc.html",
        html2 = f"qc/{SAMPLE}_2_fastqc.html"
    threads: THREADS
    shell:
        "fastqc {input.r1} {input.r2} -o qc/ -t {threads}"

rule trimming:
    input:
        r1 = f"rohdaten/{SAMPLE}_1.fastq.gz",
        r2 = f"rohdaten/{SAMPLE}_2.fastq.gz"
    output:
        r1 = temp(f"trimmed/{SAMPLE}_1.trim.fastq.gz"),
        r2 = temp(f"trimmed/{SAMPLE}_2.trim.fastq.gz"),
        json = f"logs/{SAMPLE}.fastp.json",
        html = f"logs/{SAMPLE}.fastp.html"
    threads: THREADS
    shell:
        """
        fastp -i {input.r1} -I {input.r2} \
              -o {output.r1} -O {output.r2} \
              -j {output.json} -h {output.html} \
              -w {threads}
        """

rule alignment:
    input:
        r1 = f"trimmed/{SAMPLE}_1.trim.fastq.gz",
        r2 = f"trimmed/{SAMPLE}_2.trim.fastq.gz",
        ref = REF
    output:
        bam = temp(f"alignment/{SAMPLE}.sorted.bam")
    threads: THREADS
    shell:
        """
        bwa mem -t {threads} -R '@RG\\tID:{SAMPLE}\\tSM:{NAME}\\tPL:ILLUMINA' \
            {input.ref} {input.r1} {input.r2} \
            | samtools sort -@ {threads} -o {output.bam} -
        samtools index {output.bam}
        """

rule mark_duplicates:
    input:
        bam = f"alignment/{SAMPLE}.sorted.bam"
    output:
        bam = f"alignment/{SAMPLE}.dedup.bam",
        metrics = f"logs/{SAMPLE}.dedup_metrics.txt"
    shell:
        """
        gatk MarkDuplicates \
            -I {input.bam} \
            -O {output.bam} \
            -M {output.metrics}
        samtools index {output.bam}
        """

rule variant_calling:
    input:
        bam = f"alignment/{SAMPLE}.dedup.bam",
        ref = REF
    output:
        vcf = f"varianten/{SAMPLE}.vcf.gz"
    shell:
        """
        gatk HaplotypeCaller \
            -R {input.ref} \
            -I {input.bam} \
            -O {output.vcf}
        """

rule vep_annotate:
    input:
        vcf = f"varianten/{SAMPLE}.vcf.gz",
        clinvar = "referenz/clinvar.vcf.gz"
    output:
        vcf = f"varianten/{SAMPLE}.vep.vcf.gz"
    threads: THREADS
    shell:
        """
        conda run -n vep --no-capture-output vep \
            -i {input.vcf} \
            --cache --offline --assembly GRCh38 --dir_cache ~/.vep \
            --custom file={input.clinvar},short_name=ClinVar,format=vcf,type=exact,coords=0,fields=CLNSIG \
            --vcf --compress_output bgzip --fork {threads} \
            -o {output.vcf}
        """

rule filter_pathogenic:
    input:
        vcf = f"varianten/{SAMPLE}.vep.vcf.gz"
    output:
        vcf = f"varianten/{SAMPLE}.pathogene.vcf"
    shell:
        """
        conda run -n vep --no-capture-output filter_vep \
            -i {input.vcf} \
            --filter "ClinVar_CLNSIG match athogenic" \
            -o {output.vcf}
        """

rule kandidaten_tabelle:
    input:
        vcf = f"varianten/{SAMPLE}.pathogene.vcf",
        skript = "auswertung.py"
    output:
        tsv = f"varianten/{SAMPLE}.pathogene_tabelle.tsv"
    shell:
        "python {input.skript} {input.vcf} > {output.tsv}"

rule kandidaten_streng:
    input:
        tsv = f"varianten/{SAMPLE}.pathogene_tabelle.tsv",
        skript = "kandidaten_filter.py"
    output:
        tsv = f"varianten/{SAMPLE}.kandidaten_streng.tsv"
    shell:
        "python {input.skript} {input.tsv} > {output.tsv}"

rule qc_zusammenfassung:
    input:
        fastp = f"logs/{SAMPLE}.fastp.json",
        dedup = f"logs/{SAMPLE}.dedup_metrics.txt",
        bam = f"alignment/{SAMPLE}.dedup.bam",
        skript = "qc_zusammenfassung.py"
    output:
        flagstat = f"qc/{SAMPLE}.flagstat.txt",
        txt = f"qc/{SAMPLE}.zusammenfassung.txt"
    threads: THREADS
    shell:
        """
        samtools flagstat -@ {threads} {input.bam} > {output.flagstat}
        python {input.skript} {input.fastp} {input.dedup} {output.flagstat} > {output.txt}
        """
