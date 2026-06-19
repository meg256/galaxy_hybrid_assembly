"""
Hybrid Bacterial Genome Assembly Pipeline
==========================================
FastQC -> Unicycler -> RagTag -> Pilon
QUAST + BUSCO QC + PROKKA after every assembly-modifying stage (draft, scaffold,
pilon). Bandage graph visualization after Unicycler only (RagTag/Pilon
both emit FASTA, not a graph).
Downstream: MLST + ABRicate typing on the final
Pilon-polished assembly. MultiQC aggregates everything at the end.

DRAFT STATUS - frst attempt
"""
import os

configfile: "config.yaml"

# -------------------------------------------------------------------------
# Environment paths from config
# -------------------------------------------------------------------------
ASSEMBLY_ENV   = config["assembly_env_path"]
ANNOTATION_ENV = config["annotation_env_path"]
BWASAM_ENV     = config["bwasam_env_path"]

# -------------------------------------------------------------------------
# Convenience handles into the config
# -------------------------------------------------------------------------
SHORT_READS_DIR = config["short_reads_dir"]
LONG_READS_DIR  = config["long_reads_dir"]
REFERENCE       = config["reference_genome"]
OUTDIR          = config["output_dir"]
THREADS         = config["threads"]

# -------------------------------------------------------------------------
# Sample discovery
# -------------------------------------------------------------------------
SAMPLES = glob_wildcards(os.path.join(SHORT_READS_DIR, "{sample}_1.fastq.gz")).sample
STAGES = ["draft", "scaffold", "pilon"]

# =============================================================================
# Target rule
# =============================================================================
rule all:
    input:
        expand(os.path.join(OUTDIR, "{sample}/qc/quast/{stage}/report.tsv"),
               sample=SAMPLES, stage=STAGES),
        expand(os.path.join(OUTDIR, "{sample}/qc/busco/{stage}/short_summary.txt"),
               sample=SAMPLES, stage=STAGES),
        expand(os.path.join(OUTDIR, "{sample}/qc/bandage/{sample}_assembly_graph.jpg"),
               sample=SAMPLES),
        expand(os.path.join(OUTDIR, "{sample}/prokka_draft/{sample}.gff"),
               sample=SAMPLES),
        expand(os.path.join(OUTDIR, "{sample}/prokka_scaffold/{sample}.gff"),
               sample=SAMPLES),
        expand(os.path.join(OUTDIR, "{sample}/prokka/{sample}.gff"),
               sample=SAMPLES),
        expand(os.path.join(OUTDIR, "{sample}/mlst/{sample}_mlst.tsv"),
               sample=SAMPLES),
        expand(os.path.join(OUTDIR, "{sample}/abricate/{sample}_abricate.tsv"),
               sample=SAMPLES),
        expand(os.path.join(OUTDIR, "{sample}/fastqc/{sample}_1_fastqc.html"),
               sample=SAMPLES),
        expand(os.path.join(OUTDIR, "{sample}/fastqc/{sample}_2_fastqc.html"),
               sample=SAMPLES),
        os.path.join(OUTDIR, "multiqc/multiqc_report.html")

# =============================================================================
# FASTQC
# =============================================================================
rule fastqc:
    conda: ASSEMBLY_ENV
    input:
        r1 = os.path.join(SHORT_READS_DIR, "{sample}_1.fastq.gz"),
        r2 = os.path.join(SHORT_READS_DIR, "{sample}_2.fastq.gz"),
    output:
        html1 = os.path.join(OUTDIR, "{sample}/fastqc/{sample}_1_fastqc.html"),
        html2 = os.path.join(OUTDIR, "{sample}/fastqc/{sample}_2_fastqc.html"),
    params:
        outdir = lambda wc: os.path.join(OUTDIR, wc.sample, "fastqc"),
    threads: THREADS
    shell:
        """
        mkdir -p {params.outdir}
        fastqc {input.r1} {input.r2} --outdir {params.outdir} --threads {threads}
        """

# =============================================================================
# UNICYCLER
# =============================================================================
rule unicycler:
    conda: ASSEMBLY_ENV
    input:
        r1 = os.path.join(SHORT_READS_DIR, "{sample}_1.fastq.gz"),
        r2 = os.path.join(SHORT_READS_DIR, "{sample}_2.fastq.gz"),
        long = os.path.join(LONG_READS_DIR, "{sample}.fastq.gz"),
    output:
        fasta = os.path.join(OUTDIR, "{sample}/unicycler/assembly.fasta"),
        gfa   = os.path.join(OUTDIR, "{sample}/unicycler/assembly.gfa"),
    params:
        outdir            = lambda wc: os.path.join(OUTDIR, wc.sample, "unicycler"),
        mode              = config["unicycler"]["mode"],
        min_fasta_length  = config["unicycler"]["min_fasta_length"],
        linear_seqs       = config["unicycler"]["linear_seqs"],
        min_kmer_frac     = config["unicycler"]["spades"]["min_kmer_frac"],
        max_kmer_frac     = config["unicycler"]["spades"]["max_kmer_frac"],
        kmer_count        = config["unicycler"]["spades"]["kmer_count"],
        depth_filter      = config["unicycler"]["spades"]["depth_filter"],
        start_gene_id     = config["unicycler"]["rotation"]["start_gene_id"],
        start_gene_cov    = config["unicycler"]["rotation"]["start_gene_cov"],
        min_component_size = config["unicycler"]["graph_cleaning"]["min_component_size"],
        min_dead_end_size   = config["unicycler"]["graph_cleaning"]["min_dead_end_size"],
        scores            = config["unicycler"]["long_read_alignment"]["scores"],
        keep              = config["unicycler"]["keep"],
    threads: THREADS
    shell:
        """
        mkdir -p {params.outdir}
        unicycler \
            -1 {input.r1} -2 {input.r2} -l {input.long} \
            -o {params.outdir} \
            -t {threads} \
            --mode {params.mode} \
            --min_fasta_length {params.min_fasta_length} \
            --linear_seqs {params.linear_seqs} \
            --min_kmer_frac {params.min_kmer_frac} \
            --max_kmer_frac {params.max_kmer_frac} \
            --kmer_count {params.kmer_count} \
            --depth_filter {params.depth_filter} \
            --start_gene_id {params.start_gene_id} \
            --start_gene_cov {params.start_gene_cov} \
            --min_component_size {params.min_component_size} \
            --min_dead_end_size {params.min_dead_end_size} \
            --scores {params.scores} \
            --keep {params.keep}
        """

# =============================================================================
# QUAST + BUSCO + BANDAGE (draft)
# =============================================================================
rule quast_draft:
    conda: ASSEMBLY_ENV
    input:
        fasta = os.path.join(OUTDIR, "{sample}/unicycler/assembly.fasta"),
    output:
        report = os.path.join(OUTDIR, "{sample}/qc/quast/draft/report.tsv"),
    params:
        outdir = lambda wc: os.path.join(OUTDIR, wc.sample, "qc/quast/draft"),
        min_contig = config["quast"]["min_contig"],
        gene_finding_flag = "--gene-finding" if config["quast"]["gene_finding"] else "",
    threads: THREADS
    shell:
        """
        quast.py {input.fasta} \
            --min-contig {params.min_contig} \
            {params.gene_finding_flag} \
            --threads {threads} \
            -o {params.outdir}
        """

rule busco_draft:
    conda: ASSEMBLY_ENV
    input:
        fasta = os.path.join(OUTDIR, "{sample}/unicycler/assembly.fasta"),
    output:
        summary = os.path.join(OUTDIR, "{sample}/qc/busco/draft/short_summary.txt"),
    params:
        out_path = lambda wc: os.path.join(OUTDIR, wc.sample, "qc/busco"),
        out_name = "draft",
        lineage  = config["busco"]["lineage"],
        mode     = config["busco"]["mode"],
    threads: THREADS
    shell:
        """
        busco -i {input.fasta} \
            -l {params.lineage} \
            -m {params.mode} \
            -o {params.out_name} \
            --out_path {params.out_path} \
            -c {threads} \
            -f
        cp {params.out_path}/{params.out_name}/run_{params.lineage}/short_summary.txt {output.summary}
        """

rule bandage_draft:
    conda: ASSEMBLY_ENV
    input:
        gfa = os.path.join(OUTDIR, "{sample}/unicycler/assembly.gfa"),
    output:
        image = os.path.join(OUTDIR, "{sample}/qc/bandage/{sample}_assembly_graph.jpg"),
        info  = os.path.join(OUTDIR, "{sample}/qc/bandage/{sample}_assembly_graph_info.txt"),
    params:
        height = config["bandage"]["height"],
    shell:
        """
        mkdir -p $(dirname {output.image})
        Bandage image {input.gfa} {output.image} --height {params.height}
        Bandage info {input.gfa} > {output.info}
        """

rule prokka_draft:
    conda: ANNOTATION_ENV
    input:
        fasta = os.path.join(OUTDIR, "{sample}/unicycler/assembly.fasta"),
    output:
        gff = os.path.join(OUTDIR, "{sample}/prokka_draft/{sample}.gff"),
        fna = os.path.join(OUTDIR, "{sample}/prokka_draft/{sample}.fna"),
        log = os.path.join(OUTDIR, "{sample}/prokka_draft/prokka.log"),
        txt = os.path.join(OUTDIR, "{sample}/prokka_draft/{sample}.txt"),
    params:
        outdir = lambda wc: os.path.join(OUTDIR, wc.sample, "prokka_draft"),
        prefix = "{sample}_draft",
        kingdom = config["prokka"]["kingdom"],
        gcode = config["prokka"]["gcode"],
        mincontig = config["prokka"]["mincontig"],
        evalue = config["prokka"]["evalue"],
        gffver = config["prokka"]["gffver"],
    threads: THREADS
    shell:
        """
        mkdir -p {params.outdir}
        prokka \
            --kingdom {params.kingdom} \
            --gcode {params.gcode} \
            --mincontig {params.mincontig} \
            --evalue {params.evalue} \
            --gffver {params.gffver} \
            --outdir {params.outdir} \
            --prefix {params.prefix} \
            --force \
            --cpus {threads} \
            {input.fasta} \
            > {output.log} 2>&1

        # Prokka writes summary to PROKKA_*.txt inside outdir
        cp {params.outdir}/{params.prefix}.txt {output.txt}
        """

# =============================================================================
# RAGTAG
# =============================================================================
rule ragtag:
    conda: ASSEMBLY_ENV
    input:
        fasta = os.path.join(OUTDIR, "{sample}/unicycler/assembly.fasta"),
        ref   = REFERENCE,
    output:
        fasta = os.path.join(OUTDIR, "{sample}/ragtag/ragtag.scaffold.fasta"),
    params:
        outdir = lambda wc: os.path.join(OUTDIR, wc.sample, "ragtag"),
        extra_args = config["ragtag"]["extra_args"],
    threads: THREADS
    shell:
        """
        mkdir -p {params.outdir}
        ragtag.py scaffold {input.ref} {input.fasta} \
            -o {params.outdir} \
            -t {threads} \
            {params.extra_args}
        """

# =============================================================================
# QUAST + BUSCO (scaffold)
# =============================================================================
rule quast_scaffold:
    conda: ASSEMBLY_ENV
    input:
        fasta = os.path.join(OUTDIR, "{sample}/ragtag/ragtag.scaffold.fasta"),
    output:
        report = os.path.join(OUTDIR, "{sample}/qc/quast/scaffold/report.tsv"),
    params:
        outdir = lambda wc: os.path.join(OUTDIR, wc.sample, "qc/quast/scaffold"),
        min_contig = config["quast"]["min_contig"],
        gene_finding_flag = "--gene-finding" if config["quast"]["gene_finding"] else "",
    threads: THREADS
    shell:
        """
        quast.py {input.fasta} \
            --min-contig {params.min_contig} \
            {params.gene_finding_flag} \
            --threads {threads} \
            -o {params.outdir}
        """

rule busco_scaffold:
    conda: ASSEMBLY_ENV
    input:
        fasta = os.path.join(OUTDIR, "{sample}/ragtag/ragtag.scaffold.fasta"),
    output:
        summary = os.path.join(OUTDIR, "{sample}/qc/busco/scaffold/short_summary.txt"),
    params:
        out_path = lambda wc: os.path.join(OUTDIR, wc.sample, "qc/busco"),
        out_name = "scaffold",
        lineage  = config["busco"]["lineage"],
        mode     = config["busco"]["mode"],
    threads: THREADS
    shell:
        """
        busco -i {input.fasta} \
            -l {params.lineage} \
            -m {params.mode} \
            -o {params.out_name} \
            --out_path {params.out_path} \
            -c {threads} \
            -f
        cp {params.out_path}/{params.out_name}/run_{params.lineage}/short_summary.txt {output.summary}
        """

rule prokka_scaffold:
    conda: ANNOTATION_ENV
    input:
        fasta = os.path.join(OUTDIR, "{sample}/ragtag/ragtag.scaffold.fasta"),
    output:
        gff = os.path.join(OUTDIR, "{sample}/prokka_scaffold/{sample}.gff"),
        fna = os.path.join(OUTDIR, "{sample}/prokka_scaffold/{sample}.fna"),
        log = os.path.join(OUTDIR, "{sample}/prokka_scaffold/prokka.log"),
        txt = os.path.join(OUTDIR, "{sample}/prokka_scaffold/{sample}.txt"),
    params:
        outdir = lambda wc: os.path.join(OUTDIR, wc.sample, "prokka_scaffold"),
        prefix = "{sample}_scaffold",
        kingdom = config["prokka"]["kingdom"],
        gcode = config["prokka"]["gcode"],
        mincontig = config["prokka"]["mincontig"],
        evalue = config["prokka"]["evalue"],
        gffver = config["prokka"]["gffver"],
    threads: THREADS
    shell:
        """
        mkdir -p {params.outdir}
        prokka \
            --kingdom {params.kingdom} \
            --gcode {params.gcode} \
            --mincontig {params.mincontig} \
            --evalue {params.evalue} \
            --gffver {params.gffver} \
            --outdir {params.outdir} \
            --prefix {params.prefix} \
            --force \
            --cpus {threads} \
            {input.fasta} \
            > {output.log} 2>&1

        cp {params.outdir}/{params.prefix}.txt {output.txt}
        """

# =============================================================================
# PILON (bwa index, bwa mem, pilon)
# =============================================================================
rule bwa_index_for_pilon:
    conda: BWASAM_ENV
    input:
        fasta = os.path.join(OUTDIR, "{sample}/ragtag/ragtag.scaffold.fasta"),
    output:
        bwt = os.path.join(OUTDIR, "{sample}/ragtag/ragtag.scaffold.fasta.bwt"),
    shell:
        "bwa index {input.fasta}"

rule bwa_mem_for_pilon:
    conda: BWASAM_ENV
    input:
        fasta = os.path.join(OUTDIR, "{sample}/ragtag/ragtag.scaffold.fasta"),
        bwt   = os.path.join(OUTDIR, "{sample}/ragtag/ragtag.scaffold.fasta.bwt"),
        r1    = os.path.join(SHORT_READS_DIR, "{sample}_1.fastq.gz"),
        r2    = os.path.join(SHORT_READS_DIR, "{sample}_2.fastq.gz"),
    output:
        bam = os.path.join(OUTDIR, "{sample}/pilon/illumina.sorted.bam"),
    threads: THREADS
    shell:
        """
        bwa mem -t {threads} {input.fasta} {input.r1} {input.r2} \
            | samtools sort -@ {threads} -o {output.bam}
        samtools index {output.bam}
        """

rule pilon:
    conda: ANNOTATION_ENV
    input:
        fasta = os.path.join(OUTDIR, "{sample}/ragtag/ragtag.scaffold.fasta"),
        bam   = os.path.join(OUTDIR, "{sample}/pilon/illumina.sorted.bam"),
    output:
        fasta = os.path.join(OUTDIR, "{sample}/pilon/pilon_polished.fasta"),
    params:
        outdir   = lambda wc: os.path.join(OUTDIR, wc.sample, "pilon"),
        java_mem = config["pilon"]["java_mem"],
        fix      = config["pilon"]["fix"],
    threads: THREADS
    shell:
        """
        pilon -Xmx{params.java_mem} \
            --genome {input.fasta} \
            --frags {input.bam} \
            --output pilon_polished \
            --outdir {params.outdir} \
            --fix {params.fix} \
            --threads {threads} \
            --changes
        """

# =============================================================================
# QUAST + BUSCO (pilon)
# =============================================================================
rule quast_pilon:
    conda: ASSEMBLY_ENV
    input:
        fasta = os.path.join(OUTDIR, "{sample}/pilon/pilon_polished.fasta"),
    output:
        report = os.path.join(OUTDIR, "{sample}/qc/quast/pilon/report.tsv"),
    params:
        outdir = lambda wc: os.path.join(OUTDIR, wc.sample, "qc/quast/pilon"),
        min_contig = config["quast"]["min_contig"],
        gene_finding_flag = "--gene-finding" if config["quast"]["gene_finding"] else "",
    threads: THREADS
    shell:
        """
        quast.py {input.fasta} \
            --min-contig {params.min_contig} \
            {params.gene_finding_flag} \
            --threads {threads} \
            -o {params.outdir}
        """

rule busco_pilon:
    conda: ASSEMBLY_ENV
    input:
        fasta = os.path.join(OUTDIR, "{sample}/pilon/pilon_polished.fasta"),
    output:
        summary = os.path.join(OUTDIR, "{sample}/qc/busco/pilon/short_summary.txt"),
    params:
        out_path = lambda wc: os.path.join(OUTDIR, wc.sample, "qc/busco"),
        out_name = "pilon",
        lineage  = config["busco"]["lineage"],
        mode     = config["busco"]["mode"],
    threads: THREADS
    shell:
        """
        busco -i {input.fasta} \
            -l {params.lineage} \
            -m {params.mode} \
            -o {params.out_name} \
            --out_path {params.out_path} \
            -c {threads} \
            -f
        cp {params.out_path}/{params.out_name}/run_{params.lineage}/short_summary.txt {output.summary}
        """

# =============================================================================
# ANNOTATION (Prokka, MLST, Abricate)
# =============================================================================
rule prokka:
    conda: ANNOTATION_ENV
    input:
        fasta = os.path.join(OUTDIR, "{sample}/pilon/pilon_polished.fasta"),
    output:
        gff = os.path.join(OUTDIR, "{sample}/prokka/{sample}.gff"),
        fna = os.path.join(OUTDIR, "{sample}/prokka/{sample}.fna"),
        log = os.path.join(OUTDIR, "{sample}/prokka/prokka.log"),
        txt = os.path.join(OUTDIR, "{sample}/prokka/{sample}.txt"),
    params:
        outdir   = lambda wc: os.path.join(OUTDIR, wc.sample, "prokka"),
        prefix   = "{sample}",
        kingdom  = config["prokka"]["kingdom"],
        gcode    = config["prokka"]["gcode"],
        mincontig = config["prokka"]["mincontig"],
        evalue   = config["prokka"]["evalue"],
        gffver   = config["prokka"]["gffver"],
    threads: THREADS
    shell:
        """
        mkdir -p {params.outdir}
        prokka \
            --kingdom {params.kingdom} \
            --gcode {params.gcode} \
            --mincontig {params.mincontig} \
            --evalue {params.evalue} \
            --gffver {params.gffver} \
            --outdir {params.outdir} \
            --prefix {params.prefix} \
            --force \
            --cpus {threads} \
            {input.fasta} \
            > {output.log} 2>&1

        cp {params.outdir}/{params.prefix}.txt {output.txt}
        """

rule mlst:
    conda: ANNOTATION_ENV
    input:
        fna = os.path.join(OUTDIR, "{sample}/prokka/{sample}.fna"),
    output:
        tsv = os.path.join(OUTDIR, "{sample}/mlst/{sample}_mlst.tsv"),
    params:
        scheme_flag = lambda wc: (
            "" if config["mlst"]["scheme"] == ""
            else f"--scheme {config['mlst']['scheme']}"
        ),
    shell:
        """
        mkdir -p $(dirname {output.tsv})
        mlst {params.scheme_flag} {input.fna} > {output.tsv}
        """

rule abricate:
    conda: ANNOTATION_ENV
    input:
        fna = os.path.join(OUTDIR, "{sample}/prokka/{sample}.fna"),
    output:
        tsv = os.path.join(OUTDIR, "{sample}/abricate/{sample}_abricate.tsv"),
    params:
        db     = config["abricate"]["db"],
        minid  = config["abricate"]["minid"],
        mincov = config["abricate"]["mincov"],
    shell:
        """
        mkdir -p $(dirname {output.tsv})
        abricate --db {params.db} --minid {params.minid} --mincov {params.mincov} \
            {input.fna} > {output.tsv}
        """


rule multiqc:
    conda: ASSEMBLY_ENV
    input:
        # QUAST reports (all stages)
        expand(os.path.join(OUTDIR, "{sample}/qc/quast/{stage}/report.tsv"), sample=SAMPLES, stage=STAGES),
        # BUSCO summaries (all stages)
        expand(os.path.join(OUTDIR, "{sample}/qc/busco/{stage}/short_summary.txt"), sample=SAMPLES, stage=STAGES),
        # Bandage images (draft)
        expand(os.path.join(OUTDIR, "{sample}/qc/bandage/{sample}_assembly_graph.jpg"), sample=SAMPLES),
        # Prokka output txts (all stages) – ensures annotation is done
        expand(os.path.join(OUTDIR, "{sample}/prokka_draft/{sample}.txt"), sample=SAMPLES),
        expand(os.path.join(OUTDIR, "{sample}/prokka_scaffold/{sample}.txt"), sample=SAMPLES),
        expand(os.path.join(OUTDIR, "{sample}/prokka/{sample}.txt"), sample=SAMPLES),
        # FastQC HTMLs
        expand(os.path.join(OUTDIR, "{sample}/fastqc/{sample}_1_fastqc.html"), sample=SAMPLES),
        expand(os.path.join(OUTDIR, "{sample}/fastqc/{sample}_2_fastqc.html"), sample=SAMPLES),
    output:
        html = os.path.join(OUTDIR, "multiqc/multiqc_report.html"),
    params:
        outdir = os.path.join(OUTDIR, "multiqc"),
        search_root = OUTDIR,
    shell:
        """
        mkdir -p {params.outdir}
        multiqc {params.search_root} --outdir {params.outdir} --force
        """
