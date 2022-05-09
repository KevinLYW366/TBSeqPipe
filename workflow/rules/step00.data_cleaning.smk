##### Step 00.data_cleaning (OPTIONAL) #####

# clean raw fastq data
if config["data_cleaning"]:
    rule fastp_clean:
        input:
            fq1=config["data_dir"] + "%s/{sample}_%s1.%s" % (DATADIRINDEX, FASTQREAD, FASTQSUFFIX),
            fq2=config["data_dir"] + "%s/{sample}_%s2.%s" % (DATADIRINDEX, FASTQREAD, FASTQSUFFIX)
        output:
            cfq1 = "results/00.data_cleaning/{sample}/{sample}_1.clean.fastq.gz",
            cfq2 = "results/00.data_cleaning/{sample}/{sample}_2.clean.fastq.gz",
            json = "results/00.data_cleaning/{sample}/{sample}_fastp.json",
            html = "results/00.data_cleaning/{sample}/{sample}_fastp.html"
        threads:
            config["threads"]
        log:
            "logs/00.data_cleaning/fastp_clean/{sample}.log"
        conda:
            "../envs/MTBAnalysis.yaml"
        shell:
            "fastp --detect_adapter_for_pe -w {threads} -i {input.fq1} -I {input.fq2} "
            "-o {output.cfq1} -O {output.cfq2} -j {output.json} -h {output.html} 2> {log}"

# Use kraken to detect the percentage of reads mapping to MTBC
rule kraken_detect:
    input:
        fq1 = lambda wildcards: get_fastq_reads(wildcards.sample)[0],
        fq2 = lambda wildcards: get_fastq_reads(wildcards.sample)[1]
    output:
        reads = temp("results/00.data_cleaning/{sample}/{sample}.kraken"),
        report = "results/00.data_cleaning/{sample}/{sample}.kraken.report"
    log:
        "logs/00.data_cleaning/kraken_detect/{sample}.log"
    conda:
        "../envs/kraken.yaml"
    threads:
        config["threads"]
    params:
        kraken_db = config["kraken_db"]
    shell:
        """
        kraken --db {params.kraken_db} --threads {threads} --gzip-compressed --fastq-input \
        --paired {input.fq1} {input.fq2} > {output.reads} 2> {log}
        kraken-report --db {params.kraken_db} {output.reads} > {output.report} 2>> {log}  
        """

rule kraken_filter:
    input:
        "results/00.data_cleaning/{sample}/{sample}.kraken.report"
    output:
        "results/00.data_cleaning/{sample}/{sample}.kraken.result"
    log:
        "logs/00.data_cleaning/kraken_filter/{sample}.log"
    params:
        sample_name = "{sample}",
        script =  "workflow/scripts/kraken_report_extract.py",
    conda:
        "../envs/MTBAnalysis.yaml"
    shell:
        "python {params.script} {input} {output} {params.sample_name} 2> {log}"

if config["kraken_filter"]:
    checkpoint kraken_qc:
        input:
            expand("results/00.data_cleaning/{sample}/{sample}.kraken.result", sample=SAMPLES)
        output:
            qc_result = "results/00.data_cleaning/kraken.result",
            # Updated sample list which passed kraken qc
            pass_sample_list = "results/00.data_cleaning/kraken.samples"
        params:
            qc_threshold = config["kraken_threshold"]
        log:
            "logs/00.data_cleaning/kraken_result_merge.log"
        shell:
            """
            # Generate kraken qc result output
            echo -e "##The percentage of reads mapping to MTBC and the category detected by kraken" > {output.qc_result}
            echo -e "##Kraken QC threshold: {params.qc_threshold}% to MTBC (samples with lower percentage won't run the rest of workflow)" >> {output.qc_result}
            echo -e "#sample\tpercentage_MTBC\tcategory_detected" >> {output.qc_result}
            cat {input} >> {output.qc_result}
            TMP=`head -n 3 {output.qc_result} && grep -v "^#" {output.qc_result} | sort -k2n` 
            echo "$TMP" > {output.qc_result}
            # Extract passing sample list from kraken qc result output
            ## the braces used in awk have to be {{}} to distinguished from snakemake variable names
            grep -v '^#' {output.qc_result} | awk -F '\t' '$2 >= {params.qc_threshold} {{print $1}}' | sort -k1 > {output.pass_sample_list}
            """
else:
    rule kraken_qc:
        input:
            expand("results/00.data_cleaning/{sample}/{sample}.kraken.result",sample=SAMPLES)
        output:
            qc_result="results/00.data_cleaning/kraken.result",
            # Updated sample list which passed kraken qc
            pass_sample_list="results/00.data_cleaning/kraken.samples"
        log:
            "logs/00.data_cleaning/kraken_result_merge.log"
        shell:
            """
            # Generate kraken qc result output
            echo -e "##The percentage of reads mapping to MTBC and the category detected by kraken" > {output.qc_result}
            echo -e "#sample\tpercentage_MTBC\tcategory_detected" >> {output.qc_result}
            cat {input} >> {output.qc_result}
            TMP=`head -n 2 {output.qc_result} && grep -v "^#" {output.qc_result} | sort -k2n` 
            echo "$TMP" > {output.qc_result}
            # Extract passing sample list from kraken qc result output
            ## the braces used in awk have to be {{}} to distinguished from snakemake variable names
            grep -v '^#' {output.qc_result} | awk -F '\t' '$2 >= 0 {{print $1}}' | sort -k1 > {output.pass_sample_list}
            """

# generate a quality control report of fastq reads by fastqc
rule fastqc:
    input:
        fq1 = lambda wildcards: get_fastq_reads(wildcards.sample)[0],
        fq2 = lambda wildcards: get_fastq_reads(wildcards.sample)[1]
    output:
        fq1_html = "results/00.data_cleaning/fastqc/{sample}_1_fastqc.html",
        fq1_zip = "results/00.data_cleaning/fastqc/{sample}_1_fastqc.zip",
        fq2_html = "results/00.data_cleaning/fastqc/{sample}_2_fastqc.html",
        fq2_zip = "results/00.data_cleaning/fastqc/{sample}_2_fastqc.zip",
    log:
        "logs/00.data_cleaning/fastqc/{sample}.log"
    params:
        outdir = "results/00.data_cleaning/fastqc",
        fq1_prefix = lambda wildcards: get_fastq_reads(wildcards.sample)[0].split("/")[-1].replace(".fastq.gz", ""),
        fq2_prefix = lambda wildcards: get_fastq_reads(wildcards.sample)[1].split("/")[-1].replace(".fastq.gz", "")
    threads:
        config["threads"]
    conda:
        "../envs/qc.yaml"
    shell:
        """
        # run fastqc
        fastqc -t {threads} -o {params.outdir} {input.fq1} {input.fq2} > {log} 2>&1
        # modify fastqc result file names
        mv {params.outdir}/{params.fq1_prefix}_fastqc.html {params.outdir}/{wildcards.sample}_1_fastqc.html 2>> {log}
        mv {params.outdir}/{params.fq1_prefix}_fastqc.zip {params.outdir}/{wildcards.sample}_1_fastqc.zip 2>> {log}
        mv {params.outdir}/{params.fq2_prefix}_fastqc.html {params.outdir}/{wildcards.sample}_2_fastqc.html 2>> {log}
        mv {params.outdir}/{params.fq2_prefix}_fastqc.zip {params.outdir}/{wildcards.sample}_2_fastqc.zip 2>> {log}
        """

# merge all fastqc reports using multiqc
rule multiqc_fastqc:
    input:
        expand("results/00.data_cleaning/fastqc/{sample}_{read}_fastqc.html", sample=SAMPLES, read=[1,2])
    output:
        "results/00.data_cleaning/multiqc/multiqc_report.html"
    log:
        "logs/00.data_cleaning/multiqc.log"
    params:
        indir = "results/00.data_cleaning/fastqc",
        outdir = "results/00.data_cleaning/multiqc/"
    conda:
        "../envs/qc.yaml"
    shell:
        """
        multiqc {params.indir} -m fastqc -f -o {params.outdir} > {log} 2>&1
        """