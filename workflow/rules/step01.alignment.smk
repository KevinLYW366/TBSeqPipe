##### Step 01.alignment #####

# Map sequencing reads to MTB reference genome by bwa mem
## Note: this step was done by bowtie2 before until extremely low alignment rate (bowtie2 < 1% vs. bwa men 5-10%)
##       was found for some contaminated samples.
rule bwa_map:
    input:
        fq1 = lambda wildcards: get_fastq_reads(wildcards.sample)[0],
        fq2 = lambda wildcards: get_fastq_reads(wildcards.sample)[1]
    output:
        "results/01.alignment/{sample}/{sample}.bam"
    params:
        ## bwa mem params
        # -R STR: read group header line such as '@RG\tID:foo\tSM:bar' [null]
        rg="'@RG\\tID:{sample}\\tSM:{sample}\\tPL:illumina'",
        genome=config["ref_genome"],
        # -M: mark shorter split hits as secondary
        short_split="-M",
        prefix="results/01.alignment/{sample}/{sample}"
    threads:
        config["threads"]
    log:
        "logs/01.alignment/bwa_map/{sample}.log"
    conda:
        "../envs/MTBAnalysis.yaml"
    shell:
        """
        (bwa mem -t {threads} -R {params.rg} {params.short_split} {params.genome} {input.fq1} {input.fq2} | \
        samtools view -@ {threads} -bSh | \
        samtools sort -@ {threads} -T {params.prefix} -O bam -o {output}) 2> {log}
        """

# generate quality control report of bam files by qualimap
rule qualimap:
    input:
        "results/01.alignment/{sample}/{sample}.bam"
    output:
        "results/01.alignment/qualimap/{sample}/qualimapReport.html"
    params:
        outdir = "results/01.alignment/qualimap/{sample}"
    threads:
        config["threads"]
    log:
        "logs/01.alignment/qualimap/{sample}.log"
    conda:
        "../envs/qc.yaml"
    shell:
        """
        # generate bam qc report using qualimap
        qualimap bamqc -nt {threads} -outdir {params.outdir} -bam {input} > {log} 2>&1
        """

# merge all qualimap reports using multiqc
rule multiqc_qualimap:
    input:
        lambda wildcards: get_qualified_results("results/01.alignment/qualimap/{sample}/qualimapReport.html", wildcards)
    output:
        "results/01.alignment/multiqc/multiqc_report.html"
    log:
        "logs/01.alignment/multiqc.log",
    params:
        indir = "results/01.alignment/qualimap",
        outdir = "results/01.alignment/multiqc/"
    conda:
        "../envs/qc.yaml"
    shell:
        """
        multiqc {params.indir} -m qualimap -f -o {params.outdir} > {log} 2>&1
        """

# Duplicate reads are marked and dedupliacted by SAMtools
rule samtools_dedup_sort:
    input:
        "results/01.alignment/{sample}/{sample}.bam"
    output:
        bam = temp("results/01.alignment/{sample}/{sample}_sd.bam"),
        bai = "results/01.alignment/{sample}/{sample}_sd.bam.bai"
    params:
        prefix="results/01.alignment/{sample}/{sample}_sd"
    threads:
        config["threads"]
    log:
        "logs/01.alignment/samtools_sort_index/{sample}.log"
    conda:
        "../envs/MTBAnalysis.yaml"
    shell:
        """
        # remove PCR duplicates
        (samtools sort -m 768M -n -@ {threads} {input} | \
        samtools fixmate -@ {threads} -m - - | samtools sort -m 768M -@ {threads} - | \
        samtools markdup -@ {threads} - {output.bam}) 2> {log}
        # index
        samtools index -@ {threads} -b {output.bam} {output.bai} >> {log} 2>&1
        """

# realign around InDels using GATK3
## first generate the FASTA sequence dictionary file for reference genome if not done, using following command:
## picard CreateSequenceDictionary -R ${ref_genome}
rule gatk_realinger:
    input:
        "results/01.alignment/{sample}/{sample}_sd.bam"
    output:
        temp("results/01.alignment/{sample}/{sample}_sdr.bam")
    log:
        "logs/01.alignment/gatk_realinger/{sample}.log",
    params:
        genome = config["ref_genome"],
        outdir = "results/01.alignment/{sample}"
    threads:
        config["threads"]
    conda:
        "../envs/gatk3.yaml"
    shell:
        """
        export _JAVA_OPTIONS=-Xmx32g > {log} 2>&1
        gatk3 -T RealignerTargetCreator -I {input} -R {params.genome} -o {params.outdir}/GATK.intervals \
        -nt {threads} >> {log} 2>&1  
        
        export _JAVA_OPTIONS=-Xmx4g >> {log} 2>&1
        gatk3 -T IndelRealigner -I {input} -R {params.genome} -targetIntervals {params.outdir}/GATK.intervals \
        -o {output} >> {log} 2>&1
        
        rm -f {params.outdir}/GATK.intervals >> {log} 2>&1
        """

# base quality recalibrator using GATK3
rule gatk_base_recalibrator:
    input:
        "results/01.alignment/{sample}/{sample}_sdr.bam"
    output:
        temp("results/01.alignment/{sample}/{sample}_sdrb.bam")
    log:
        "logs/01.alignment/gatk_base_recalibrator/{sample}.log",
    params:
        genome = config["ref_genome"],
        known_site = config["gatk_base_recalibrator_vcf"],
        outdir = "results/01.alignment/{sample}"
    threads:
        config["threads"]
    conda:
        "../envs/gatk3.yaml"
    shell:
        """
        export _JAVA_OPTIONS=-Xmx4g > {log} 2>&1
        gatk3 -T BaseRecalibrator -I {input} -R {params.genome} --knownSites {params.known_site} \
        -o {params.outdir}/GATK_Resilist.grp -nct {threads} >> {log} 2>&1
        
        export _JAVA_OPTIONS=-Xmx4g >> {log} 2>&1
        gatk3 -T PrintReads -I {input} -R {params.genome} -BQSR {params.outdir}/GATK_Resilist.grp \
        -o {output} -nct {threads} >> {log} 2>&1
        
        rm -f {params.outdir}/GATK_Resilist.grp >> {log} 2>&1
        """

# sort and index bam file using samtools
rule samtools_sort_index:
    input:
        "results/01.alignment/{sample}/{sample}_sdrb.bam"
    output:
        bam = temp("results/01.alignment/{sample}/{sample}_sdrbs.bam"),
        bai = "results/01.alignment/{sample}/{sample}_sdrbs.bam.bai"
    log:
        "logs/01.alignment/samtools_sort_index/{sample}.log",
    params:
        prefix = "results/01.alignment/{sample}/{sample}_sdrbs"
    threads:
        config["threads"]
    conda:
        "../envs/MTBAnalysis.yaml"
    shell:
        """
        samtools sort -@ {threads} -T {params.prefix} -O bam {input} > {output.bam} 2> {log}
        samtools index -@ {threads} -b {output.bam} {output.bai} >> {log} 2>&1
        """

# filter out unmapped reads using samtools view
rule samtools_remove_unmapped:
    input:
        "results/01.alignment/{sample}/{sample}_sdrbs.bam"
    output:
        bam = "results/01.alignment/{sample}/{sample}_final.bam",
        bai = "results/01.alignment/{sample}/{sample}_final.bam.bai"
    log:
        "logs/01.alignment/samtools_remove_unmapped/{sample}.log",
    params:
        # -b       output BAM
        # -h       include header in SAM output
        # -F INT   only include reads with none of the FLAGS in INT present [0]
        ## FLAGS 4: Read unmapped
        samtools_view_params = "-bhF 4"
    threads:
        config["threads"]
    conda:
        "../envs/MTBAnalysis.yaml"
    shell:
        """
        # filter out unmapped reads
        samtools view {params.samtools_view_params} -o {output.bam} {input} 2> {log}
        # index
        samtools index -@ {threads} -b {output.bam} {output.bai} >> {log} 2>&1
        # remove tmp files
        rm -f results/01.alignment/{wildcards.sample}/{wildcards.sample}_sd*.bai >> {log} 2>&1
        """