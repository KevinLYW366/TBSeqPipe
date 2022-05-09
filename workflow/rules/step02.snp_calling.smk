##### Step 02.snp_calling #####
# Note: Variants would be called as diplotypes and two different types of SNPs calling (fSNP for phylogenetic analysis and vSNP for DR and mixed infection)

# SNP calling using bcftools mpileup and bcftools call
rule bcftools_call:
    input:
        "results/01.alignment/{sample}/{sample}_final.bam"
    output:
        "results/02.snp_calling/{sample}/{sample}_raw.vcf"
    params:
        ## bcftools mpileup params
        genome = config["ref_genome"],
        # -O, --output-type TYPE: 'b' compressed BCF; 'u' uncompressed BCF;
        #                         'z' compressed VCF; 'v' uncompressed VCF [v]
        output_type = "-Ou",
        # -q, --min-MQ INT: skip alignments with mapQ smaller than INT [0]
        minMQ = "-q 20",
        # -Q, --min-BQ INT: skip bases with baseQ/BAQ smaller than INT [1]
        minBQ = "-Q 20",
        # -a, --annotate LIST: optional tags to output; '?' to list available tags []
        annotate = "-a FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP,FORMAT/SP,INFO/AD,INFO/ADF,INFO/ADR",
        ## bcftools call params
        # --ploidy ASSEMBLY[?]: Predefined ploidy, 'list' to print available settings, append '?' for details [2]
        ploidy = "--ploidy 2",
        # -m, --multiallelic-caller: Alternative model for multiallelic and rare-variant calling (conflicts with -c)
        multiallelic = "-m",
        # -v, --variants-only: Output variant sites only
        variants_only = "-v",
        # -O, --output-type b|u|z|v: Output type: 'b' compressed BCF; 'u' uncompressed BCF;
        #                                         'z' compressed VCF; 'v' uncompressed VCF [v]
        output_type2 = "-Ov"
    threads:
        config["threads"]
    log:
        "logs/02.snp_calling/bcftools_call/{sample}.log"
    conda:
        "../envs/MTBAnalysis.yaml"
    shell:
        # bcftools mpileup: multi-way pileup producing genotype likelihoods
        "(bcftools mpileup --threads {threads} {params.output_type} "
        "{params.minMQ} {params.minBQ} -f {params.genome} {params.annotate} "
        "{input} | " 
        # bcftools call: SNP/indel calling (former "view")
        "bcftools call --threads {threads} {params.ploidy} "
        "{params.multiallelic} {params.variants_only} "
        "{params.output_type2} -o {output}) 2> {log}"

# filter SNP calling results using bcftools filter
rule bcftools_filter:
    input:
        "results/02.snp_calling/{sample}/{sample}_raw.vcf"
    output:
        fsnp = temp("results/02.snp_calling/{sample}/{sample}_fSNP.flt1.vcf"),
        vsnp = temp("results/02.snp_calling/{sample}/{sample}_vSNP.flt1.vcf")
    params:
        # -i, --include <expr>: include only sites for which the expression is true (see man page for details)
        ## (FORMAT/AD[0:1])/(FORMAT/DP) > 0.9: Fixed SNPs: frequency > 90%
        ## FORMAT/ADF[0:1] >=5 && FORMAT/ADR[0:1] >=5: Alternative basecalls were supported by at least
        ##                                             five reads without strand bias
        ##                                             (both strands were mapped by the reads)
        ## QUAL > 20: Variant call quality > 20
        ## FORMAT/DP > 20: Depth > 20
        include_fsnp = "-i '(FORMAT/AD[0:1])/(FORMAT/DP) > 0.9 && FORMAT/ADF[0:1] >=5 "
                 "&& FORMAT/ADR[0:1] >=5 && QUAL > 20 && FORMAT/DP > 20 && FORMAT/AD > 20'",
        include_vsnp = "-i '(FORMAT/AD[0:1])/(FORMAT/DP) > 0.1 && FORMAT/ADF[0:1] >=5 "
                "&& FORMAT/ADR[0:1] >=5 && QUAL > 20 && FORMAT/DP > 20'",
        # -O, --output-type <b|u|z|v>:   b: compressed BCF, u: uncompressed BCF,
        #                                z: compressed VCF, v: uncompressed VCF [v]
        output_type = "-O v",
        # -g, --SnpGap <int>[:type]     filter SNPs within <int> base pairs of an indel (the default)
        #                               or any combination of indel,mnp,bnd,other,overlap
        snpgap = "-g 4"
    threads:
        config["threads"]
    log:
        "logs/02.snp_calling/bcftools_filter/{sample}.log"
    conda:
        "../envs/MTBAnalysis.yaml"
    shell:
        """
        # bcftools filter
        # fSNP
        bcftools filter --threads {threads} {params.include_fsnp} {params.output_type} {params.snpgap} -o {output.fsnp} {input} 2> {log}
        # vSNP
        bcftools filter --threads {threads} {params.include_vsnp} {params.output_type} {params.snpgap} -o {output.vsnp} {input} 2>> {log}
        """

# filter SNP calling results using vcftools
## 1. Discard SNPs if a threshold of SNPs in a 10bp sliding window
## 2. At same time exclude SNPs in our final list of excluding regions, including high-GC regions (PE/PPE/PGRS) and
##    repetitive regions.
rule vcftools_filter:
    input:
        fsnp = "results/02.snp_calling/{sample}/{sample}_fSNP.flt1.vcf",
        vsnp = "results/02.snp_calling/{sample}/{sample}_vSNP.flt1.vcf"
    output:
        fsnp = "results/02.snp_calling/{sample}/{sample}_fSNP.flt.vcf",
        vsnp = "results/02.snp_calling/{sample}/{sample}_vSNP.flt.vcf",
        indel = "results/02.snp_calling/{sample}/{sample}_indel.flt.vcf"
    params:
        # threshold of SNPs in a 10bp sliding window
        thres = 3,
        # final list of excluding regions which includes high-GC regions (PE/PPE/PGRS) and repetitive regions
        exclude = config["excluding_list"]
    log:
        "logs/02.snp_calling/vcftools_filter/{sample}.log"
    conda:
        "../envs/MTBAnalysis.yaml"
    shell:
        """
        # vcftools filter
        # fSNP
        (cat {input.fsnp} | vcf-annotate --filter SnpCluster={params.thres},10 | \
        vcftools --vcf - --stdout --exclude-bed {params.exclude} --remove-filtered-all --recode --recode-INFO-all | \
        bcftools view --types snps > {output.fsnp}) 2> {log}
        # vSNP
        (cat {input.vsnp} | vcf-annotate --filter SnpCluster={params.thres},10 | \
        vcftools --vcf - --stdout --exclude-bed {params.exclude} --remove-filtered-all --recode --recode-INFO-all | \
        bcftools view --types snps > {output.vsnp}) 2>> {log}
        # InDel
        (cat {input.vsnp} | vcf-annotate --filter SnpCluster={params.thres},10 | \
        vcftools --vcf - --stdout --exclude-bed {params.exclude} --remove-filtered-all --recode --recode-INFO-all | \
        bcftools view --types indels > {output.indel}) 2>> {log}
        """

# Annotate filtered VCF using snpEff
rule snpEff_annotate:
    input:
        fsnp = "results/02.snp_calling/{sample}/{sample}_fSNP.flt.vcf",
        vsnp = "results/02.snp_calling/{sample}/{sample}_vSNP.flt.vcf",
        indel = "results/02.snp_calling/{sample}/{sample}_indel.flt.vcf"
    output:
        fsnp_vcf_ann = "results/02.snp_calling/{sample}/{sample}_fSNP.flt.ann.vcf",
        vsnp_vcf_ann = "results/02.snp_calling/{sample}/{sample}_vSNP.flt.ann.vcf",
        indel_vcf_ann = "results/02.snp_calling/{sample}/{sample}_indel.flt.ann.vcf",
        # HTML summary report
        fsnp_summary_html = "results/02.snp_calling/{sample}/{sample}_fSNP_snpEff_summary.html",
        vsnp_summary_html = "results/02.snp_calling/{sample}/{sample}_vSNP_snpEff_summary.html",
        indel_summary_html = "results/02.snp_calling/{sample}/{sample}_indel_snpEff_summary.html",
        # CSV summary report
        fsnp_summary_csv = "results/02.snp_calling/{sample}/{sample}_fSNP_snpEff_summary.csv",
        vsnp_summary_csv = "results/02.snp_calling/{sample}/{sample}_vSNP_snpEff_summary.csv",
        indel_summary_csv = "results/02.snp_calling/{sample}/{sample}_indel_snpEff_summary.csv"
    log:
        "logs/02.snp_calling/snpEff_annotate/{sample}.log"
    params:
        # snpEff config file to use local MTB database
        # -c , -config: Specify config file
        snpEff_config = config["snpEff_config"]
    conda:
        "../envs/MTBAnalysis.yaml"
    shell:
        """
        # Modify the chromosome name to match SnpEff database format and Annotation by snpEff ann
        # fSNP
        cat {input.fsnp} | sed 's/^NC_000962.3/Chromosome/g' | \
        snpEff ann -v -stats {output.fsnp_summary_html} -c {params.snpEff_config} \
        -csvStats {output.fsnp_summary_csv} Mycobacterium_tuberculosis_h37rv > {output.fsnp_vcf_ann} 2> {log}
        # vSNP
        cat {input.vsnp} | sed 's/^NC_000962.3/Chromosome/g' | \
        snpEff ann -v -stats {output.vsnp_summary_html} -c {params.snpEff_config} \
        -csvStats {output.vsnp_summary_csv} Mycobacterium_tuberculosis_h37rv > {output.vsnp_vcf_ann} 2>> {log}
        # InDel
        cat {input.indel} | sed 's/^NC_000962.3/Chromosome/g' | \
        snpEff ann -v -stats {output.indel_summary_html} -c {params.snpEff_config} \
        -csvStats {output.indel_summary_csv} Mycobacterium_tuberculosis_h37rv > {output.indel_vcf_ann} 2>> {log}
        """