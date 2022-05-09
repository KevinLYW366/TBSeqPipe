##### Step 04.snp_distance #####

# Extract chromosome, position, ref allele, alternate allele and translated genotype from the merged VCF file
rule bcftools_query:
    input:
        # use fSNP in phylogenetic analysis
        "results/03.vcf_merging/merged_fsnp.vcf"
    output:
        "results/04.snp_distance/all_samples_snp_genotype.txt"
    params:
        dr_genes = config["dr_genes_list"],
        outdir = "results/04.snp_distance"
    log:
        "logs/04.snp_distance/bcftools_query.log"
    threads:
        config["threads"]
    conda:
        "../envs/MTBAnalysis.yaml"
    shell:
        """
        # remove SNPs in drug-resistance associated genes before SNP distance calculation
        bedtools intersect -v -a {input} -b {params.dr_genes} -wa -header | bgzip -c -f -@ {threads} > {params.outdir}/tmp.vcf.gz 2> {log}
        
        # bcftools query: Extracts fields from VCF/BCF file and prints them in user-defined format
        bcftools query -H -f '%CHROM\\t%POS\\t%REF\\t%ALT[\\t%TGT]\\n' {params.outdir}/tmp.vcf.gz > {output} 2>> {log}
        
        rm -f {params.outdir}/tmp.vcf.gz 2>> {log}
        """

# Output a fasta file with all SNP concatenated and
# calculate the number of different SNP between each sample pair
# using in-house python script
rule snp_diff_calculate:
    input:
        "results/04.snp_distance/all_samples_snp_genotype.txt"
    output:
        # a fasta file with all SNP concatenated
        snp_fasta = "results/04.snp_distance/all_samples_snp_genotype.fasta",
        snp_fasta_canettii_added = "results/04.snp_distance/all_samples_snp_genotype_canettii_added.fasta",
        # the number of different SNP between each sample pair
        snp_diff = "results/04.snp_distance/all_samples_snp_diff.txt",
        snp_diff_canettii_added = "results/04.snp_distance/all_samples_snp_diff_canettii_added.txt"
    params:
        out_dir = "results/04.snp_distance",
        data_dir = config["data_dir"],
        canettii_flt_snps = config["canettii_flt_snps"],
        script = "workflow/scripts/check_snp_diff.py",
        script_canettii_added = "workflow/scripts/check_snp_diff_canettii_added.py"
    log:
        "logs/04.snp_distance/snp_diff_calculate.log"
    conda:
        "../envs/MTBAnalysis.yaml"
    shell:
        """
        python {params.script} {params.out_dir} {params.data_dir} 2 2> {log}
        python {params.script_canettii_added} {params.out_dir} {params.canettii_flt_snps} {params.data_dir} 2 2>> {log}
        """

rule snp_diff_heatmap:
    input:
        snp_diff = "results/04.snp_distance/all_samples_snp_diff.txt",
        tbprofiler_all = "results/05.lineage_dr/all.txt"
    output:
        pdf = "results/04.snp_distance/all_samples_snp_diff_heatmap.pdf",
        png = "results/04.snp_distance/all_samples_snp_diff_heatmap.png"
    params:
        script = "workflow/scripts/heatmap.R",
        show_rownames = "TRUE",
        show_colnames = "TRUE",
        cellwidth = 15,
        cellheight = 15
    log:
        "logs/04.snp_distance/snp_diff_heatmap.log"
    conda:
        "../envs/R_pheatmap.yaml"
    shell:
        "Rscript --vanilla {params.script} {input.snp_diff} {input.tbprofiler_all} {output.pdf} "
        "{params.show_colnames} {params.show_rownames} {params.cellwidth} {params.cellheight} {output.png} 2> {log}"