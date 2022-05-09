##### Step 03.vcf_merging #####

# Compress VCF files using bgzip and index compressed VCF files by bcftools index
rule vcf_bgzip_index:
    input:
        fsnp = "results/02.snp_calling/{sample}/{sample}_fSNP.flt.vcf",
        vsnp = "results/02.snp_calling/{sample}/{sample}_vSNP.flt.vcf",
        indel = "results/02.snp_calling/{sample}/{sample}_indel.flt.vcf"
    output:
        fsnp_vcf_gz = "results/02.snp_calling/{sample}/{sample}_fSNP.flt.vcf.gz",
        fsnp_vcf_gz_index = "results/02.snp_calling/{sample}/{sample}_fSNP.flt.vcf.gz.tbi",
        vsnp_vcf_gz = "results/02.snp_calling/{sample}/{sample}_vSNP.flt.vcf.gz",
        vsnp_vcf_gz_index ="results/02.snp_calling/{sample}/{sample}_vSNP.flt.vcf.gz.tbi",
        indel_vcf_gz = "results/02.snp_calling/{sample}/{sample}_indel.flt.vcf.gz",
        indel_vcf_gz_index = "results/02.snp_calling/{sample}/{sample}_indel.flt.vcf.gz.tbi"
    threads:
        config["threads"]
    log:
        "logs/03.vcf_merging/vcf_bgzip_index/{sample}.log"
    conda:
        "../envs/MTBAnalysis.yaml"
    shell:
        """
        # Compress VCF files using bgzip 
        bgzip -c -f -@ {threads} {input.fsnp} > {output.fsnp_vcf_gz} 2> {log}
        bgzip -c -f -@ {threads} {input.vsnp} > {output.vsnp_vcf_gz} 2>> {log}
        bgzip -c -f -@ {threads} {input.indel} > {output.indel_vcf_gz} 2>> {log}
        # index compressed VCF files by bcftools index
        ## bcftools index: Index bgzip compressed VCF/BCF files for random access.
        bcftools index --tbi --threads {threads} -o {output.fsnp_vcf_gz_index} {output.fsnp_vcf_gz} 2>> {log}
        bcftools index --tbi --threads {threads} -o {output.vsnp_vcf_gz_index} {output.vsnp_vcf_gz} 2>> {log}
        bcftools index --tbi --threads {threads} -o {output.indel_vcf_gz_index} {output.indel_vcf_gz} 2>> {log}
        """

# Write a VCF files list
rule vcf_list_write:
    input:
        fsnp = lambda wildcards: get_qualified_results("results/02.snp_calling/{sample}/{sample}_fSNP.flt.vcf.gz",wildcards),
        vsnp = lambda wildcards: get_qualified_results("results/02.snp_calling/{sample}/{sample}_vSNP.flt.vcf.gz",wildcards),
        indel = lambda wildcards: get_qualified_results("results/02.snp_calling/{sample}/{sample}_indel.flt.vcf.gz",wildcards)
    output:
        fsnp = "results/03.vcf_merging/vcf_list_fsnp.txt",
        vsnp = "results/03.vcf_merging/vcf_list_vsnp.txt",
        indel = "results/03.vcf_merging/vcf_list_indel.txt"
    log:
        "logs/03.vcf_merging/vcf_list_write.log"
    shell:
        """
        echo {input.fsnp} | sed 's# #\\n#g' > {output.fsnp} 2> {log}
        echo {input.vsnp} | sed 's# #\\n#g' > {output.vsnp} 2>> {log}
        echo {input.indel} | sed 's# #\\n#g' > {output.indel} 2>> {log}
        """

# Merge VCF files using bcftools merge
rule bcftools_merge_fsnp:
    input:
        "results/03.vcf_merging/vcf_list_fsnp.txt"
    output:
        merged_vcf = "results/03.vcf_merging/merged_fsnp.vcf",
        merged_vcf_gz = "results/03.vcf_merging/merged_fsnp.vcf.gz"
    params:
        # -0  --missing-to-ref: assume genotypes at missing sites are 0/0
        missing_to_ref = "--missing-to-ref",
        # -l, --file-list <file>: read file names from the file
        # -O, --output-type <b|u|z|v>: 'b' compressed BCF; 'u' uncompressed BCF;
        #                              'z' compressed VCF; 'v' uncompressed VCF [v]
        output_dir = lambda wildcards, output: extract_dirname(wildcards, output[0])
    threads:
        config["threads"]
    log:
        "logs/03.vcf_merging/bcftools_merge_fsnp.log"
    conda:
        "../envs/MTBAnalysis.yaml"
    shell:
        """
        # split large sample size into several 500-sample batch to prevent a bcftool merge command issue
        mkdir -p {params.output_dir}/tmp_fsnp 2> {log}
        cat {input} | split -l 500 - {params.output_dir}/tmp_fsnp/subset_vcfs 2>> {log}
        
        for i in {params.output_dir}/tmp_fsnp/subset_vcfs*
        do
		  sort -n -o $i $i  2>> {log}
          # bcftools merge: Merge multiple VCF/BCF files from non-overlapping sample sets to 
          #                 create one multi-sample file.
          #                 Note that only records from different files can be merged, never from the same file. For
          #                 "vertical" merge take a look at "bcftools norm" instead.
          bcftools merge {params.missing_to_ref} -l $i -O v -o $i.merged.vcf 2>> {log}
          bgzip -c -f -@ 16 $i.merged.vcf > $i.merged.vcf.gz 2>> {log}
          bcftools index --tbi --threads 16 -o $i.merged.vcf.gz.tbi $i.merged.vcf.gz 2>> {log}
        done
        
        ls {params.output_dir}/tmp_fsnp/*.merged.vcf.gz > {params.output_dir}/tmp_fsnp/vcf_list.txt 2>> {log}
		sort -n -o {params.output_dir}/tmp_fsnp/vcf_list.txt {params.output_dir}/tmp_fsnp/vcf_list.txt 2>> {log}
        
        bcftools merge {params.missing_to_ref} -l {params.output_dir}/tmp_fsnp/vcf_list.txt -O v -o {output.merged_vcf} 2>> {log}
        bgzip -c -f -@ {threads} {output.merged_vcf} > {output.merged_vcf_gz} 2>> {log}
        
        rm -rf {params.output_dir}/tmp_fsnp 2>> {log}
        """

rule bcftools_merge_vsnp:
    input:
        "results/03.vcf_merging/vcf_list_vsnp.txt"
    output:
        merged_vcf = "results/03.vcf_merging/merged_vsnp.vcf",
        merged_vcf_gz = "results/03.vcf_merging/merged_vsnp.vcf.gz"
    params:
        # -0  --missing-to-ref: assume genotypes at missing sites are 0/0
        missing_to_ref = "--missing-to-ref",
        # -l, --file-list <file>: read file names from the file
        # -O, --output-type <b|u|z|v>: 'b' compressed BCF; 'u' uncompressed BCF;
        #                              'z' compressed VCF; 'v' uncompressed VCF [v]
        output_dir = lambda wildcards, output: extract_dirname(wildcards,output[0])
    threads:
        config["threads"]
    log:
        "logs/03.vcf_merging/bcftools_merge_vsnp.log"
    conda:
        "../envs/MTBAnalysis.yaml"
    shell:
        """
        # split large sample size into several 500-sample batch to prevent a bcftool merge command issue
        mkdir -p {params.output_dir}/tmp_vsnp 2> {log}
        cat {input} | split -l 500 - {params.output_dir}/tmp_vsnp/subset_vcfs 2>> {log}

        for i in {params.output_dir}/tmp_vsnp/subset_vcfs*
        do
		  sort -n -o $i $i  2>> {log}
          # bcftools merge: Merge multiple VCF/BCF files from non-overlapping sample sets to 
          #                 create one multi-sample file.
          #                 Note that only records from different files can be merged, never from the same file. For
          #                 "vertical" merge take a look at "bcftools norm" instead.
          bcftools merge {params.missing_to_ref} -l $i -O v -o $i.merged.vcf 2>> {log}
          bgzip -c -f -@ 16 $i.merged.vcf > $i.merged.vcf.gz 2>> {log}
          bcftools index --tbi --threads 16 -o $i.merged.vcf.gz.tbi $i.merged.vcf.gz 2>> {log}
        done

        ls {params.output_dir}/tmp_vsnp/*.merged.vcf.gz > {params.output_dir}/tmp_vsnp/vcf_list.txt 2>> {log}
		sort -n -o {params.output_dir}/tmp_vsnp/vcf_list.txt {params.output_dir}/tmp_vsnp/vcf_list.txt 2>> {log}

        bcftools merge {params.missing_to_ref} -l {params.output_dir}/tmp_vsnp/vcf_list.txt -O v -o {output.merged_vcf} 2>> {log}
        bgzip -c -f -@ {threads} {output.merged_vcf} > {output.merged_vcf_gz} 2>> {log}

        rm -rf {params.output_dir}/tmp_vsnp 2>> {log}
        """

rule bcftools_merge_indel:
    input:
        "results/03.vcf_merging/vcf_list_indel.txt"
    output:
        merged_vcf = "results/03.vcf_merging/merged_indel.vcf",
        merged_vcf_gz = "results/03.vcf_merging/merged_indel.vcf.gz"
    params:
        # -0  --missing-to-ref: assume genotypes at missing sites are 0/0
        missing_to_ref = "--missing-to-ref",
        # -l, --file-list <file>: read file names from the file
        # -O, --output-type <b|u|z|v>: 'b' compressed BCF; 'u' uncompressed BCF;
        #                              'z' compressed VCF; 'v' uncompressed VCF [v]
        output_dir = lambda wildcards, output: extract_dirname(wildcards,output[0])
    threads:
        config["threads"]
    log:
        "logs/03.vcf_merging/bcftools_merge_indel.log"
    conda:
        "../envs/MTBAnalysis.yaml"
    shell:
        """
        # split large sample size into several 500-sample batch to prevent a bcftool merge command issue
        mkdir -p {params.output_dir}/tmp_indel 2> {log}
        cat {input} | split -l 500 - {params.output_dir}/tmp_indel/subset_vcfs 2>> {log}

        for i in {params.output_dir}/tmp_indel/subset_vcfs*
        do
		  sort -n -o $i $i  2>> {log}
          # bcftools merge: Merge multiple VCF/BCF files from non-overlapping sample sets to 
          #                 create one multi-sample file.
          #                 Note that only records from different files can be merged, never from the same file. For
          #                 "vertical" merge take a look at "bcftools norm" instead.
          bcftools merge {params.missing_to_ref} -l $i -O v -o $i.merged.vcf 2>> {log}
          bgzip -c -f -@ 16 $i.merged.vcf > $i.merged.vcf.gz 2>> {log}
          bcftools index --tbi --threads 16 -o $i.merged.vcf.gz.tbi $i.merged.vcf.gz 2>> {log}
        done

        ls {params.output_dir}/tmp_indel/*.merged.vcf.gz > {params.output_dir}/tmp_indel/vcf_list.txt 2>> {log}
		sort -n -o {params.output_dir}/tmp_indel/vcf_list.txt {params.output_dir}/tmp_indel/vcf_list.txt 2>> {log}

        bcftools merge {params.missing_to_ref} -l {params.output_dir}/tmp_indel/vcf_list.txt -O v -o {output.merged_vcf} 2>> {log}
        bgzip -c -f -@ {threads} {output.merged_vcf} > {output.merged_vcf_gz} 2>> {log}

        rm -rf {params.output_dir}/tmp_indel 2>> {log}
        """

# Annotate merged VCF file using snpEff
rule vcf_merge_annotate:
    input:
        fsnp = "results/03.vcf_merging/merged_fsnp.vcf",
        vsnp = "results/03.vcf_merging/merged_vsnp.vcf",
        indel = "results/03.vcf_merging/merged_indel.vcf"
    output:
        fsnp_vcf_merge_ann = "results/03.vcf_merging/merged_fsnp.ann.vcf",
        vsnp_vcf_merge_ann = "results/03.vcf_merging/merged_vsnp.ann.vcf",
        indel_vcf_merge_ann = "results/03.vcf_merging/merged_indel.ann.vcf",
        # HTML summary report
        fsnp_summary_html = "results/03.vcf_merging/merged_fsnp_snpEff_summary.html",
        vsnp_summary_html = "results/03.vcf_merging/merged_vsnp_snpEff_summary.html",
        indel_summary_html = "results/03.vcf_merging/merged_indel_snpEff_summary.html",
        # CSV summary report
        fsnp_summary_csv = "results/03.vcf_merging/merged_fsnp_snpEff_summary.csv",
        vsnp_summary_csv = "results/03.vcf_merging/merged_vsnp_snpEff_summary.csv",
        indel_summary_csv = "results/03.vcf_merging/merged_indel_snpEff_summary.csv"
    log:
        "logs/03.vcf_merging/vcf_merge_annotate.log"
    params:
        # snpEff config file to use local MTB database
        # -c , -config: Specify config file
        snpEff_config = config["snpEff_config"]
    conda:
        "../envs/MTBAnalysis.yaml"
    shell:
        """
        # Modify the chromosome name to match SnpEff database format and Annotation by snpEff ann
        # fsnp
        (cat {input.fsnp} | sed 's/^NC_000962.3/Chromosome/g' | \
        snpEff ann -v -stats {output.fsnp_summary_html} -c {params.snpEff_config} \
        -csvStats {output.fsnp_summary_csv} Mycobacterium_tuberculosis_h37rv > {output.fsnp_vcf_merge_ann}) 2> {log}
        # vsnp
        (cat {input.vsnp} | sed 's/^NC_000962.3/Chromosome/g' | \
        snpEff ann -v -stats {output.vsnp_summary_html} -c {params.snpEff_config} \
        -csvStats {output.vsnp_summary_csv} Mycobacterium_tuberculosis_h37rv > {output.vsnp_vcf_merge_ann}) 2>> {log}
        # vsnp
        (cat {input.indel} | sed 's/^NC_000962.3/Chromosome/g' | \
        snpEff ann -v -stats {output.indel_summary_html} -c {params.snpEff_config} \
        -csvStats {output.indel_summary_csv} Mycobacterium_tuberculosis_h37rv > {output.indel_vcf_merge_ann}) 2>> {log}
        """

