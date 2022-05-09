##### Step 07.mix_infection #####

# Detect mix infection from diploid type version vSNP VCF file using MixInfect (https://github.com/bensobkowiak/MixInfect)

# filter out false postive SNPs as much as possible to increase the accuracy of mix infection detection
# refined regions of low confidence list is from 
## Marin, Maximillian et al. “Benchmarking the empirical accuracy of short-read sequencing across the M. tuberculosis genome.” 
## Bioinformatics (Oxford, England), btac023. 10 Jan. 2022, doi:10.1093/bioinformatics/btac023
rule mixinfect_filter:
    input:
        "results/02.snp_calling/{sample}/{sample}_vSNP.flt.vcf"
    output:
        "results/07.mix_infection/{sample}/{sample}_vSNP.flt.mixinfect.vcf"
    params:
        regions = config["refined_low_confidence_regions"]
    log:
        "logs/07.mix_infection/mixinfect_filter/{sample}.log"
    conda:
        "../envs/MTBAnalysis.yaml"
    shell:
        """
        (cat {input} | vcf-annotate --filter MinMQ=20 | \
        vcftools --vcf - --stdout --exclude-bed {params.regions} \
        --remove-filtered-all --recode --recode-INFO-all > {output}) 2> {log}
        """

# Detect mix infection
rule mixinfect_detect:
    input:
        #"results/07.mix_infection/{sample}/{sample}_SNP_diploid.flt.vcf"
        "results/07.mix_infection/{sample}/{sample}_vSNP.flt.mixinfect.vcf"
    output:
        "results/07.mix_infection/{sample}/{sample}_MixedSampleSummary.txt"
    params:
        script = "workflow/scripts/MixInfectVCF.R"
    log:
        "logs/07.mix_infection/mixinfect_detect/{sample}.log"
    conda:
        "../envs/MixInfect.yaml"
    shell:
        "Rscript {params.script} {input} {output} > {log} 2>&1"

# Merge mix infection results
rule mixinfect_merge:
    input:
        lambda wildcards: get_qualified_results("results/07.mix_infection/{sample}/{sample}_MixedSampleSummary.txt", wildcards)
    output:
        "results/07.mix_infection/ALL_MixedSampleSummary.txt"
    log:
        "logs/07.mix_infection/mixinfect_merge.log"
    shell:
        """
        echo -e "##Potential mixed infection detected by MixInect for each sample" > {output}
        echo -e "#sample\tmix\tmix_snps\ttot_snps\thet_snps_proportion\tn_mixed_strains\tmajor_strain_proportion" >> {output}
        for result in {input}
        do
          cp {output} {output}.tmp
          tail -n1 ${{result}} | cat {output}.tmp - > {output}
          rm {output}.tmp
        done
        TMP=`head -n 2 {output} && grep -v "^#" {output} | sort -k1` 
        echo "$TMP" > {output}
        """

