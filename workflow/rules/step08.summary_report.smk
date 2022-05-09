##### Step 08.summary_report #####
"""
Note:
For now, phyloTree visualization results are generated by R:ggtree which makes the whole workflow ONLY run once.
Some issues might happen in R::ggtree with long sample names or large sample size.

Use iTOL website version (batch mode needs subscription) ONLY if you want to generate more sophisticated plots or if you
meet issues with R:ggtree. Using iTOL website version makes the workflow run twice to get the final results.
Please refer to following instructions:
1. Since phyloTree has to be generated from a website app iTol (no command line version),
you have to run snakemake twice in order to get summary report.
2. Firstly, run snakemake the 1st time to get all result files.
3. Delete "results/06.phylotree/Phylogenetic_Tree.png" and "results/07.summary_report/summary_report.html".
4. Use iTol to generate Phylogenetic_Tree.png and put it as "results/06.phyloTree/Phylogenetic_Tree.png".
5. Rerun snakemake the 2nd time to get the final summary report.
"""

# Generate summary report using a Rmarkdown script
rule summary_report:
    input:
        # TB-Profiler merged result file from step05.lineage_dr
        tbprofiler_all = "results/05.lineage_dr/all.txt",
        # All samples SNP distance matrix from step04.snp_distance
        snp_diff_matrix = "results/04.snp_distance/all_samples_snp_diff.txt",
        # SNP distance heatmap from step04.snp_distance
        snp_diff_heatmap= "results/04.snp_distance/all_samples_snp_diff_heatmap.png",
        # phyloTree from step06.phyloTree
        phylotree = "results/06.phylotree/Phylogenetic_Tree.png",
        # Kraken result from step00.data_cleaning
        kraken_result = "results/00.data_cleaning/kraken.result",
        # MixInfect result from step07.mixed_infection
        mixinfect_result = "results/07.mix_infection/ALL_MixedSampleSummary.txt"
    output:
        "results/08.summary_report/summary_report.html"
    params:
        script = "workflow/scripts/run_rmd.R",
        rmd_script = "workflow/scripts/generate_summary_html.Rmd",
        percentage_MTBC_threshold = config["kraken_threshold"]
    log:
        "logs/08.summary_report/summary_report.log"
    conda:
        "../envs/summary_report.yaml"
    shell:
        # Run a R script to pass parameters to Rmd script
        """
        Rscript --vanilla {params.script} {params.rmd_script} {output} {input.tbprofiler_all} {input.snp_diff_matrix} \
        {input.snp_diff_heatmap} {input.phylotree} {input.kraken_result} {params.percentage_MTBC_threshold} \
        {input.mixinfect_result} > {log} 2>&1
        """


        