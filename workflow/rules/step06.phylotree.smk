##### Step 06.phylotree #####

# Phylogenetic tree construction by RAxML-NG with M.canettii added as the tree root
rule phylotree_raxml_all:
    input:
        # SNPs concatenated fasta file generated by check_snp_diff.py
        fasta = "results/04.snp_distance/all_samples_snp_genotype_canettii_added.fasta"
    output:
        "results/06.phylotree/RAxML/phyloTree.raxml.rba",
        "results/06.phylotree/RAxML/phyloTree.raxml.startTree",
        "results/06.phylotree/RAxML/phyloTree.raxml.bestTree",
        "results/06.phylotree/RAxML/phyloTree.raxml.mlTrees",
        "results/06.phylotree/RAxML/phyloTree.raxml.support",
        "results/06.phylotree/RAxML/phyloTree.raxml.bestModel",
        "results/06.phylotree/RAxML/phyloTree.raxml.bootstraps",
        "results/06.phylotree/RAxML/phyloTree.raxml.log"
    params:
        model = config["phylotree_model"],
        out_dir = "results/06.phylotree/RAxML",
        # --seed: for pseudo-random number generator (default: current time)
        seed = 12345,
        threads = config["threads_phylotree"],
        # --bs-trees: number of bootstraps replicates
        bs_tree = 100
    log:
        "logs/06.phylotree/phylotree_raxml_all.log"
    conda:
        "../envs/phyloTree.yaml"
    shell:
        "raxml-ng --all --msa {input.fasta} -model {params.model} --prefix {params.out_dir}/phyloTree "
        "--seed {params.seed} --threads {params.threads} --bs-trees {params.bs_tree} --outgroup canettii 2>&1 > {log}"

# Mapping the bootstrap support values onto the ML tree
rule phylotree_raxml_support:
    input:
        best_tree = "results/06.phylotree/RAxML/phyloTree.raxml.bestTree",
        bs_tree = "results/06.phylotree/RAxML/phyloTree.raxml.bootstraps"
    output:
        "results/06.phylotree/RAxML/phyloTree.raxml.bestTree.raxml.support",
        "results/06.phylotree/RAxML/phyloTree.raxml.bestTree.raxml.log"
    log:
        "logs/06.phylotree/phylotree_raxml_support.log"
    conda:
        "../envs/phyloTree.yaml"
    shell:
        "raxml-ng --support --tree {input.best_tree} --bs-trees {input.bs_tree} 2>&1 > {log}"

# rule phylotree_ete3 to 1. reroot on M.canettii and 2. prune node canettii and node REF
## input tree file phyloTree.raxml.bestTree.raxml.support generated by RAxML
## output modified tree file
rule phylotree_ete3:
    input:
        "results/06.phylotree/RAxML/phyloTree.raxml.bestTree.raxml.support"
    output:
        "results/06.phylotree/RAxML/phyloTree.final"
    params:
        script = "workflow/scripts/phylotree_ete3.py"
    log:
        "logs/06.phylotree/phylotree_ete3.log"
    conda:
        "../envs/phyloTree.yaml"
    shell:
        "python {params.script} {input} {output} 2>&1 > {log}"


# rule phylotree_ggtree to generate phyloTree visulization result
## input modified tree file by ETE3
## output phyloTree visulization result
rule phylotree_ggtree:
    input:
        tree = "results/06.phylotree/RAxML/phyloTree.final",
        ggtree_dr_indiv = "results/05.lineage_dr/all.dr.indiv.ggtree.txt",
        ggtree_dr = "results/05.lineage_dr/all.dr.ggtree.txt",
        ggtree_lineage = "results/05.lineage_dr/all.lineage.ggtree.txt"
    output:
        "results/06.phylotree/Phylogenetic_Tree.png",
        "results/06.phylotree/Phylogenetic_Tree.pdf"
    params:
        script = "workflow/scripts/ggtree.R"
    log:
        "logs/06.phylotree/phylotree_ggtree.log"
    conda:
        "../envs/R_ggtree.yaml"
    shell:
        """
        Rscript {params.script} {input.tree} {input.ggtree_lineage} \
        {input.ggtree_dr} {input.ggtree_dr_indiv} {output} 2> {log}
        rm -f Rplots.pdf
        """

