# pass parameters
args <- commandArgs(trailingOnly=TRUE)
rmd_script <- args[1]
output <- args[2]
tbprofiler <- args[3]
snp_diff_matrix <- args[4]
snp_diff_heatmap <- args[5]
phylotree <- args[6]
kraken <- args[7]
percentage_MTBC_threshold <- args[8]
mixinfect <- args[9]

# absolute path for these files otherwise Rmarkdown will raise an error
tbprofiler <- paste0(normalizePath(dirname(tbprofiler)), '/',  basename(tbprofiler))
snp_diff_matrix <- paste0(normalizePath(dirname(snp_diff_matrix)), '/',  basename(snp_diff_matrix))
snp_diff_heatmap <- paste0(normalizePath(dirname(snp_diff_heatmap)), '/',  basename(snp_diff_heatmap))
phylotree <- paste0(normalizePath(dirname(phylotree)), '/',  basename(phylotree))
kraken <- paste0(normalizePath(dirname(kraken)), '/',  basename(kraken))
mixinfect <- paste0(normalizePath(dirname(mixinfect)), '/',  basename(mixinfect))

# run Rmarkdown script
rmarkdown::render(
  rmd_script,
  output_file = basename(output),
  output_dir = dirname(output),
  params=list(
    tbprofiler = tbprofiler,
    snp_diff_matrix = snp_diff_matrix,
    snp_diff_heatmap = snp_diff_heatmap,
    phylotree = phylotree,
    kraken = kraken,
    percentage_MTBC_threshold = percentage_MTBC_threshold,
    mixinfect = mixinfect
  )
)
