# log <- file(snakemake@log[[1]], open="wt")
# sink(log)

library(pheatmap)

# input_snp_diff <- snakemake@input[["snp_diff"]]
# input_tbprofiler_all <- snakemake@input[["tbprofiler_all"]]
# output_heatmap <- snakemake@output[1]

args <- commandArgs(trailingOnly=TRUE)
input_snp_diff <- args[1]
input_tbprofiler_all <- args[2]
output_heatmap_pdf <- args[3]
if (args[4] == "TRUE") {
  show_colnames_arg <- TRUE
} else if (args[4] == "FALSE") {
  show_colnames_arg <- FALSE
}
if (args[5] == "TRUE") {
  show_rownames_arg <- TRUE
} else if (args[5] == "FALSE") {
  show_rownames_arg <- FALSE
}
cellwidth_arg <- as.numeric(args[6])
cellheight_arg <- as.numeric(args[7])
output_heatmap_png <- args[8]

snp_diff <- read.table(input_snp_diff, header = T, check.names = F)
snp_diff$REF <- NULL
snp_diff <- snp_diff[-c(1),]
#head(snp_diff)

tbprofiler <- read.table(input_tbprofiler_all, sep="\t", na.strings="NA", header = T, row.names = 1)
#head(tbprofiler)
tbprofiler <- tbprofiler[, c(1,3)]
#head(tbprofiler)
dr_cat_order <- c("Sensitive", "RR-TB", "HR-TB", "MDR-TB","Pre-XDR-TB", "XDR-TB", "Other")
dr_cat_values <- unique(tbprofiler$DR_type)
dr_cat_levels <- dr_cat_order[dr_cat_order %in% dr_cat_values]
tbprofiler$DR_type <- factor(tbprofiler$DR_type,
                           levels = dr_cat_levels)

pheatmap(snp_diff, main="SNP Distance Heatmap", annotation_col = tbprofiler,
         show_rownames = show_rownames_arg, show_colnames = show_colnames_arg,
         cellwidth = cellwidth_arg, cellheight = cellheight_arg, filename=output_heatmap_pdf)

pheatmap(snp_diff, main="SNP Distance Heatmap", annotation_col = tbprofiler,
         show_rownames = show_rownames_arg, show_colnames = show_colnames_arg,
         cellwidth = cellwidth_arg, cellheight = cellheight_arg, filename=output_heatmap_png)