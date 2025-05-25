############
# Libaries #
############

library(tidyverse)
library(VennDetail)
library(ggvenn)
library(ggpubr)

###################
# Read gene lists #
###################

#cluster row order form TCGA-HNSC heatmap 
tcga <- read.csv2("2025_03_25_limma_trim45_top_bottom_tcga_brca_er_positive_bh.csv")
col <- "logFC"
sorted_tcga <- tcga[order(tcga[[col]]), ]
subset_tcga <- head(sorted_tcga, 100)

#Cluster row order for noroc 
meta <- read.csv2("limma_metabric_er_positive_top_bottom_t45_level_BH.csv")
sorted_meta <- meta[order(meta[[col]]), ]
subset_meta <- head(sorted_meta, 100)

subset_tcga <- subset_tcga[,"hgnc_symbol"]
subset_meta <- subset_meta[,"X"]



##################
## Venn Diagram ##
##################

overlap_list <- list("TCGA-BRCA" = subset_tcga, "METABRIC" = subset_meta)
overlap_plot <- venndetail(overlap_list) #identify overlapping terms
overlap_df <- result(overlap_plot)

#fwrite(overlap_df, "2025_04_11_overlap_tcga_meta_downregulated_genes.csv", row.names = TRUE)


x = list(subset_tcga,subset_meta)

names(x) <- c("TCGA-BRCA", "METABRIC")


plot <- ggvenn(x,set_name_size = 3,
               fill_color = c("#E41A1C", "#377EB8"),stroke_size = 0.5)

plot

dev.off()

ggexport(plot, filename = "2025_04_11_tcga_meta_downregulated_genes_venndiagram.png",res=200,width = 2500, height = 2000)
ggexport(plot, filename = "2025_04_11_tcga_meta_downregulated_genes_venndiagram.pdf",width = 15, height = 15)
