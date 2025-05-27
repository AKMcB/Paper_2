library(tidyverse)
library(ggpubr)
#Upload TCGA file 
ER <- read.csv2("Raw data/gdc_brca_expr_edit.csv",sep = ";", as.is = T, check.names = F)
ER <- ER[,-1]

#Subset the gene of interest
genes <- subset(ER, ER$gene %in% c("TRIM45", "ESR1"))

rownames(genes) <- genes[,1]
genes <- genes[,-1]
genes <- as.data.frame(t(genes))

#Upload the new info file 
info <- read.csv2("Raw data/TCGA_BRCA_Updated_Clinical_Data.csv",sep = ";", as.is = T, check.names = F)
info_rec <- info[,-c(2:44)]

#Merge the genes and infor file based on patient ID 
genes <- tibble::rownames_to_column(genes, "ID")
merged <- merge(genes, info_rec, by.x= "ID", by.y = "TCGA_id")
#write.csv2(merged, "TRIM45_receptor status_TCGA_clinical.csv")

#Subset based on positive ER status 
ER <- subset(merged, merged$er_status_by_ihc %in% "Positive")

selected_genes <- colnames(merged)[2:3]  
for (gene in selected_genes) {
  cor_value <- cor(ER$TRIM45, ER[[gene]], method = "pearson")
  line_color <- ifelse(cor_value > 0, "#D73027", "#0000EE")
  
  plot <- ggscatter(ER, x = "TRIM45", y = gene,
                    title = paste("TRIM45 vs", gene, "ER+ (n =", nrow(ER), ")"),
                    color = "black", shape = 21, size = 3, fill = line_color, alpha = 0.4,
                    add = "reg.line",
                    add.params = list(color = line_color, fill = line_color),
                    conf.int = TRUE,
                    cor.coef = TRUE,
                    cor.coef.size = 5,
                    cor.coeff.args = list(method = "pearson", label.sep = "\n"),
                    xlab = "TRIM45", ylab = gene) +
    theme(legend.position = "bottom",
          plot.title = element_text(hjust = 0.5, face = "italic"),
          axis.text.x = element_text(size = 18), axis.ticks.x = element_blank(),
          axis.title.x = element_text(size = 18), axis.title.y = element_text(size = 18),
          axis.text.y = element_text(size = 18))
  
  pdf_filename <- paste0("TRIM45_", gene, "_Corr_TCGA.pdf")
  
  pdf(pdf_filename, height = 6, width = 6)
  print(plot)
  dev.off()
}

