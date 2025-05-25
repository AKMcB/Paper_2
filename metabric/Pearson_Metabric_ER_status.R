library(tidyverse)
library(ggpubr)

#Upload TCGA file 
ER <- read.csv2("Ensembl_Metabric_17724_duplicates removed .csv",sep = ";", as.is = T, check.names = F)
ER <- ER[,-c(2,3)]

#Subset the gene of interest
genes <- subset(ER, ER$Hugo_Symbol %in% c("TRIM45", "KLF4", 
                                          "MYBL1", "HES1", "ASCL1", 
                                          "ELF3", "FOS", "MYB", "BATF", 
                                          "SOX3", "MYC", "ELF1", 
                                          "BHLHE40", "RARA", 
                                          "OVOL2", "TGIF2", 
                                          "BCL11B", "FOXC1"))
rownames(genes) <- genes[,1]
genes <- genes[,-1]
genes <- t(genes)
genes <- as.data.frame(genes)


#Upload the new info file 
info <- read.csv2("brca_metabric_clinical_data.csv",sep = ";", as.is = T, check.names = F)
info_rec <- info[,-c(1,3:11)]
info_rec <- info_rec[,-c(4,5,7:29)]

#Merge the genes and infor file based on patient ID 
genes <- tibble::rownames_to_column(genes, "ID")
merged <- merge(genes, info_rec, by.x= "ID", by.y = "Patient ID")
#write.csv2(merged, "TRIM45_receptor_status_Metabric.csv")

# 1)Subset based on positive ER status 
er <- subset(merged, merged$`ER status measured by IHC` %in% "Positve")

#ER <- read.csv2("T45_ESR1_ER_pos_Meta.csv", as.is = T, check.names = F)
selected_genes <- colnames(er)[2:18]  # Replace with your specific selection criteria if needed
# Loop through each gene in the selected_genes vector
# Loop through each gene in the selected_genes vector
for (gene in selected_genes) {
  # Calculate the correlation coefficient
  cor_value <- cor(er$TRIM45, er[[gene]], method = "pearson")
  
  # Determine the color based on the correlation value
  line_color <- ifelse(cor_value > 0, "#D73027", "#0000EE")
  
  # Create the scatter plot for each gene
  plot <- ggscatter(er, x = "TRIM45", y = gene,
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
  
  # Define the file name for the PDF
  pdf_filename <- paste0("TRIM45_", gene, "_Corr_meta.pdf")
  
  # Save the plot as a PDF
  pdf(pdf_filename, height = 6, width = 6)
  print(plot)
  dev.off()
}

