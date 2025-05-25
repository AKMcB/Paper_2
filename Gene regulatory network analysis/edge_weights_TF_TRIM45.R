#Import the libraries
library(tidyverse)
library(data.table)
library(gson)


#Read in the edge file 
expr <- as.data.frame(fread("C:/Users/abe186/UiT Office 365/O365-Phd Anne - General/NCMM/TCGA-BRCA/lioness_output/lioness_filtered_for_genes_trim45.csv"))
head(expr)[1:5]

#Read in the gene list
list1 <- read.gmt("C:/Users/abe186/UiT Office 365/O365-Phd Anne - General/NCMM/TCGA-BRCA/gene_lists/HALLMARK_ESTROGEN_RESPONSE_EARLY.v2023.2.Hs.gmt")

list2 <- read.gmt("C:/Users/abe186/UiT Office 365/O365-Phd Anne - General/NCMM/TCGA-BRCA/gene_lists/HALLMARK_ESTROGEN_RESPONSE_LATE.v2023.2.Hs.gmt")

####################
#Combine the gene files
genes <- rbind(list1, list2)
nrow(genes)

#Check for duplicates
genes <- genes[!duplicated(genes$gene),]
genes$term<- NULL
nrow(genes) #Removed 18 genes that were duplicated
head(genes)

expr <- subset(expr, expr$TF %in% genes$gene)

expr <- expr[,-2] #all targets are affect TRIM45
rownames(expr) <- expr[,1]
expr <- expr[,-1]


# Calculate the row sums (or use rowMeans for averages)
row_means <- as.data.frame(rowMeans(expr))
row_means

row_means <- rownames_to_column(row_means, "gene")
colnames(row_means)[2] <- "value"

row_means$type <- ifelse(row_means$value >=0, "Positive", "Negative")
# Get the indices of the top 10 and bottom 10 rows
#top_10_indices <- order(row_means, decreasing = TRUE)[1:10]
#bottom_10_indices <- order(row_means, decreasing = FALSE)[1:10]

# Extract the row means for the top and bottom 10 rows
#top_10_means <- row_means[top_10_indices]
#bottom_10_means <- row_means[bottom_10_indices]

# Visualize the rows
#plot_data <- data.frame(
  #Row = c(paste0("Top_", 1:10)),
  #MeanValue = top_10_means,
  #Type = c(rep("Positive", 10)))

#plot_data <- rownames_to_column(plot_data, "gene")

p <-ggplot(row_means, aes(x = reorder(gene, value), y = value, fill = type), color= NA) +
  geom_bar(stat = "identity", )+
  scale_fill_manual(values = c("Positive" = "firebrick2", "Negative" = "steelblue2")) +  
  labs(title = "Top 10 and Bottom 10 TF affecting TRIM45",
       x = "Transcription Factor",
       y = "Mean value") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 18, color = "black"),
        axis.text.y = element_text( size = 18, color = "black"),
    panel.background = element_rect(fill = "white", colour = "white"),
    plot.background = element_rect(fill = "white", colour = "white"),
    axis.line = element_line(linewidth = 0.7, linetype = "solid",
                             colour = "black"))

p



# Save as PNG
png("edge_weight_TF_trim45_er_response_comb.png", width = 8, height = 6, units = "in", res = 300)  # Set resolution to 300 dpi for PNG
plot(p)
dev.off()

# Save as PDF
pdf("edge_weight_TF_TRIM45_er_response_comb.pdf", width = 6, height = 4)  # PDF works with dimensions in inches
plot(p)
dev.off()