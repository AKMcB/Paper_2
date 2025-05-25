#############
# Libraries #
#############
library(GEOquery)
library(limma)
library(umap)
library(data.table)
library(ggplot2)
library(ggpubr)

# load series and platform data from GEO
gset <- getGEO("GSE27473", GSEMatrix =TRUE, AnnotGPL=TRUE)[[1]]

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group membership for all samples
gsms <- "undefined"
sml <- strsplit(gsms, split="")[[1]]

# log2 transformation
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(gset) <- log2(ex) }

min(ex)
max(ex)
#fwrite(ex, "expr_log2.csv", row.names = TRUE)
#change gene id 
gene_id <- data.table::fread("GPL570-55999 (4).txt")
colnames(gene_id) <- as.character(gene_id[1,])
gene_id <- gene_id[-1,]

gene_id <- dplyr::select(gene_id, c(ID, `Gene Symbol`))
gene_id$`Gene Symbol` <- gsub("\\///.*", "", gene_id$`Gene Symbol`)

ex <- tibble::rownames_to_column(as.data.frame(ex), "ID")

merged <- merge(ex, gene_id, by = "ID")
merged[merged== ""] <- NA
merged <- na.omit(merged)

dup <- merged[duplicated(merged$`Gene Symbol`),]

library(dplyr)
merged <-  merged %>%
  group_by(`Gene Symbol`) %>%
  dplyr::summarize(across(starts_with("GSM"), mean, .names = "{.col}"), .groups = "drop")

dup <- merged[duplicated(merged$`Gene Symbol`),]

merged <- tibble::column_to_rownames(merged, "Gene Symbol")

#chnage id
pheno <- gset@phenoData@data
attach(pheno)
pheno <- pheno[,c("title", "geo_accession")]

colnames(merged) <- pheno$title[match(colnames(merged),pheno$geo_accession)] 


#boxplot 
merged <- as.data.frame(t(merged))
merged <- tibble::rownames_to_column(merged, "id")

merged <- merge(merged, pheno, by.x = "id", by.y ="title")

merged_1 <- merged[,c("STC2", "id")]
merged_1
merged_1$sample <- c("KD", "KD", "KD",
                     "WT", "WT","WT")
#control smaple is just MCF7 treated with a scrmabled siRNA and a emtpy vehicle. 

unique(merged_1$id)
merged_1$sample <- factor(merged_1$sample,
                          levels = c("WT", "KD"))

my_comparisons <- list( c("WT", "KD")) 
symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
                    symbols = c("****", "***", "**", "*", "ns"))


p <- ggplot(merged_1, aes(x = sample, y = STC2, fill = sample))+
  geom_point(alpha=1,position = position_jitter(), shape= 21, size= 6)+
  geom_boxplot(fill = "white", alpha = 0.5, outlier.shape = NA) +
  labs(y = expression(paste(italic("STC2")~"(log2)")), 
       x = "Sample Type", 
       title = expression(paste(italic("STC2"), "ESR1 KD in MCF7")))+ 
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, face ="bold"), 
        axis.text.x = element_text(size=10), axis.ticks.x=element_blank(), 
        axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 11, face ="italic"),#Italic if it is a gene. 
        axis.text.y = element_text(size = 10), 
        panel.background = element_rect(fill = "white",
                                        colour = "white"),
        axis.line = element_line(linewidth = 0.7, linetype = "solid",
                                 colour = "black")) +
  scale_fill_manual(values = c("#00CD00","firebrick2"))+ 
  geom_signif(comparisons = my_comparisons,map_signif_level = T, textsize=5, y_position = c(7.0), test = "t.test")
p

# Save as PNG
png("stc2_v2.png", width = 6, height = 6, units = "in", res = 300)  # Set resolution to 300 dpi for PNG
plot(p)
dev.off()

# Save as PDF
pdf("stc2_v2.pdf", width = 6, height = 6)  # PDF works with dimensions in inches
plot(p)
dev.off()


