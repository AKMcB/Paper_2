library(tidyr)
library(tidyverse)
library(ggpubr)
#Upload TCGA file 
ER <- read.csv2("gdc_brca_expr_edit.csv",sep = ";", as.is = T, check.names = F)
ER <- ER[,-1]

#Subset the gene of interest
genes <- subset(ER, ER$gene %in% c("TRIM45","ESR1"))
rownames(genes) <- genes[,1]
genes <- genes[,-1]
genes <- as.data.frame(t(genes))

#Upload the new info file 
info <- read.csv2("TCGA_BRCA_Updated_Clinical_Data.csv",sep = ";", as.is = T, check.names = F)
info_rec <- info[,-c(2:44)]

#Merge the genes and infor file based on patient ID 
genes <- tibble::rownames_to_column(genes, "ID")
merged <- merge(genes, info_rec, by.x= "ID", by.y = "TCGA_id")
write.csv2(merged, "TRIM45_receptor status_TCGA_clinical.csv")

# 1)Subset based on positive ER status 
ER <- subset(merged, merged$er_status_by_ihc %in% "Positive")
write.csv2(ER, "T45_ESR1_ER_pos_TCGA.csv")

# 2)Subset based on HER2 negative status 
HER <- subset(ER, ER$her2_status_by_ihc %in% "Negative")
write.csv2(HER, "T45_ESR1_ER_pos_HER_neg_TCGA.csv")

ER <- read.csv2("T45_ESR1_ER_pos_TCGA.csv", as.is = T, check.names = F)

ER_1 <- ggscatter(ER, x = "ESR1", y = "TRIM45",
          title = "ER+ (n = 764)",
          color = "black", shape = 21, size = 1,fill="#D73027",alpha=0.4, # Points color, shape and size
          add = "reg.line",  # Add regressin line
          add.params = list(color = "#D73027", fill = "#D73027"), # Customize reg. line
          conf.int = TRUE, # Add confidence interval
          cor.coef = TRUE,
          cor.coef.size = 5,# Add correlation coefficient. see ?stat_cor
          cor.coeff.args = list(method = "pearson", label.sep = "\n"),
          xlab = "TRIM45", ylab = "ESR1")+ 
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0.5, face ="italic"), #Italic if it is a gene. 
        axis.text.x = element_text(size=10), axis.ticks.x=element_blank(), 
        axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 11),
        axis.text.y = element_text(size = 10))

her_ESR1 <- ggscatter(HER, x = "ESR1", y = "TRIM45",
          title = "ER+ HER2-",
          color = "black", shape = 21, size = 1.5,fill="#D73027",alpha=0.4, # Points color, shape and size
          add = "reg.line",  # Add regressin line
          add.params = list(color = "#D73027", fill = "#D73027"), # Customize reg. line
          conf.int = TRUE, # Add confidence interval
          cor.coef = TRUE,
          cor.coef.size = 5,# Add correlation coefficient. see ?stat_cor
          cor.coeff.args = list(method = "pearson", label.sep = "\n"),
          ylab = "TRIM45", xlab = "ESR1")


pdf("TRIM45_ESR1_Corr_TCGA.pdf", height = 5, width = 5)
print(ER_1)
dev.off()

list <- list(T45,ESR1,ER_TRIM45,ER_ESR1,her_t45,her_ESR1)
ggexport(ER_1, filename = "TRIM45_ESR1_Pearson_TCGA_V2.png", height = 1204, width = 1204)
dev.off()

expr$`O-65b`
ggscatter(expr, x = "O-168a", y = "O-168b", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "O-168 Replicate A", ylab = "O-168 Replicate B")



library(Hmisc)


flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}
e <- as.matrix(lin_expr)
e <- data.frame(lin_expr)
res2<-rcorr(as.matrix(lin_expr[,1:2]))
x <- flattenCorrMatrix(res2$r, res2$P)
x
write.csv(x, "Correlation between 0-65 replicates.csv")
library(corrplot)
res <- cor(expr)
corrplot(res, type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45)

shapiro.test(lin_expr$ZEB1)

ggqqplot(expr$`O-65b`)
