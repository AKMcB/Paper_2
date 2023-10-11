#Download libraries
library(tibble)
library(ggpubr)
#-------------------------------------------------------------------------------
setwd("C:/Users/abe186/UiT Office 365/O365-Bioinformatikk TRIM27 - General/TCGA-BRCA/Raw data")
expr <- read.csv2("TMM_TCGA_BRCA_counts_log2.csv")
rownames(expr) <- expr$X
expr <- expr[,-1]

expr <- as.data.frame(t(expr))
expr <- tibble::rownames_to_column(expr, "id")
expr$id <- gsub("\\.","-", expr$id)


#Subset the df based on normal samples -11
#grepl is to search for a specific pattern within a character vector or string
normal<- subset(expr, grepl("-11", expr$id))

#Save file 
#write.csv2(normal, "TMM_TCGA_BRCA_counts_log2_only_normal.csv")

#create a df with just the tumor samples 
tumor <- expr[!grepl("-11", expr$id), ]

#Save file 
#write.csv2(tumor, "TMM_TCGA_BRCA_counts_log2_only_tumor.csv")

#Subset the dfs based on gene of interest
rownames(tumor) <- tumor[,1]
tumor$id <- NULL
tumor <- as.data.frame(t(tumor))

tumor <- tibble::rownames_to_column(tumor, "Gene")
trim_tum <- subset(tumor, tumor$Gene == "ENSG00000134253")#TRIM45 
#trim_tum <- subset(tumor, tumor$Gene == "ENSG00000091831")#ESR1

#trim_tum <- subset(tumor, tumor$Gene == "ENSG00000119888")#EPCAM
# Create a new row with the same value for all variables
new_row <- rep("Tumor", 1105)

# Add the new row to the dataframe using rbind()
trim_tum <- rbind(trim_tum, new_row)


rownames(normal) <- normal[,1]
normal$id <- NULL
normal <- as.data.frame(t(normal))

normal <- tibble::rownames_to_column(normal, "Gene")
trim_nom <- subset(normal, normal$Gene == "ENSG00000134253")#TRIM45
#trim_nom <- subset(normal, normal$Gene == "ENSG00000091831")#ESR1
#trim_nom <- subset(normal, normal$Gene == "ENSG00000119888")#EPCAM
# Create a new row with the same value for all variables
new_row <- rep("Normal", 114)

# Add the new row to the dataframe using rbind()
trim_nom <- rbind(trim_nom, new_row)


trim_nom <- trim_nom[,-1]
test <- cbind(trim_tum, trim_nom)

#-------------------------------------------------------------------------------
#Making boxplot

rownames(test) <- test[,1]
test$Gene <- NULL
test <- as.data.frame(t(test))
colnames(test)[1] <- "TRIM45"
colnames(test)[2] <- "Sample"
test$TRIM45 <- as.numeric(test$TRIM45)
str(test)

my_comparisons <- list( c("Tumor", "Normal") ) 
symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
                    symbols = c("****", "***", "**", "*", "ns"))



p <- ggboxplot(test, x="Sample", y="TRIM45",outlier.shape = NA,
                    show.legend=F,
                    palette = c("#00CD00", "#EE2C2C"), 
                    order = c("Normal", "Tumor"),
                    ylab = "TRIM45 Expression",
               ggtheme = theme_pubr(legend = "right")) +
        theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, face ="italic"), #Italic if it is a gene. 
        axis.text.x = element_text(size=10), axis.ticks.x=element_blank(), 
        axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10),
        axis.text.y = element_text(size = 10))+
  geom_signif(comparisons = my_comparisons,map_signif_level = T,y_position = c(6.5), textsize = 10)+ 
  geom_jitter(shape = 21,stroke=0.1,size=2, aes(fill= Sample, alpha=0.7), position = position_jitter(width = 0.2, height = 0.5))



p
pdf("TRIM45_exp_normal_vs_tumor.pdf", height = 5, width = 5)
print(p)
dev.off()


ggexport(p, filename = "TRIM45_exp_normal_vs_tumor_BRCA.pdf", height = 15, width = 15)
ggexport(p, filename = "TRIM45_exp_normal_vs_tumor_BRCA.png",res=200, height = 2000, width = 2000)
?ggexport
rm(p)




