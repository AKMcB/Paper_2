library(tidyr)
library(tidyverse)
library(ggpubr)
library(RColorBrewer)
install.packages("wesanderson")
library(wesanderson)
install.packages("pals")
library(pals)
pal.bands(alphabet, alphabet2, cols25, glasbey, kelly, polychrome, 
          stepped, tol, watlington,
          show.names=FALSE)
setwd("C:/Users/abe186/UiT Office 365/O365-Bioinformatikk TRIM27 - General/Cancer cell lines cyclopedia/Raw Data")
expr <- read.csv2("CCLE_exp_info_Cancer_Codes.csv", sep = ";", as.is = T, check.names = F)
bcra <- subset(expr, expr$OncotreeCode == "BRCA")
luminal <- subset(bcra, bcra$MolecularSubtype == "luminal")


luminal <- luminal[,c("CellLineName", "TRIM45", "ESR1")]
#rownames(luminal) <- luminal[,1]
#luminal$CellLineName <- NULL

luminal$TRIM45 <- as.numeric(luminal$TRIM45)
luminal$ESR1 <- as.numeric(luminal$ESR1)
colnames(luminal)[1] <- "Cell Lines"
# Create a data frame for plotting


p <- ggscatter(luminal, x = "ESR1", y = "TRIM45", label = "Cell Lines", 
          conf.int =F , shape = 21, fill="Cell Lines", size=3,
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "ESR1", ylab = "TRIM45",cor.coef.size = 4,
          title = paste("TRIM45 vs ESR1"))+
theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0.5, face ="italic"), #Italic if it is a gene. 
        axis.text.x = element_text(size=10), axis.ticks.x=element_blank(), 
        axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 11),
        axis.text.y = element_text(size = 10))+
          scale_fill_manual(values=as.vector(glasbey(15)))


pdf("TRIM45_ESR1_corr_CCLE.pdf", height = 8, width = 10)
print(p)
dev.off()

list <- list(t45, esr1)
ggexport(list, filename = "Correlation_TRIM45_ESR1_luminal_CCLs.pdf", height = 10, width = 12)
ggexport(esr1, filename = "Correlation_ESR1_TRIM45_luminal_CCLs.png",res=200, height = 1800, width = 2100)





