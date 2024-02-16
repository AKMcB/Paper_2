#Libraries
library(tibble)
library(survival)
library(survminer)
setwd("C:/Users/abe186/UiT Office 365/O365-Bioinformatikk TRIM27 - General/TCGA-BRCA/Kaplan-Meier/Receptor Status TRIM45")
#Upload TCGA file 
ER <- read.csv2("gdc_brca_expr_edit.csv",sep = ";", as.is = T, check.names = F)
ER <- ER[,-2]

#Subset the gene of interest
genes <- subset(ER, ER$id %in% "ENSG00000134253")
rownames(genes) <- genes[,1]
genes <- genes[,-1]
genes <- t(genes)
genes <- as.data.frame(genes)
colnames(genes) <- "TRIM45"

#Upload the new info file 
info <- read.csv2("TCGA_BRCA_Updated_Clinical_Data.csv",sep = ";", as.is = T, check.names = F)
info_rec <- info[,-c(2:44)]

#Merge the genes and infor file based on patient ID 
genes <- tibble::rownames_to_column(genes, "ID")
merged <- merge(genes, info_rec, by.x= "ID", by.y = "TCGA_id")
#write.csv2(merged, "TRIM45_receptor status_TCGA_clinical.csv")
#merged<- read.csv2("TRIM45_receptor status_TCGA_clinical.csv",sep = ";", as.is = T, check.names = F)
#merged <- merged[,-1]
#rownames(merged) <- merged[,1]
#merged <- merged[,-1]


#Subset info file to obtain DSS and PFI 
info_dss <- info[,c(1,29,30,33,34)]

#Merge the merged file and the info_dss file 
merged <- tibble::rownames_to_column(merged, "ID")
ER <- merge(merged, info_dss, by.x= "ID", by.y = "TCGA_id")


#Subset df based on either DSS or PFI
ER_DSS <- ER[,-c(9,10)]


#ER_DSS_pos <- subset(ER_DSS, ER_DSS$DSS.time > 0)

# 1)Subset ER_DSS and ER_PFI based on positive ER status 
ER_DSS_pos <- subset(ER_DSS, ER_DSS$er_status_by_ihc %in% "Positive")

#Define TRIM45 expression after subsetting based on receptor status
ER_DSS_pos$TRIM45_expression <- ifelse(ER_DSS_pos$TRIM45 >= median(ER_DSS_pos$TRIM45), 'High', "Low")

#Remove timepoints that start at 0 
ER_DSS_pos <- subset(ER_DSS_pos, ER_DSS_pos$DSS.time > 0)

#Convert DSS/PFI.time into years 
ER_DSS_pos$years <- ER_DSS_pos$DSS.time/356

#Remove NA form the df
ER_DSS_pos <- ER_DSS_pos[complete.cases(ER_DSS_pos),]

#Define survival
survival_DSS_pos = Surv(time= ER_DSS_pos$years, event = ER_DSS_pos$DSS)

survival_fit_DSS_pos<- survfit(formula = survival_DSS_pos ~ ER_DSS_pos$TRIM45_expression, data = ER_DSS_pos)



p <- ggsurvplot(fit = survival_fit_DSS_pos, 
           pval = TRUE, 
           legend = "right",
           title = "ER-Positive Patients",
           xlab = "Years", 
           ylab = "DSS Probability", 
           ylim=c(0.3,1), xlim=c(0,25), 
           pval.coord=c(0,0.32), 
           conf.int = T, risk.table = T, 
           legend.labs=c("High TRIM45", "Low TRIM45"),
           legend.title="TRIM45 Expression",
           palette = c("#E69F00","#0072B2"),
           risk.table.y.text = TRUE, risk.table.title="",
           risk.table.height=0.15, 
           risk.table.fontsize=3.5, 
           tables.theme = theme(axis.text.x = element_blank(), axis.ticks.x=element_blank(), 
                                axis.title.x = element_blank(), axis.title.y = element_blank(),
                                axis.text.y = element_text(size = 10)),font.x=10, font.y=11, font.tickslab=10)

p

pdf("TRIM45_ER_pos_DSS_TCGA_2.pdf", width = 8, height = 7, onefile = F)
print(p)
dev.off()




