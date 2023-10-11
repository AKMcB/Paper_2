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
write.csv2(merged, "TRIM45_receptor status_TCGA_clinical.csv")
#merged<- read.csv2("TRIM45_receptor status_TCGA_clinical.csv",sep = ";", as.is = T, check.names = F)
#merged <- merged[,-1]
rownames(merged) <- merged[,1]
merged <- merged[,-1]


#Subset info file to obtain DSS and PFI 
info_dss <- info[,c(1,29,30,33,34)]

#Merge the merged file and the info_dss file 
ER <- merge(merged, info_dss, by.x= "ID", by.y = "TCGA_id")
write.csv2(ER, "TRIM45_receptor status_TCGA_clinical_DSS_PFI_info.csv")

ER<- read.csv2("TRIM45_receptor status_TCGA_clinical_DSS_PFI_info.csv",sep = ";", as.is = T, check.names = F)
ER <- ER[,-1]
rownames(ER) <- ER[,1]
ER <- ER[,-1]


#Subset df based on either DSS or PFI
ER_DSS <- ER[,-c(3,7,8)]
ER_PFI <- ER[,-c(3,5,6)]

ER_DSS_pos <- subset(ER_DSS, ER_DSS$DSS.time > 0)
ER_DSS_her <- subset(ER_DSS, ER_DSS$DSS.time > 0)

ER_PFI_pos <- subset(ER_PFI, ER_PFI$PFI.time > 0)
ER_PFI_her <- subset(ER_PFI, ER_PFI$PFI.time > 0)

# 1)Subset ER_DSS and ER_PFI based on positive ER status 
ER_DSS_pos <- subset(ER_DSS, ER_DSS$er_status_by_ihc %in% "Positive")
ER_PFI_pos <- subset(ER_PFI, ER_PFI$er_status_by_ihc %in% "Positive")


# 2) You can also subset the df based on HER2 negative status 
ER_DSS_her <- subset(ER_DSS_pos, ER_DSS_pos$her2_status_by_ihc %in% "Negative")
ER_PFI_her <- subset(ER_PFI_pos, ER_PFI_pos$her2_status_by_ihc %in% "Negative")


#Define TRIM45 expression after subsetting based on receptor status
ER_DSS_pos$TRIM45_expression <- ifelse(ER_DSS_pos$TRIM45 >= median(ER_DSS_pos$TRIM45), 'High', "Low")
ER_DSS_her$TRIM45_expression <- ifelse(ER_DSS_her$TRIM45 >= median(ER_DSS_her$TRIM45), 'High', "Low")

ER_PFI_pos$TRIM45_expression <- ifelse(ER_PFI_pos$TRIM45 >= median(ER_PFI_pos$TRIM45), 'High', "Low")
ER_PFI_her$TRIM45_expression <- ifelse(ER_PFI_her$TRIM45 >= median(ER_PFI_her$TRIM45), 'High', "Low")


#Remove timepoints that start at 0 
#ER_DSS_pos <- subset(ER_DSS_pos, ER_DSS_pos$DSS.time > 0)
#ER_DSS_her <- subset(ER_DSS_her, ER_DSS_her$DSS.time > 0)

#ER_PFI_pos <- subset(ER_PFI_pos, ER_PFI_pos$PFI.time > 0)
#ER_PFI_her <- subset(ER_PFI_her, ER_PFI_her$PFI.time > 0)


#Convert DSS/PFI.time into years 
ER_DSS_pos$years <- ER_DSS_pos$DSS.time/356
ER_DSS_her$years <- ER_DSS_her$DSS.time/356 

ER_PFI_pos$years <- ER_PFI_pos$PFI.time/365
ER_PFI_her$years <- ER_PFI_her$PFI.time/365


#Remove NA form the df
ER_DSS_pos <- ER_DSS_pos[complete.cases(ER_DSS_pos),]
ER_DSS_her <- ER_DSS_her[complete.cases(ER_DSS_her),]

ER_PFI_pos <- ER_PFI_pos[complete.cases(ER_PFI_pos),]
ER_PFI_her <- ER_PFI_her[complete.cases(ER_PFI_her),]


#Define survival
survival_DSS_pos = Surv(time= ER_DSS_pos$years, event = ER_DSS_pos$DSS)
survival_DSS_her = Surv(time= ER_DSS_her$years, event = ER_DSS_her$DSS)

survival_PFI_pos = Surv(time = ER_PFI_pos$years, event= ER_PFI_pos$PFI)
survival_PFI_her = Surv(time = ER_PFI_her$years, event= ER_PFI_her$PFI)


survival_fit_DSS_pos<- survfit(formula = survival_DSS_pos ~ ER_DSS_pos$TRIM45_expression, data = ER_DSS_pos)
survival_fit_DSS_her<- survfit(formula = survival_DSS_her ~ ER_DSS_her$TRIM45_expression, data = ER_DSS_her)

survival_fit_PFI_pos<- survfit(formula = survival_PFI_pos ~ ER_PFI_pos$TRIM45_expression, data = ER_PFI_pos)
survival_fit_PFI_her<- survfit(formula = survival_PFI_her ~ ER_PFI_her$TRIM45_expression, data = ER_PFI_her)


ggsurv <- ggsurvplot(fit= survival_fit_DSS_pos, 
           pval = TRUE, 
           surv.median.line = "hv", 
           xlab = "DSS (Years)", 
           ylab = "DSS Probability",
           ylim=c(0.0,1), 
           xlim=c(0,25),
           palette = c("#E69F00","#0072B2", "#009E73", "#CC79A7"),
           pval.coord=c(0.1,0.1), 
           break.x.by= 5,         
           conf.int = T, 
           conf.int.alpha = c(0.3), 
           conf.int.style="ribbon",
           risk.table = T,
           risk.table.height = 0.15,
           ncensor.plot = F,
           ncensor.plot.height = 0.20,
           #legend = c(0.5, 0.95),
           legend.labs=c("High TRIM45","Low TRIM45"),
           legend.title= "ER+", 
           font.x=c("12"), 
           font.y= c("12"), 
           font.legend=c("12"), 
           pval.font= c("12"))
#?ggsurvplot
ggsurv$plot <- ggsurv$plot+
  theme(axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size=12))

ggsurv$table <- ggsurv$table+
  theme(plot.title = element_text(size = 12),
        axis.text.x = element_text(size=12),
        axis.title.x.bottom = element_text(size=12), 
        axis.title.y.left = element_text(size=12), 
        legend.title = element_text(size=12),
        legend.text = element_text(size=12))


ggsurv$ncensor.plot <-ggsurv$ncensor.plot+
  theme( plot.title = element_text(size = 12),
         axis.text.x = element_text(size=12),
         axis.text.y = element_text(size=12), 
         axis.title.x.bottom = element_text(size=12),
         axis.title.y.left = element_text(size=12))

ggsurv

pdf("TRIM45_ER_pos_DSS_TCGA.pdf", height = 8, width = 8)
print(ggsurv)
dev.off()

pfi_her <- ggsurvplot(fit= survival_fit_PFI_her, 
           pval = TRUE, 
           surv.median.line = "hv", 
           xlab = "PFI (Years)", 
           ylab = "PFI Probability",
           ylim=c(0.25,1), 
           xlim=c(0,15),
           palette = c("#E69F00","#0072B2", "#009E73", "#CC79A7"),
           pval.coord=c(0.1,0.27), 
           break.x.by= 5,         
           conf.int = T, 
           conf.int.alpha = c(1), 
           conf.int.style="step",
           risk.table = T,
           risk.table.height = 0.15,
           ncensor.plot = TRUE,
           ncensor.plot.height = 0.15,
           #legend = c(0.5, 0.95),
           legend.labs=c("High TRIM45", "Low TRIM45"),
           legend.title= "ER+ HER2-")




list <- list(dss_er, dss_er_5, dss_her, pfi_er, pfi_er_5, pfi_her)
ggexport(list, filename = "TRIM45_ER_pos_HER_neg_DSS_PFI_V4.pdf", height= 15, width = 15)

