#############
# Libraries #
#############
library(tidyverse)
library(data.table)
library(survival)
library(survminer)
setwd("C:/Users/abe186/UiT Office 365/O365-Phd Anne - General/TCGA-BRCA/Kaplan-Meier/Receptor Status TRIM45")

########################
# Read expression file #
########################
ER <- fread("gdc_brca_expr_edit.csv",sep = ";")
ER <- ER[,-2]

dup<-ER[duplicated(ER$id),]

#Subset the gene of interest
genes <- subset(ER, ER$id %in% "ENSG00000134253")
rownames(genes) <- genes[,1]
genes <- genes[,-1]
genes <- as.data.frame(t(genes))
colnames(genes) <- "TRIM45"

genes <- tibble::rownames_to_column(genes, "ID")
dup <- genes[duplicated(genes$ID),]

#Upload the new info file 
info <- read.csv2("TCGA_BRCA_Updated_Clinical_Data.csv",sep = ";", as.is = T, check.names = F)
info_rec <- info[,-c(2:44)]

#Subset info file to obtain DSS and PFI 
info_dss <- info[,c(1,3,29,30)]

info_pfi <- info[, c(1,3,33, 34)]
#Merge the merged file and the info_dss file 
#merged <- tibble::rownames_to_column(merged, "ID")
info_dss <- merge(genes, info_dss, by.x= "ID", by.y = "TCGA_id")

info_pfi <- merge(genes, info_pfi, by.x = "ID", by.y = "TCGA_id")

shapiro.test(info_dss$TRIM45)
hist(info_dss$TRIM45)

info_dss <- subset(info_dss, info_dss$BRCA_Subtype_PAM50 == "Normal")
hist(info_dss$TRIM45)
median(info_dss$TRIM45)
mean(info_dss$TRIM45)

info_pfi <- subset(info_pfi, info_pfi$BRCA_Subtype_PAM50 == "Normal")
hist(info_pfi$TRIM45)
median(info_pfi$TRIM45)
mean(info_pfi$TRIM45)

#Define TRIM45 expression after subsetting based on receptor status
info_dss$TRIM45_expression <- ifelse(info_dss$TRIM45 >= median(info_dss$TRIM45), 'High', "Low")

info_pfi$TRIM45_expression <- ifelse(info_pfi$TRIM45 >= median(info_pfi$TRIM45), 'High', "Low")

#Remove timepoints that start at 0 
info_dss <- subset(info_dss, info_dss$DSS.time > 0) 

info_pfi <- subset(info_pfi, info_pfi$PFI.time > 0) 

#Convert DSS/PFI.time into years 
info_dss$years <- info_dss$DSS.time/356

info_pfi$years <- info_pfi$PFI.time/356

#Remove NA form the df
info_dss <- info_dss[complete.cases(info_dss),] #1024 -> 1006

info_pfi <- info_pfi[complete.cases(info_pfi),]

#Define survival
survival_dss = Surv(time= info_dss$years, event = info_dss$DSS)

survival_pfi = Surv(time = info_pfi$years, event = info_pfi$PFI)

survival_fit_dss<- survfit(formula = survival_dss ~ info_dss$TRIM45_expression, data = info_dss)

survival_fit_pfi<- survfit(formula = survival_pfi ~ info_pfi$TRIM45_expression, data = info_pfi)


p <- ggsurvplot(fit = survival_fit_dss, 
                pval = TRUE, 
                legend = "right",
                title = "Normal-like Patients",
                xlab = "Years", 
                ylab = "DSS Probability", 
                ylim=c(0.7,1), xlim=c(0,10), 
                pval.coord=c(0,0.72), 
                conf.int = T, risk.table = T, 
                legend.labs=c("High", "Low"),
                legend.title="TRIM45 Expression",
                palette = c("#E69F00","#0072B2"),
                risk.table.y.text = TRUE, risk.table.title="",
                risk.table.height=0.15, 
                risk.table.fontsize=3.5, 
                ncensor.plot = TRUE,
                ncensor.plot.height = 0.25,
                break.x.by= 5,
                tables.theme = theme(axis.text.x = element_blank(), axis.ticks.x=element_blank(), 
                                     axis.title.x = element_blank(), axis.title.y = element_blank(),
                                     axis.text.y = element_text(size = 10)),font.x=10, font.y=11, font.tickslab=10)

p


pdf("2025_03_14_trim45_normal_dss.pdf", width = 5, height = 5, onefile = F)
print(p)
dev.off()

b <- ggsurvplot(fit = survival_fit_pfi, 
                pval = TRUE, 
                legend = "right",
                title = "Normal-like Patients",
                xlab = "Years", 
                ylab = "PFI Probability", 
                ylim=c(0.35,1), xlim=c(0,12), 
                pval.coord=c(0,0.41), 
                conf.int = T, risk.table = T, 
                legend.labs=c("High", "Low"),
                legend.title="TRIM45 Expression",
                palette = c("#E69F00","#0072B2"),
                risk.table.y.text = TRUE, risk.table.title="",
                risk.table.height=0.15, 
                risk.table.fontsize=3.5, 
                ncensor.plot = TRUE,
                ncensor.plot.height = 0.25,
                break.x.by= 5,
                tables.theme = theme(axis.text.x = element_blank(), axis.ticks.x=element_blank(), 
                                     axis.title.x = element_blank(), axis.title.y = element_blank(),
                                     axis.text.y = element_text(size = 10)),font.x=10, font.y=11, font.tickslab=10)

b

pdf("2025_03_14_trim45_normal_pfi.pdf", width = 5, height = 5, onefile = F)
print(b)
dev.off()

