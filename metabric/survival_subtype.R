#Libraries
library(tidyverse)
library(survival)
library(survminer)
library(data.table)
setwd("C:/Users/abe186/UiT Office 365/O365-Phd Anne - General/Metabric")
#Upload METABRIC file 
ER <- fread("Ensembl_Metabric_17724_duplicates removed .csv",sep = ";")
ER <- ER[,-c(1,3)]

#Subset the gene of interest
genes <- subset(ER, ER$Ensembl %in% "ENSG00000134253")
rownames(genes) <- genes[,1]
genes <- genes[,-1]
genes <- t(genes)
genes <- as.data.frame(genes)
colnames(genes) <- "TRIM45"
genes <- rownames_to_column(genes, "id")

#Upload the new info file 
info <- read.csv2("brca_metabric_clinical_data.csv",sep = ";", as.is = T, check.names = F)

info_dss <- info[,c(2,10,26,27,39)]
unique(info_dss$`Patient's Vital Status`)
mapping <- c("Living" = 0, "Died of Disease" = 1, "Died of Other Causes"= 0)
info_dss$dss <- ifelse(is.na(info_dss$`Patient's Vital Status`), NA, mapping[info_dss$`Patient's Vital Status`])


info_pfi <- info[, c(2,10,30,31)]
info_pfi$`Relapse Free Status` <- gsub("\\:.*", "",info_pfi$`Relapse Free Status`)

#merge with expression
info_dss <- merge(info_dss, genes, by.x = "Patient ID", by.y = "id")

info_pfi <- merge(info_pfi, genes, by.x = "Patient ID", by.y = "id")

shapiro.test(info_dss$TRIM45)
hist(info_dss$TRIM45)


info_dss <- subset(info_dss, info_dss$`Pam50 + Claudin-low subtype` == "Normal")
hist(info_dss$TRIM45)
median(info_dss$TRIM45)
mean(info_dss$TRIM45)

info_pfi <- subset(info_pfi, info_pfi$`Pam50 + Claudin-low subtype` == "Normal")
hist(info_pfi$TRIM45)
median(info_pfi$TRIM45)
mean(info_pfi$TRIM45)


#Define TRIM45 expression after subsetting based on receptor status
info_dss$TRIM45_expression <- ifelse(info_dss$TRIM45 >= median(info_dss$TRIM45), 'High', "Low")

info_pfi$TRIM45_expression <- ifelse(info_pfi$TRIM45 >= median(info_pfi$TRIM45), 'High', "Low")

#Remove timepoints that start at 0 
info_dss <- subset(info_dss, info_dss$`Overall Survival (Months)` > 0)
info_pfi <- subset(info_pfi, info_pfi$`Relapse Free Status (Months)` > 0)

#Convert DSS/PFI.time into years 
info_dss$years <- info_dss$`Overall Survival (Months)`/12

info_pfi$years <- info_pfi$`Relapse Free Status (Months)`/12
#Define survival
survival_dss = Surv(time= info_dss$years, event = info_dss$dss)

info_pfi$`Relapse Free Status` <- as.numeric(info_pfi$`Relapse Free Status`)
survival_pfi = Surv(time= info_pfi$years, event = info_pfi$`Relapse Free Status`)

survival_fit_dss<- survfit(formula = survival_dss ~ info_dss$TRIM45_expression, data = info_dss)

survival_fit_pfi<- survfit(formula = survival_pfi ~ info_pfi$TRIM45_expression, data = info_pfi)

p <- ggsurvplot(fit = survival_fit_dss, 
                pval = TRUE, 
                legend = "right",
                title = "Normal-like Patients",
                xlab = "Years", 
                ylab = "DSS Probability", 
                ylim=c(0.0,1), xlim=c(0,28), 
                pval.coord=c(0,0.12), 
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


pdf("2025_03_14_trim45_normal_dss_metabric.pdf", width = 5, height = 5, onefile = F)
print(p)
dev.off()

b <- ggsurvplot(fit = survival_fit_pfi, 
                pval = TRUE, 
                legend = "right",
                title = "Normal-like Patients",
                xlab = "Years", 
                ylab = "PFI Probability", 
                ylim=c(0.0,1), xlim=c(0,28), 
                pval.coord=c(0,0.12), 
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

pdf("2025_03_14_trim45_normal_pfi_metabric.pdf", width = 5, height = 5, onefile = F)
print(b)
dev.off()

