#############
# Libraries #
#############
library(tidyverse)
library(survival)
library(survminer)
library(data.table)

########################
# Read expression data #
########################

setwd("C:/Users/abe186/UiT Office 365/O365-Phd Anne - General/Metabric")
#Upload METABRIC file 
ER <- fread("Ensembl_Metabric_17724_duplicates removed .csv",sep = ";")
ER <- ER[,-c(1,3)]

#Subset the gene of interest
genes <- subset(ER, ER$Ensembl %in% c("ENSG00000134253"))
genes <- column_to_rownames(genes, "Ensembl")

#rownames(genes) <- genes[,1]
#genes <- genes[,-1]
genes <- t(genes)
genes <- as.data.frame(genes)
colnames(genes) <- c("TRIM45")
genes <- rownames_to_column(genes, "id")

#Upload the new info file 
info <- read.csv2("brca_metabric_clinical_data.csv",sep = ";", as.is = T, check.names = F)

info_dss <- info[,c(2,10,12,26,27,39)]
unique(info_dss$`Patient's Vital Status`)
mapping <- c("Living" = 0, "Died of Disease" = 1, "Died of Other Causes"= 0)
info_dss$dss <- ifelse(is.na(info_dss$`Patient's Vital Status`), NA, mapping[info_dss$`Patient's Vital Status`])


info_pfi <- info[, c(2,10, 12,30,31)]
info_pfi$`Relapse Free Status` <- gsub("\\:.*", "",info_pfi$`Relapse Free Status`)

#merge with expression
info_dss <- merge(info_dss, genes, by.x = "Patient ID", by.y = "id")

info_pfi <- merge(info_pfi, genes, by.x = "Patient ID", by.y = "id")

#subset based on ER IHC
info_dss <- subset(info_dss, info_dss$`ER status measured by IHC` == "Positve")

info_pfi <- subset(info_pfi, info_pfi$`ER status measured by IHC` == "Positve")

################################################################################
#GMM as cutoff
#line 52-106 was created by Silje Berg-Henry (Bachelor student at Tumor Biology research group)
data_dss <- info_dss %>%
  dplyr::pull(S100A9) 

fit_dss <- Mclust(data_dss, G = 2) #G is number of clusters 

means <- fit_dss$parameters$mean 
v <- median(data_dss)
shapiro.test(data_dss)
v_dss <- mean(means) 

plot <- ggplot(info_dss, aes(x = S100A9)) +
  geom_histogram(aes(y = ..density..), fill = "#9B111E", color = "black", alpha = 0.5, bins = 30) +
  geom_density(color = "black", size = 1) +
  geom_vline(xintercept = v_dss, color = "black", linetype = "dashed", size = 1.5) + 
  labs(y = "Frequency", x = "Expression (log2+1)") +
  theme(legend.position = "none",
        plot.title = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 11),
        axis.title.y = element_text(size = 11),
        panel.background = element_rect(fill = "white", colour = "white"),
        axis.line = element_line(linewidth = 0.7, linetype = "solid", colour = "black"))

plot
pdf("2025_05_21_S100A9_distribution_meta.pdf", width = 10, height = 10, onefile = F)
print(plot)
dev.off()


data_pfi <- info_pfi %>%
  dplyr::pull(STC2) 

fit_pfi <- Mclust(data_pfi, G = 2) #G is number of clusters 

means <- fit_pfi$parameters$mean 

v_pfi <- mean(means)

ggplot(info_pfi, aes(x = STC2)) +
  geom_histogram(aes(y = ..density..), fill = "#9B111E", color = "black", alpha = 0.5, bins = 30) +
  geom_density(color = "#40E0D0", size = 1) +
  geom_vline(xintercept = v_pfi, color = "#40E0D0", linetype = "dashed", size = 1.5) + 
  labs(y = "Frequency", x = "Expression (log2+1)") +
  theme(legend.position = "none",
        plot.title = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 11),
        axis.title.y = element_text(size = 11),
        panel.background = element_rect(fill = "white", colour = "white"),
        axis.line = element_line(linewidth = 0.7, linetype = "solid", colour = "black"))

#################################################################################

#create trim45 level column 
info_dss$stc2_expression <- ifelse(info_dss$STC2 >= v_dss, 'High', "Low")

info_pfi$stc2_expression <- ifelse(info_pfi$STC2 >= v_pfi, 'High', "Low")

info_dss$TRIM45_expression <- ifelse(info_dss$TRIM45 >= median(info_dss$TRIM45), 'High', "Low")

info_pfi$TRIM45_expression <- ifelse(info_pfi$TRIM45 >= median(info_pfi$TRIM45), 'High', "Low")

#Remove timepoints that start at 0 
info_dss <- subset(info_dss, info_dss$`Overall Survival (Months)` > 0) #1445
info_pfi <- subset(info_pfi, info_pfi$`Relapse Free Status (Months)` > 0) #1445

#Convert DSS/PFI.time into years 
info_dss$years <- info_dss$`Overall Survival (Months)`/12

info_pfi$years <- info_pfi$`Relapse Free Status (Months)`/12
#Define survival
survival_dss = Surv(time= info_dss$years, event = info_dss$dss)

info_pfi$`Relapse Free Status` <- as.numeric(info_pfi$`Relapse Free Status`)
survival_pfi = Surv(time= info_pfi$years, event = info_pfi$`Relapse Free Status`)

survival_fit_dss<- survfit(formula = survival_dss ~ info_dss$TRIM45_expression , data = info_dss)

survival_fit_pfi<- survfit(formula = survival_pfi ~ info_pfi$TRIM45_expression , data = info_pfi)

table(TRIM45=info_pfi$TRIM45_expression, stc2=info_pfi$stc2_expression)

p <- ggsurvplot(fit = survival_fit_dss, 
                pval = TRUE, 
                legend = "right",
                title = "METABRIC ER+",
                xlab = "Years", 
                ylab = "DSS Probability", 
                ylim=c(0.0,1), xlim=c(0,30), 
                pval.coord=c(0,0.12), 
                conf.int = T, risk.table = T, 
                legend.labs=c("High TRIM45", "Low TRIM45"),
                legend.title="Expression Level",
                palette = c("#E69F00","#0072B2","#009E73","#CC79A7"),
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


pdf("2025_05_06_stc2_trim45_dss_metabric.pdf", width = 5, height = 5, onefile = F)
print(p)
dev.off()

p <- ggsurvplot(fit = survival_fit_pfi, 
                pval = TRUE, 
                legend = "right",
                title = "METABRIC ER+",
                xlab = "Years", 
                ylab = "PFI Probability", 
                ylim=c(0.0,1), xlim=c(0,29), 
                pval.coord=c(0,0.12), 
                conf.int = T, risk.table = T, 
                legend.labs=c("High TRIM45", "Low TRIM45"),
                legend.title="Expression Level",
                palette = c("#E69F00","#0072B2","#009E73","#CC79A7"),
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

pdf("2025_05_21_stc2_trim45_pfi_metabric.pdf", width = 5, height = 5, onefile = F)
print(p)
dev.off()