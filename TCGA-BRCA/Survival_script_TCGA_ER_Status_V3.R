#############
# Libraries #
#############
library(tibble)
library(survival)
library(survminer)
library(ggpubr)
library(mclust)

########################
# Read expression data #
########################

setwd("C:/Users/abe186/UiT Office 365/O365-Phd Anne - General/TCGA-BRCA/Kaplan-Meier/Receptor Status TRIM45")
#Upload TCGA file 
ER <- read.csv2("gdc_brca_expr_edit.csv",sep = ";", as.is = T, check.names = F)
ER <- ER[,-2]

#Subset the gene of interest
genes <- subset(ER, ER$id %in% c("ENSG00000134253", "ENSG00000113739"))
rownames(genes) <- genes[,1]
genes <- genes[,-1]
genes <- as.data.frame(t(genes))
colnames(genes) <- c("STC2","TRIM45")

#Upload the new info file 
info <- read.csv2("TCGA_BRCA_Updated_Clinical_Data.csv",sep = ";", as.is = T, check.names = F)
info_rec <- info[,-c(2:44)]

#Merge the genes and info file based on patient ID 
genes <- tibble::rownames_to_column(genes, "ID")
merged <- merge(genes, info_rec, by.x= "ID", by.y = "TCGA_id")
#write.csv2(merged, "TRIM45_receptor status_TCGA_clinical.csv")

#Subset info file to obtain DSS and PFI 
info_dss <- info[,c(1,29,30)]
info_pfi <- info[,c(1, 33,34)]
#Merge the merged file and the info_dss file 
#merged <- tibble::rownames_to_column(merged, "ID")
info_dss <- merge(merged, info_dss, by.x= "ID", by.y = "TCGA_id")

info_pfi <- merge(merged, info_pfi, by.x= "ID", by.y = "TCGA_id")

#Subset ER_DSS and ER_PFI based on positive ER status 
info_dss <- subset(info_dss, info_dss$er_status_by_ihc %in% "Positive")


info_pfi <- subset(info_pfi, info_pfi$er_status_by_ihc %in% "Positive")

################################################################################
#GMM as cutoff
#line 51-106 was created by Silje Berg-Henry (Bachelor student at Tumor Biology research group)
data_dss <- info_dss %>%
  dplyr::pull(STC2) 

fit_dss <- Mclust(data_dss, G = 2) # G is number of clusters 

means <- fit_dss$parameters$mean 
v <- median(data_dss)
shapiro.test(data_dss)
v_dss <- mean(means) 

plot <- ggplot(info_dss, aes(x = STC2)) +
  geom_histogram(aes(y = ..density..), fill = "#9B111E", color = "black", alpha = 0.5, bins = 30) +
  geom_density(color = "black", size = 1) +
  geom_vline(xintercept = v_dss, color = "black", linetype = "dashed", size = 1.5) +  # Add vertical line for v
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

pdf("2025_05_21_S100A9_distribution_tcga.pdf", width = 10, height = 10, onefile = F)
print(plot)
dev.off()

# GMM as cutoff: 
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

################################################################################
#Define TRIM45 expression after subsetting based on receptor status
info_dss$stc2_expression <- ifelse(info_dss$STC2 >= v_dss, 'High', "Low")

info_pfi$stc2_expression <- ifelse(info_pfi$STC2 >= v_pfi, 'High', "Low")

info_dss$TRIM45_expression <- ifelse(info_dss$TRIM45 >= median(info_dss$TRIM45), 'High', "Low")

info_pfi$TRIM45_expression <- ifelse(info_pfi$TRIM45 >= median(info_pfi$TRIM45), 'High', "Low")

#Remove timepoints that start at 0 
info_dss <- subset(info_dss, info_dss$DSS.time > 0)

info_pfi <- subset(info_pfi, info_pfi$PFI.time > 0)

#Convert DSS/PFI.time into years 
info_dss$years <- info_dss$DSS.time/356

info_pfi$years <- info_pfi$PFI.time/365

#Remove NA form the df
info_dss <- info_dss[complete.cases(info_dss),]

info_pfi <- info_pfi[complete.cases(info_pfi),]

#Define survival
survival_dss = Surv(time= info_dss$years, event = info_dss$DSS)

survival_pfi = Surv(time = info_pfi$years, event = info_pfi$PFI)

survival_fit_dss<- survfit(formula = survival_dss ~ info_dss$TRIM45_expression , data = info_dss)

survival_fit_pfi<- survfit(formula = survival_pfi ~ info_pfi$TRIM45_expression , data = info_pfi)

table(TRIM45=info_pfi$TRIM45_expression, stc2=info_pfi$stc2_expression)

p <- ggsurvplot(fit = survival_fit_dss, 
                pval = TRUE, surv.median.line = "none",
                legend = "right",
                title = "TCGA-BRCA ER+",
                xlab = "Years", 
                ylab = "DSS Probability", 
                ylim=c(0.0,1), xlim=c(0,25), 
                pval.coord=c(0,0.37), 
                conf.int = T, risk.table = T, 
                legend.labs=c("High TRIM45", "Low TRIM45"),
                legend.title="Expression Level",
                palette = c("#E69F00","#0072B2","#009E73","#CC79A7"),
                risk.table.y.text = TRUE, risk.table.title="",
                risk.table.height=0.15, 
                risk.table.fontsize=3.5, 
                conf.int.style="step",
                ncensor.plot = TRUE,
                ncensor.plot.height = 0.25,
                break.x.by= 5,
                tables.theme = theme(axis.text.x = element_blank(), axis.ticks.x=element_blank(), 
                                     axis.title.x = element_blank(), axis.title.y = element_blank(),
                                     axis.text.y = element_text(size = 10)),font.x=10, font.y=11, font.tickslab=10)

p


pdf("2025_05_21_trim45_dss_tcga_25.pdf", width = 10, height = 10, onefile = F)
print(p)
dev.off()

p <- ggsurvplot(fit = survival_fit_pfi, 
                pval = TRUE, surv.median.line = "none",
                legend = "right",
                title = "TCGA-BRCA ER+",
                xlab = "Years", 
                ylab = "PFI Probability", 
                ylim=c(0.6,1), xlim=c(0,5.5), 
                pval.coord=c(0,0.62), 
                conf.int = T, risk.table = T, 
                legend.labs=c("High TRIM45", "Low TRIM45"),
                legend.title="Expression Level",
                palette = c("#E69F00","#0072B2","#009E73", "#CC79A7"),
                risk.table.y.text = TRUE, risk.table.title="",
                risk.table.height=0.15, 
                risk.table.fontsize=3.5, 
                conf.int.style="step",
                ncensor.plot = TRUE,
                ncensor.plot.height = 0.25,
                break.x.by= 1,
                tables.theme = theme(axis.text.x = element_blank(), axis.ticks.x=element_blank(), 
                                     axis.title.x = element_blank(), axis.title.y = element_blank(),
                                     axis.text.y = element_text(size = 10)),font.x=10, font.y=11, font.tickslab=10)


p

pdf("2025_05_21_trim45_pfi_tcga_5.pdf", width = 10, height = 10, onefile = F)
print(p)
dev.off()