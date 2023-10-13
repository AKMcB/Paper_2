#Libraries
library(tibble)
library(survival)
library(survminer)
setwd("C:/Users/abe186/UiT Office 365/O365-Bioinformatikk TRIM27 - General/Metabric")
#Upload METABRIC file 
ER <- read.csv2("Ensembl_Metabric_17724_duplicates removed .csv",sep = ";", as.is = T, check.names = F)
ER <- ER[,-c(1,3)]

#Subset the gene of interest
genes <- subset(ER, ER$Ensembl %in% "ENSG00000134253")
rownames(genes) <- genes[,1]
genes <- genes[,-1]
genes <- t(genes)
genes <- as.data.frame(genes)
colnames(genes) <- "TRIM45"

#Upload the new info file 
info <- read.csv2("brca_metabric_clinical_data.csv",sep = ";", as.is = T, check.names = F)
info_rec <- info[,-c(1,3:11)]
info_rec <- info_rec[,-c(2,4,5,7:29)]

#Merge the genes and info file based on patient ID 
genes <- tibble::rownames_to_column(genes, "ID")
merged <- merge(genes, info_rec, by.x= "ID", by.y = "Patient ID")
#write.csv2(merged, "TRIM45_receptor_status_Metabric.csv")
#merged<- read.csv2("TRIM45_receptor status_TCGA_clinical.csv",sep = ";", as.is = T, check.names = F)
#merged <- merged[,-1]

#Subset info file to obtain DSS and PFI 
info_dss <- info[,c(2,26, 27,30,31,39)]
write.csv2(info_dss, "DSS_PFI_info_Metabric.csv")

#merged<- read.csv2("TRIM45_receptor_status_Metabric.csv",sep = ";", as.is = T, check.names = F)
#dss<- read.csv2("DSS_Metabric.csv",sep = ";", as.is = T, check.names = F)


#Merge the merged file and the dss/pfi file 
ER_DSS <- merge(merged, dss, by.x= "ID", by.y = "Patient.ID")

# 1)Subset ER_DSS and ER_PFI based on positive ER status 
ER_DSS_pos <- subset(ER_DSS, ER_DSS$`ER Status` %in% "Positive")

#Define TRIM45 expression after subsetting based on receptor status
ER_DSS_pos$TRIM45_expression <- ifelse(ER_DSS_pos$TRIM45 >= median(ER_DSS_pos$TRIM45), 'High', "Low")

#Remove timepoints that start at 0 
ER_DSS_pos <- subset(ER_DSS_pos, ER_DSS_pos$DSS.months > 0)

#Convert DSS/PFI.time into years 
ER_DSS_pos$years <- ER_DSS_pos$DSS.months/12

#Define survival
survival_DSS_pos = Surv(time= ER_DSS_pos$years, event = ER_DSS_pos$DSS)

survival_fit_DSS_pos<- survfit(formula = survival_DSS_pos ~ ER_DSS_pos$TRIM45_expression, data = ER_DSS_pos)

#Make plot
p <- ggsurvplot(fit = survival_fit_DSS_pos, 
                pval = TRUE, 
                legend = "right",
                title = "ER-Positive Patients",
                xlab = "Years", 
                ylab = "DSS Probability", 
                ylim=c(0.1,1), xlim=c(0,30), 
                pval.coord=c(0,0.12), 
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

pdf("TRIM45_ER_pos_DSS_metabric_2.pdf", width = 8, height = 7, onefile = F)
print(p)
dev.off()


