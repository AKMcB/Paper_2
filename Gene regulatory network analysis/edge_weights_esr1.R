#############
# Libraries #
#############
library(tidyverse)
library(data.table)

#########################
# Read in the edge file #
#########################
expr <- as.data.frame(fread("C:/Users/abe186/UiT Office 365/O365-Phd Anne - General/NCMM/TCGA-BRCA/lioness_output/lioness_filtered_for_TFs_ESR1.csv"))
head(expr)[1:5]

#Filter gene of interest 
expr <- subset(expr, expr$Target == "GREB1")

expr <- expr[,-1] #all targets are affected by ESR1
rownames(expr) <- expr[,1]
expr <- expr[,-1]

#Replace the sample with TCGA-IDs
id <- read.csv2("C:/Users/abe186/UiT Office 365/O365-Phd Anne - General/NCMM/TCGA-BRCA/raw_files/tcga_id.csv")
id
id$extra <- "c"
rownames(id) <- id[, 1]

#Change the colnames of indegrees to tgca_id
expr <- setnames(expr, rownames(id))
head(expr)[1:5]

#Load expression data 
expr_1 <- fread("C:/Users/abe186/UiT Office 365/O365-PhD Anne - General/NCMM/TCGA-BRCA/raw_files/tcga_brca_er_positive_expr.csv",sep = ",")
ann1 <- expr_1[,c("id", "ESR1", "TRIM45", "GREB1")]
ann1 <- column_to_rownames(ann1, "id")

#The patients are in columns and genes in rownames
all(rownames(ann1) == colnames(expr))
#all(rownames(clin) == colnames(expr))

test <- expr

test <- as.data.frame(t(test))
test <- rownames_to_column(test, "id")
test$GREB1 <- format(test$GREB1, scientific = FALSE)
test$GREB1 <- as.numeric(test$GREB1)

ann1 <- rownames_to_column(ann1, "id")
colnames(ann1) <- c("id", "ESR1_exp", "TRIM45_exp", "GREB1_exp")

test <- merge(test,ann1 , by = "id")

# Find the position where TRIM45 transitions from negative to positive
separation_position <- which(sort(test$GREB1) >= 0)[1]  # First position where TRIM45 >= 0
# Create the plot
number <- as.data.frame(test$GREB1 >0)

p <-ggplot(test, aes(x = reorder(id, GREB1), y = GREB1, fill = GREB1), color= NA) +
  geom_bar(stat = "identity", )+
  geom_vline(xintercept = separation_position, linetype = "dashed", color = "black", size = 0.8, alpha = 0.8) +  # Add vertical line
  scale_fill_gradient(low = "#00CD00", high = "firebrick2") +  
  labs(title = "Comparison of GREB1 and Expression Values with Separation Line", 
       x = "Samples (Ordered by GREB1)", 
       y = "Value") +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.background = element_rect(fill = "white", colour = "white"),
    plot.background = element_rect(fill = "white", colour = "white"),
    axis.line = element_line(linewidth = 0.7, linetype = "solid",
                             colour = "black"),
    axis.line.x = element_blank())

p

# Save as PNG
png("edge_weight_esr1_greb1.png", width = 8, height = 6, units = "in", res = 300)  # Set resolution to 300 dpi for PNG
plot(p)
dev.off()

# Save as PDF
pdf("edge_weight_esr1_greb1.pdf", width = 6, height = 4)  # PDF works with dimensions in inches
plot(p)
dev.off()


test$group <- ifelse(test$GREB1 >= 0, "Positive", "Negative")
my_comparisons <- list( c("Negative", "Positive")) 
symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
                    symbols = c("****", "***", "**", "*", "ns"))

p <- ggplot(test, aes(x = group, y = GREB1_exp, fill = group))+
  geom_point(alpha=0.5,position = position_jitter(), shape= 21, size= 6)+
  geom_boxplot(fill = "white", alpha = 0.8, outlier.shape = NA) +
  labs(y = expression(paste(italic("ESR1")~"(log2+1)")), 
       x = "Sample Type", 
       title = expression(paste(italic("ESR1"), "expression between ESR1 regulation groups")))+ 
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, face ="bold"), 
        axis.text.x = element_text(size=10), axis.ticks.x=element_blank(), 
        axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 11, face ="italic"),#Italic if it is a gene. 
        axis.text.y = element_text(size = 10), 
        panel.background = element_rect(fill = "white",
                                        colour = "white"),
        axis.line = element_line(linewidth = 0.7, linetype = "solid",
                                 colour = "black")) +
  scale_fill_manual(values = c("#00CD00","firebrick2"))+ 
ggpubr::geom_signif(comparisons = my_comparisons,map_signif_level = T, textsize=5, y_position = c(8.5))

p
    
# Save as PNG
png("esr1_greb1_expr_edge_weight_groups.png", width = 6, height = 6, units = "in", res = 300)  # Set resolution to 300 dpi for PNG
plot(p)
dev.off()

# Save as PDF
pdf("esr1_greb1_expr_edge_weight_groups.pdf", width = 6, height = 6)  # PDF works with dimensions in inches
plot(p)
dev.off()  
