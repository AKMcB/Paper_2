#######################
## Loading libraries ##
#######################

library(data.table)
library(dplyr)
library(table1)
library(tidyverse)

#############################
## Reading expression file ##
#############################

setwd("C:/Users/abe186/UiT Office 365/O365-Phd Anne - General/Metabric")
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
genes <- rownames_to_column(genes, "Patient ID")

#Upload the new info file 
info <- read.csv2("brca_metabric_clinical_data.csv",sep = ";", as.is = T, check.names = F)
head(info)
unique(info$`Sample Type`)

info <- info[, !names(info) %in% c(
  "Study ID",
  "Sample ID",
  "Breast Surgery",
  "Cancer Type",
  "Cellularity",
  "Chemotherapy",
  "Cohort",
  "Hormone Therapy",
  "Integrative Cluster",
  "Primary Tumor Laterality",
  "Mutation Count",
  "Nottingham prognostic index",
  "Oncotree Code",
  "Overall Survival (Months)",
  "Radio Therapy",
  "Relapse Free Status (Months)",
  "Number of Samples Per Patient",
  "Sample Type",
  "TMB (nonsynonymous)"
)]

info[info == ""] <- NA
str(info)

merged <- merge (info, genes, by = "Patient ID")
merged <- subset(merged, merged$`ER status measured by IHC` %in% "Positve") #764
merged$TRIM45_expression <- ifelse(merged$TRIM45 >= median(merged$TRIM45), 'High', "Low")

merged$`Tumor Stage` <- as.factor(merged$`Tumor Stage`)
merged$`Neoplasm Histologic Grade` <- as.factor(merged$`Neoplasm Histologic Grade`)
##################
# Clinical table #
##################

#check_normality function
check_normality <- function(data) {
  numerical_columns <- sapply(data, is.numeric)
  normality_results <- lapply(names(data)[numerical_columns], function(col_name) {
    col <- data[[col_name]]
    shapiro_test <- shapiro.test(col)
    return(list(
      Column = col_name,
      W = shapiro_test$statistic,
      p_value = shapiro_test$p.value,
      Normal = shapiro_test$p.value > 0.05))
  })
  vars.normal <- names(data)[numerical_columns][sapply(normality_results, function(result) result$Normal)]
  return(vars.normal)
}

vars.normal <- check_normality(merged)

rndr <- function(x, name, ...) {
  cont <- ifelse(name %in% vars.normal, "Mean (SD)", "Median (Min, Max)")
  render.default(x, name, render.continuous=c("", cont), ...)
}


pvalue <- function(x, ...) {
  x <- x[-length(x)] 
  y <- unlist(x)
  g <- factor(rep(1:length(x), times=sapply(x, length))) 
  if (is.numeric(y)) {
    swtest <- shapiro.test(y)$p.value
    vartest <- var.test(y ~ g)$p.value
    if (all(swtest > 0.05, vartest > 0.05)) {
      p <- t.test(y ~ g, var.equal=TRUE)$p.value
      symbol <- "T-Equal"
    } else {
      p <- t.test(y ~ g,var.equal=FALSE)$p.value
      symbol <- "T-Unequal"
    } 
    if (all(swtest <= 0.05)) {
      p <- wilcox.test(y ~ g, exact=FALSE)$p.value
      symbol <- "Wil"
    }} 
  if (is.factor(y)) {
    p <- tryCatch({
      chisq.test(table(y, g))$p.value
    }, warning = function(w) {
      print(paste("Warning message:", w))
      return(NULL)
    })
    if (!is.null(p)) {
      p <- chisq.test(table(y, g))$p.value
      symbol <- "Chi"
    } else
    {
      p <- fisher.test(table(y, g), simulate.p.value=TRUE)$p.value
      symbol <- "Fish"
    }}
  c(paste0(sub("<", "&lt;", format.pval(p, digits=3, eps=0.001)), "<sup>", symbol, "</sup>"))
}


caption <- c("Table 1: Association of clinical parameters 
between ER+ TRIM45 expression groups")

footnote <- c("T-Equal = t.test with equal variance
               T-Unequal = t.test with unequal variance
               ns= not significant
               SD = standard deviation 
               OS = Overall Survival 
               DSS = Death Specific Survival")


head(merged)

table <- table1(~ Sex+
                  `Age at Diagnosis`+
                  `Pam50 + Claudin-low subtype`+
                  `Cancer Type Detailed`+
                  `Tumor Stage`+
                  `Tumor Size` +
                  `Neoplasm Histologic Grade`+
                  `Inferred Menopausal State`+
                  `Patient's Vital Status` +
                  `Overall Survival Status`+
                  `Relapse Free Status`+
                  `PR Status`+
                  `HER2 Status`
                |TRIM45_expression, data = merged,
                extra.col = list(`p-value` = pvalue), caption = caption,
                footnote = footnote, extra.col.pos = 3,
                render = rndr, render.categorical = "FREQ (PCTnoNA%)",
                topclass = " Rtable1-times ")

print(table)

merged$`Tumor Stage` <- factor(merged$`Tumor Stage` ,
                               levels = c("0","1", "2", "3", "4"),
                               labels = c("Stage 0", "Stage I", "Stage II", "Stage III", "Stage IV"))

merged$`Overall Survival Status` <- factor(merged$`Overall Survival Status` ,
                                           levels = c("0:LIVING","1:DECEASED"),
                                           labels = c("Living", "Deceased"))


merged$`Relapse Free Status` <- factor(merged$`Relapse Free Status` ,
                                       levels = c("0:Not Recurred","1:Recurred"),
                                       labels = c("Not Recurred", "Recurred"))

merged$`Pam50 + Claudin-low subtype` <- factor(merged$`Pam50 + Claudin-low subtype` ,
                                               levels = c("LumA","LumB", "Her2", "Basal", "claudin-low", "Normal"),
                                               labels = c("LumA", "LumB", "HER2-enriched", "Basal-like", "Basal-like", "Normal"))


#p-value function
symbols <- c("****", "***", "**", "*", "ns")
cutoffs <- c(0.0001, 0.001, 0.01, 0.05)
lookup_table <- setNames(symbols, cutoffs)

footnote <- c("SD = Standard Deviation &nbsp;
               OS = Overall survival &nbsp;
               DSS = Death Specific Survival &nbsp;
              DFI = Disease Free Interval &nbsp;
              PFI = Progression Free Interval",
              "ns= Not Significant &nbsp
               <0.0001 = ****  &nbsp;
              <0.001 = ***  &nbsp;
              <0.01 = **  &nbsp;
              <=0.05 = *")

#Change the p-value function so that the symbol are not included
pvalue <- function(x, ...) {
  x <- x[-length(x)] 
  y <- unlist(x)
  g <- factor(rep(1:length(x), times=sapply(x, length)))
  if (is.numeric(y)) {
    swtest <- shapiro.test(y)$p.value
    vartest <- var.test(y ~ g)$p.value
    if (all(swtest > 0.05, vartest > 0.05)) {
      p <- t.test(y ~ g, var.equal=TRUE)$p.value
      symbol <- lookup_table[findInterval(p, cutoffs)+1]
    } else {
      p <- t.test(y ~ g, var.equal=FALSE)$p.value
      symbol <- lookup_table[findInterval(p, cutoffs)+1]
    } 
    if (all(swtest <= 0.05)) {
      p <- wilcox.test(y ~ g, exact=FALSE)$p.value
      symbol <- lookup_table[findInterval(p, cutoffs)+1]
    }} 
  if (is.factor(y)) {
    p <- tryCatch({
      chisq.test(table(y, g))$p.value
    }, warning = function(w) {
      print(paste("Warning message:", w))
      return(NULL)
    })
    if (!is.null(p)) {
      p <- chisq.test(table(y, g))$p.value
      symbol <- lookup_table[findInterval(p, cutoffs)+1]
    } else 
    {
      p <- fisher.test(table(y, g), simulate.p.value=TRUE)$p.value
      symbol <- lookup_table[findInterval(p, cutoffs)+1]
    }} 
  c(paste0(sub("<", "&lt;", format.pval(p, digits=3, eps=0.001)), "<sup>", symbol, "</sup>"))
}



table <- table1(~ Sex+
                  `Age at Diagnosis`+
                  `Cancer Type Detailed`+
                  `Tumor Stage`+
                  `Tumor Size` +
                  `Neoplasm Histologic Grade`+
                  `Inferred Menopausal State`+
                  `PR Status`+
                  `HER2 Status`
                |TRIM45_expression, data = merged,
                extra.col = list(`p-value` = pvalue), caption = caption,
                footnote = footnote, extra.col.pos = 3,
                render = rndr, render.categorical = "FREQ (PCTnoNA%)",
                topclass = " Rtable1-times ")

print(table)


#########################
##Write the output file##
#########################
#https://github.com/benjaminrich/table1/issues/33 

write_table1 <- function(x,                   # a table1 object
                         file,                # path to output .pdf file
                         landscape = FALSE,   # landscape print?
                         scale = 1,           # scaling factor
                         width = 8.5,         # width of resulting pdf (in.)
                         height = 15,         # height of resulting pdf (in.)
                         dump_html = TRUE) {  
  
  x <- htmltools::HTML(x)
  src <- system.file(package = "table1", "table1_defaults_1.0")
  default.style <- htmltools::htmlDependency("table1", "1.0",
                                             src = src,
                                             stylesheet = "table1_defaults.css")
  
  x <- htmltools::div(class = "Rtable1", default.style, x)
  x <- htmltools::browsable(x)
  
  htmltools::save_html(x, file = gsub(".pdf", ".html", file))
  pagedown::chrome_print(input = gsub(".pdf", ".html", file),
                         output = file,
                         options = list(landscape = landscape,
                                        scale = scale,
                                        paperWidth = width,
                                        paperHeight = height))
  if(dump_html) unlink(gsub(".pdf", ".html", file))
}

write_table1(table, "2025_03_25_clinical_table_metabric_trim45_high_vs_low.pdf")
