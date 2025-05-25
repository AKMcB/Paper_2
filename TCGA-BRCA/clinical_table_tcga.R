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

ER <- read.csv2("gdc_brca_expr_edit.csv",sep = ";", as.is = T, check.names = F)
ER <- ER[,-2]

#Subset the gene of interest
genes <- subset(ER, ER$id %in% "ENSG00000134253")
rownames(genes) <- genes[,1]
genes <- genes[,-1]
genes <- as.data.frame(t(genes))
colnames(genes) <- "TRIM45"
genes <- rownames_to_column(genes, "TCGA_id")

#Upload the new info file 
info <- read.csv2("TCGA_BRCA_Updated_Clinical_Data.csv",sep = ";", as.is = T, check.names = F)
head(info)

info <- info[, !names(info) %in% c(
  "type",
  "race",
  "clinical_stage",
  "histological_grade",
  "initial_pathologic_dx_year",
  "birth_days_to",
  "last_contact_days_to",
  "death_days_to",
  "cause_of_death",
  "new_tumor_event_dx_days_to",
  "treatment_outcome_first_course",
  "margin_status",
  "residual_tumor",
  "OS.time",
  "DSS.time",
  "DFI.time",
  "PFI.time",
  "CNV.Clusters",
  "Mutation.Clusters",
  "DNA.Methylation.Clusters",
  "mRNA.Clusters",
  "miRNA.Clusters",
  "Protein.Clusters",
  "PARADIGM.Clusters",
  "Pan.Gyn.Clusters",
  "Redaction", 
  "subtype"
)]

info[info == ""] <- NA
str(info)

merged <- merge (info, genes, by = "TCGA_id")
merged <- subset(merged, merged$er_status_by_ihc %in% "Positive") #764
merged$TRIM45_expression <- ifelse(merged$TRIM45 >= median(merged$TRIM45), 'High', "Low")

merged <- merged %>%
  mutate(new_tumor_site_comb = if_else(str_trim(new_tumor_event_type) == "" |  new_tumor_event_site == "Other, specify",
                              str_trim(new_tumor_event_type) %>% paste(str_trim(new_tumor_event_site_other), sep = " "),
                              str_trim(new_tumor_event_type) %>% paste(str_trim(new_tumor_event_site), sep = " ")))

##################
# Clinical table #
##################

# Define the check_normality function
check_normality <- function(data) {
  # Step 2: Identify numerical columns
  numerical_columns <- sapply(data, is.numeric)
  
  # Step 3: Perform Shapiro-Wilk test on each numerical column
  normality_results <- lapply(names(data)[numerical_columns], function(col_name) {
    col <- data[[col_name]]
    shapiro_test <- shapiro.test(col)
    return(list(
      Column = col_name,
      W = shapiro_test$statistic,
      p_value = shapiro_test$p.value,
      Normal = shapiro_test$p.value > 0.05
    ))
  })
  
  # Step 4: Identify normally distributed columns
  vars.normal <- names(data)[numerical_columns][sapply(normality_results, function(result) result$Normal)]
  
  return(vars.normal)
}

# Run the normality test to determine normally distributed columns
vars.normal <- check_normality(merged)

rndr <- function(x, name, ...) {
  cont <- ifelse(name %in% vars.normal, "Mean (SD)", "Median (Min, Max)")
  render.default(x, name, render.continuous=c("", cont), ...)
}


pvalue <- function(x, ...) {
  x <- x[-length(x)] # Remove "overall" group
  # Construct vectors of data y, and groups (strata) g
  y <- unlist(x)
  g <- factor(rep(1:length(x), times=sapply(x, length))) # nolint
  if (is.numeric(y)) {
    # check for normality using Shapiro-Wilk test and variance equality using var.test
    swtest <- shapiro.test(y)$p.value
    vartest <- var.test(y ~ g)$p.value
    if (all(swtest > 0.05, vartest > 0.05)) {
      # If both p-values are >0.05 we assume normality and use t.test with equal variances
      p <- t.test(y ~ g, var.equal=TRUE)$p.value
      symbol <- "T-Equal"
    } else {
      # Otherwise, we assume normality with unequal variances and use t.test with var.equal=F
      p <- t.test(y ~ g,var.equal=FALSE)$p.value
      symbol <- "T-Unequal"
    } 
    if (all(swtest <= 0.05)) {
      # Otherwise, reject normality and use Wilcoxon test
      p <- wilcox.test(y ~ g, exact=FALSE)$p.value
      symbol <- "Wil"
    }} 
  if (is.factor(y)) {
    # For categorical variables, perform a chi-squared test of independence
    p <- tryCatch({
      # Si warning chi-test, rC)aliser un Fisher
      chisq.test(table(y, g))$p.value
    }, warning = function(w) {
      # If a warning appears (due to expected cell counts <5 or too few observations), use Fisher's Exact test
      # If there is a warning message, print it and change the condition
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
    }
  }
  # Format the p-value, using an HTML entity for the less-than sign.
  # The initial empty string places the output on the line below the variable label.
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

table <- table1(~ gender+
                  age_at_initial_pathologic_diagnosis.x+
                  BRCA_Subtype_PAM50+
                  ajcc_pathologic_tumor_stage+
                  histological_type+
                  new_tumor_event_type+
                  new_tumor_site_comb+
                  menopause_status+
                  pr_status_by_ihc+ 
                  her2_status_by_ihc
                |TRIM45_expression, data = merged,
                extra.col = list(`p-value` = pvalue), caption = caption,
                footnote = footnote, extra.col.pos = 3,
                render = rndr, render.categorical = "FREQ (PCTnoNA%)",
                topclass = " Rtable1-times ")

print(table)


unique(merged$histological_type)

merged <- merged %>% 
  mutate(across(where(is.factor), as.character)) %>%  
  mutate(across(where(is.character), ~ na_if(., "[Not Evaluated]")))

merged <- merged %>% 
  mutate(across(where(is.factor), as.character)) %>%  
  mutate(across(where(is.character), ~na_if(.,"[Unknown]")))

merged <- merged %>% 
  mutate(across(where(is.factor), as.character)) %>%  
  mutate(across(where(is.character), ~ na_if(., "[Not Available]")))



merged$ajcc_pathologic_tumor_stage <- factor(merged$ajcc_pathologic_tumor_stage,
                                             levels = c("Stage I","Stage IA", "Stage IB", "Stage II", "Stage IIA", "Stage IIB", "Stage III", "Stage IIIA", "Stage IIIB", "Stage IIIC", "Stage X", "[Discrepancy]" ),
                                             labels = c("Stage I", "Stage I", "Stage I", "Stage II", "Stage II", "Stage II", "Stage III", "Stage III", "Stage III", "Stage III", "Stage X", "Discrepancy"))

merged$histological_type <- factor(merged$histological_type,
                                       levels = c("Mixed Histology (please specify)","Other, specify", "Infiltrating Lobular Carcinoma", "Infiltrating Ductal Carcinoma", "Mucinous Carcinoma", "Metaplastic Carcinoma", "Medullary Carcinoma"),
                                       labels = c("Mixed Histology", "Other", "Infiltrating Lobular Carcinoma", "Infiltrating Ductal Carcinoma", "Mucinous Carcinoma", "Metaplastic Carcinoma", "Medullary Carcinoma"))

merged$OS <- factor(merged$OS ,
                    levels = c("0","1"),
                    labels = c("Alive", "Deceased"))

merged$DSS <- factor(merged$DSS ,
                    levels = c("0","1"),
                    labels = c("Died of other Causes", "Died of Disease"))


merged$DFI <- factor(merged$DFI ,
                    levels = c("0","1"),
                    labels = c("Disease Free", "Recurrence"))

merged$PFI <- factor(merged$PFI ,
                     levels = c("0","1"),
                     labels = c("No Prpgression", "Progression"))



#Make lookup table for the p-value function
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
  x <- x[-length(x)] # Remove "overall" group
  # Construct vectors of data y, and groups (strata) g
  y <- unlist(x)
  g <- factor(rep(1:length(x), times=sapply(x, length)))
  if (is.numeric(y)) {
    # check for normality using Shapiro-Wilk test and variance equality using var.test
    swtest <- shapiro.test(y)$p.value
    vartest <- var.test(y ~ g)$p.value
    if (all(swtest > 0.05, vartest > 0.05)) {
      # If both p-values are >0.05 we assume normality and use t.test with equal variances
      p <- t.test(y ~ g, var.equal=TRUE)$p.value
      symbol <- lookup_table[findInterval(p, cutoffs)+1]
    } else {
      # Otherwise, we assume normality with unequal variances and use t.test with var.equal=F
      p <- t.test(y ~ g, var.equal=FALSE)$p.value
      symbol <- lookup_table[findInterval(p, cutoffs)+1]
    } 
    if (all(swtest <= 0.05)) {
      # Otherwise, reject normality and use Wilcoxon test
      p <- wilcox.test(y ~ g, exact=FALSE)$p.value
      symbol <- lookup_table[findInterval(p, cutoffs)+1]
    }} 
  if (is.factor(y)) {
    # For categorical variables, perform a chi-squared test of independence
    p <- tryCatch({
      # Si warning chi-test, rC)aliser un Fisher
      chisq.test(table(y, g))$p.value
    }, warning = function(w) {
      # If a warning appears (due to expected cell counts <5 or too few observations), use Fisher's Exact test
      # If there is a warning message, print it and change the condition
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
    }
  } 
  # Format the p-value, using an HTML entity for the less-than sign.
  # The initial empty string places the output on the line below the variable label.
  c(paste0(sub("<", "&lt;", format.pval(p, digits=3, eps=0.001)), "<sup>", symbol, "</sup>"))
}


table <- table1(~ gender+
                  age_at_initial_pathologic_diagnosis.x+
                  ajcc_pathologic_tumor_stage+
                  histological_type+
                  new_tumor_event_type+
                  pr_status_by_ihc+ 
                  her2_status_by_ihc
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
                         dump_html = TRUE) {  # delete intermediate html files?
  
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

write_table1(table, "2025_03_25_clinical_table_tcga_trim45_high_vs_low.pdf")
