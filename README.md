# Script for Paper II

#### The folders are divided based on what data has been used: 
- **CCLE**: Scripts for the data from the cancer cell line encyclopedia 
- **PANCAN**: Scripts for the data from the Pan Cancer Atlas project 
- **TCGA-BRCA**: Scripts for the data from the The Cancer Genome Atlas program 
- **METABRIC**: Scripts for the data from METABRIC
- **Gene regulatory network**: Scripts for handling edge weights

##### CCLE 
- **CCLE_preprocessing.R**
    - Preprocessing of CCLE data.
- **CCLE_boxplots_all_lines.R**
    - This script explore the CCLE data and the expression of a specific gene, and all cell lines are included.
- **Breast cancer cell lines_script_V4_boxplot.R**
    - This scripts explores breast cancer cell lines from CCLE and the expression of a specific gene.

##### PANCAN 
- **PanCan_Exp_boxplot.R** 
    - This script explore the PanCan data and the expression of a specific gene, and include all the 33 cancer types. 

##### TCGA-BRCA 
- **Survival_script_TCGA_ER_Status_V3.R** 
   - This script explores the different survival parameters such as disease-specific survival and progression free survival in TCGA-BRCA, based on the expression values of specific genes. It also includes calculations for GMM. 
- **survival_subtype.R**
    -This script explores the different survival parameters such as disease-specific survival and progression free survival in TCGA-BRCA PAM50 subtypes, based on the expression values of specific genes. 
- **Pearson_TCGA_ER_status.R** 
    - This script caluclates the pearson correlation between two genes in ER-positive TCGA-BRCA patients.
- **BRCA_normal_vs_tumor_exp.R**
    - This scripts explores TRIM45 expression in normal versus tumor tissue. 
- **making_boxplots_TCGA_V3.R**
    - Explores the expression of a specific gene in the TCGA-BRCA PAM50 subtypes. 
- **clinical_table_tcga.R**
    -Create a clinical table based on TRIM45 expression and clinical features.  
- **limma_brca.R**
    - LIMMA analysis of TCGA-BRCA based on TRIM45 expression.
- **go_analysis_limma.R**
    - Performs GO enrichemnt of LIMMA results. 

#### METABRIC
- **Survival_script_Metabric_ER_Status.R** 
   - This script explores the different survival parameters such as disease-specific survival and progression free survival in METABRIC, based on the expression values of specific genes. It also includes calculations for GMM. 
- **survival_subtype.R**
    -This script explores the different survival parameters such as disease-specific survival and progression free survival in METABRIC PAM50 subtypes, based on the expression values of specific genes. 
- **Pearson_Metabric_ER_status.R** 
    - This script caluclates the pearson correlation between two genes in ER-positive METABRIC patients.
- **making_boxplots_metabric_V2.R**
    - Explores the expression of a specific gene in the METABRIC PAM50 subtypes. 
- **clinical_table_metabric.R**
    -Create a clinical table based on TRIM45 expression and clinical features. 
- **limma_metabric.R**
    - LIMMA analysis of METABRIC based on TRIM45 expression.
- **go_analysis_metabric.R**
    - Performs GO enrichemnt of LIMMA results. 

##### Gene regulatory network analysis
- **edge_weights_TF_TRIM45.R**
    - Explores the edge weights between ESR1 and TRIM45 in TCGA-BRCA ER-positive patients. 
- **edge_weights_TF_TRIM45.R**
    - Explores the edge weights between ER-response transcription factors and TRIM45 in TCGA-BRCA ER-positive patients. 

- **venndiagram_tcga-meta.R**
    - Creates a venndiagram comparing genes between TCGA-BRCA and METABRIC LIMMA analysis. 
- **boxplot_GSE27473.R**
    - Process and creates a boxplots based on the data from GSE27473. 
- **boxplot_GSE74032.R**
    - Process and creates a boxplots based on the data from GSE74032. 
