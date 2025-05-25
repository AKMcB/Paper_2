##################
# Load libraries #
##################

library(clusterProfiler)
library(org.Hs.eg.db)
library(data.table)
library(dplyr)
library(ggpubr)

##################
# Load Gene List #
##################
df <- read.csv2("2025_03_25_limma_trim45_top_bottom_tcga_brca_er_positive_bh.csv")

df <- df %>% mutate_at(vars(logFC, AveExpr, t, P.Value, adj.P.Val, B), as.numeric)

#Upregulated/downregulated genes 
df <- df[df$logFC < 0,]

#only significant results
df <- df[df$adj.P.Val <= 0.05,]


#Check for duplicates
df <- distinct(df)


#Convert the gene symbols to ENTREZ ID
df <- clusterProfiler::bitr(geneID = df$hgnc_symbol,fromType="SYMBOL", toType="ENTREZID",
                            OrgDb="org.Hs.eg.db")

#Perform the GO enrichment
set.seed(123)
unique_GO <- enrichGO(gene = df$ENTREZID, #Your genes of interest
                      keyType = "ENTREZID",
                      ont = "BP", #Biological pathway
                      OrgDb = org.Hs.eg.db, #If no background is used, the function will use this instead
                      readable = T,
                      pvalueCutoff = 0.05, #It should be 0,05
                      qvalueCutoff = 0.25) #It should be 0,25.


#result <- unique_GO@result
#write.csv2(result, "upregulated_genes_go_term_high_vs_low_t45_results.csv")

p <- barplot(unique_GO, 
             drop = FALSE, 
             showCategory = 11, 
             title = "GO Enrichment for downregulated genes in high vs low TRIM45:BP",font.size = 15)
p

#Save the files 
pdf("GO_enrichment_metabric_for_downregulated_genes_high_vs_low_trim45_MF.pdf", width = 8, height = 10, onefile = F)
print(p)
dev.off()


png("GO_enrichment_for_downregulated_genes_high_vs_low_trim45_CC.png", res= 200, width = 1800, height = 2200)
print(p)
dev.off()

