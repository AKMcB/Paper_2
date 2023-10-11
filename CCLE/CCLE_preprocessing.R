library(tibble)
library(dplyr)

expr1 <- read.csv2("Expression_Public_22Q4.csv", sep = ",", as.is = T, check.names = F)

info <- read.csv2("Model.csv", sep = ",", as.is = T, check.names = F)
info <- info[,c(1,3,11,13,17,21:24)]

rownames(expr1) <- expr1[,1]
expr1 <- expr1[,-1]
expr <- as.data.frame(t(expr1))
expr <- tibble::rownames_to_column(expr,"ID")
dup<-expr[duplicated(expr$ID),]

rownames(expr1) <- expr1[,1]
expr1 <- expr1[,-1]
expr <- as.data.frame(t(expr1))

merged <- merge(expr1, info, by.x = "ID", by.y = "ModelID")
merged <- merged %>% relocate(c(19146:19153), .after = 1)

merged <- merged[!merged$OncotreePrimaryDisease == "Non-Cancerous", ]

#Subset on Primary Disease
#code <- merged[,8]
#write.csv2(code, "Cancer_Code_CCLE.csv")

code <- read.csv2("Cancer_Code_CCLE.csv", sep = ";", as.is = T, check.names = F)


                                     merged_1 <- cbind(merged, code)
merged_1 <- merged_1 %>% relocate(c(19154:19155), .after = 1)
merged_1 <- merged_1[,-c(2,5:8,10)]
write.csv2(merged_1, "CCLE_exp_info_csv")

#Minimize the df
merged_2 <- merged_1 %>% relocate("TRIM5", .after = 1)
merged_2 <- merged_2[,-c(7:19149)]
T5 <- merged_2
rm(merged_2, dup, code, expr, expr1, merged, merged_1, TRIM5)

write.csv2(T5, "CCLE_TRIM5_exp_info_csv")




