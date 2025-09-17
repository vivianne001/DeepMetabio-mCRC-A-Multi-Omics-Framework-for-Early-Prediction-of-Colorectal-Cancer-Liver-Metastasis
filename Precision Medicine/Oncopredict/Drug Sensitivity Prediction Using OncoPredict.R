###################### Drug Sensitivity Prediction Using OncoPredict in GSE204805 #########################

rm(list = ls())

options(stringsAsFactors = FALSE)
library(oncoPredict)
library(data.table)
library(gtools)
library(reshape2)
library(ggpubr)
library(stringr)

# Load OncoPredict drug sensitivity dataset
drug = read.csv("drug1 info of oncopredict training set.csv", header = TRUE)
drug = drug[, -1]

# Clean target and pathway columns
drug$target <- iconv(drug$target, to = "UTF-8")
drug$target <- gsub("[^a-zA-Z0-9 ]", ",", drug$target)
drug$target <- gsub(",+", ",", drug$target)
drug$target <- gsub("^,|,$", "", drug$target)

drug$pathway <- iconv(drug$pathway, to = "UTF-8")
drug$pathway <- gsub("[^a-zA-Z0-9 ]", ",", drug$pathway)
drug$pathway <- gsub(",+", ",", drug$pathway)
drug$pathway <- gsub("^,|,$", "", drug$pathway)

# Example: filter drugs by target or pathway
drug_target1 = drug[grepl("Cetuximab", drug$drug), ]

# other drug_targets include dabrafenib, axitinib,brivanib,cabozantinib,linifanib,ponatinib,sunitinib,tivozanib

# Load training datasets
exp_CTRP2 = log2(`CTRP2_Expr (RPKM, not log transformed)` + 1)
exp_GDSC1 = `GDSC1_Expr (RMA Normalized and Log Transformed)`
exp_GDSC2 = `GDSC2_Expr (RMA Normalized and Log Transformed)`

# Standardize column names for exp_CTRP2
new_colnames <- str_extract(colnames(exp_CTRP2), "CVCL_[^)]+")
colnames(exp_CTRP2) <- new_colnames
exp_CTRP2 = as.matrix(exp_CTRP2)

filename <- c("`CTRP2_Expr (TPM, not log transformed)`", "`GDSC1_Expr (RMA Normalized and Log Transformed)`", "`GDSC2_Expr (RMA Normalized and Log Transformed)`")
filename_no_backticks <- gsub("`", "", filename)
rm(list = c(filename_no_backticks))

save(CTRP2_Res, GDSC1_Res, GDSC2_Res, exp_CTRP2, exp_GDSC1, exp_GDSC2, file = "Training Datasets.Rdata")

# Load testing data
load("ACMSD group_high and low.Rdata")
exp = exp[c(1:43933), ]
exp_filtered <- exp[rowSums(exp == 0) < ncol(exp), ]
index = rownames(exp_filtered)
exp = exp[index, ]

table(clin$SUBTYPE)
class(exp)
range(exp)

# Prepare training and testing datasets
k = new_data$id
exp = exp[, k]
exp = as.matrix(exp)

# Drug category information
drug_category = colnames(GDSC1_Res)
drug_category = as.data.frame(drug_category)
drug_category$group = c("GDSC1")
colnames(drug_category)[1] = "drug_name"

drug_category1 = colnames(GDSC2_Res)
drug_category1 = as.data.frame(drug_category1)
drug_category1$group = c("GDSC2")
colnames(drug_category1)[1] = "drug_name"

drug_category2 = colnames(CTRP2_Res)
drug_category2 = as.data.frame(drug_category2)
drug_category2$group = c("CTRP2")
colnames(drug_category2)[1] = "drug_name"

drug_category = rbind(drug_category, drug_category1)
drug_category = rbind(drug_category, drug_category2)

save(drug, drug_category, file = "oncopredict-drug and drug category.Rdata")

CDK = drug_category[drug_category$drug_name %in% index, ]
folate = drug_category[drug_category$drug_name %in% index, ]
RAS = drug_category[drug_category$drug_name %in% ras, ]

# Example: select specific drugs from response matrices
column_index = which(colnames(GDSC1_Res) == "Navitoclax_1011")
GDSC1_drug = GDSC1_Res[, c(313, 228)]
column_index = which(colnames(GDSC2_Res) == "Navitoclax_1011")
GDSC2_drug = GDSC2_Res[, c(66, 48, 7)]
column_index = which(colnames(CTRP2_Res) == "Navitoclax_1011")
CTRP2_drug = CTRP2_Res[, c(531, 2)]

# Batch query for RAS drugs
drug_names <- RAS$drug_name[3:9]
column_indices <- match(drug_names, colnames(CTRP2_Res))
CTRP2_drug <- CTRP2_Res[, column_indices]

##################### OncoPredict Prediction ###########################
Sys.setenv(R_MAX_PPSTACK = '500000')
exp = as.data.frame(exp)

calcPhenotype(
  trainingExprData = exp_GDSC2,
  trainingPtype = GDSC2_drug,
  testExprData = exp,
  batchCorrect = 'eb',
  powerTransformPhenotype = TRUE,
  removeLowVaryingGenes = 0.2,
  minNumSamples = 10,
  printOutput = TRUE,
  removeLowVaringGenesFrom = 'rawData'
)

testPtype <- fread('./calcPhenotype_Output/DrugPredictions.csv', data.table = FALSE)
range(testPtype$Sorafenib_1085)
testPtype[is.na(testPtype)] = 0

testPtype_CTRP2 = testPtype
testPtype_GDSC1 = testPtype
testPtype_GDSC2 = testPtype

save(testPtype_CTRP2, testPtype_GDSC1, testPtype_GDSC2, file = "Oncopredict-GSE204805 drugs.Rdata")
testPtype = testPtype_CTRP2
##################### Correlation Analysis and Heatmap Visualization ###########################