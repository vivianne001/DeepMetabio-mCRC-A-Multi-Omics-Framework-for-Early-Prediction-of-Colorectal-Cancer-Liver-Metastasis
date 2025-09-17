################### Merge datasets and define sample groups ###################
rm(list = ls())
################### Merge datasets and expression matrices ##################
## Merge RMA datasets (GSE81986, GSE41568, GSE71222, GSE51244, GSE21510, GSE18105)
## Import RMA datasets

# load("GSE81986_Adjusted Data_Gene Symbol.Rdata")
exp1 = exp
exp1$ID = rownames(exp1)
ph1 = ph
ph1$group = grouplist
ph1 = ph1[, c(2,35)]
ph1$dataset = "GSE81986"

##
# load("GSE27854_Raw_Data.Rdata")
exp2 = exp
exp2$ID = rownames(exp2)
ph2 = ph
ph2$group = grouplist
ph2 = ph2[, c(2,40)]
ph2$dataset = "GSE27854"

##
# load("GSE71222_Raw_Data.Rdata")
exp3 = exp
exp3$ID = rownames(exp3)
ph3 = ph
ph3$group = grouplist
ph3 = ph3[, c(2,39)]
ph3$dataset = "GSE71222"

##
# load("GSE41568_Raw_Data.Rdata")
exp4 = exp
exp4$ID = rownames(exp4)
ph4 = ph
ph4$group = grouplist
ph4 = ph4[, c(2,35)]
ph4$dataset = "GSE41568"

##
# load("GSE21510_Raw_Data.Rdata")
exp5 = exp
exp5$ID = rownames(exp5)
ph5 = ph
ph5$group = grouplist
ph5 = ph5[, c(2,41)]
ph5$dataset = "GSE21510"

##
# load(GSE18105_Raw_Data.Rdata")
exp6 = exp
exp6$ID = rownames(exp6)
ph6 = ph
ph6$group = grouplist
ph6 = ph6[, c(2,40)]
ph6$dataset = "GSE18105"

## Concatenate datasets
## `exp` expression matrix
exp = merge(exp1, exp2, by = "ID")
exp = merge(exp, exp3, by = "ID")
exp = merge(exp, exp4, by = "ID")
exp = merge(exp, exp5, by = "ID")
exp = merge(exp, exp6, by = "ID")
colnames(exp)[1:10]
exp[1:10,1:10]
rownames(exp) = exp$ID
exp = exp[,-1]

## `ph` clinical metadata
ph = rbind(ph1, ph2)
ph = rbind(ph, ph3)
ph = rbind(ph, ph4)
ph = rbind(ph, ph5)
ph = rbind(ph, ph6)

##
save(exp, ph, file = "GEO_mCRC_datasets.rdata")

###############
