
rm(list = ls())
#########################################################################
## GSE41568 dataset preprocessing
## set groups

## GEO data download and analysis
rm(list=ls())
library(BiocManager)
library(GEOquery)
options(stringsAsFactors = F)
library(limma)
library(GEOquery) # load GEOquery R Package

# set working path to "...\raw data\Model Traning\GSE41568"
# upload gene expression matrix
# upload from the download file "series_matrix.txt"
eSet=getGEO(filename = "GSE41568_series_matrix.txt.gz",getGPL = F) # read in GEO Format

##
# get gene expression matrix from eSet 
exp=exprs(eSet) 
LL=head(exp) 
str(exp) 
exp=as.matrix(exp)

# normalized between arrays
boxplot(exp) # check the intensity of expression between arrays
range(exp)
exp=normalizeBetweenArrays(exp) 
range(exp)

# access the need to normalize expression matrix
# whether to log data or not
qx <- as.numeric(quantile(exp, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { exp[which(exp <= 0)] <- NaN
exp <- log2(exp+1) }

range(exp) 

exp=as.data.frame(exp) # transform matrix to dataframe

# acquire sample clinical info matrix 
ph=pData(eSet)#获得临床信息
ph=as.data.frame(ph)

# save gene expression matrix and samples clinical infos as rdata
save(ph,exp,file = "GSE41568_raw data.Rdata")
# -----------------------------------------------------------------------------

## Gene Symbol labels match to gene IDs 
# upload list of gene symbols from website "GPL570.xlsx"

library(readxl)
ID=read_excel("GPL570.xlsx",sheet = 1,col_names = T) # 
ID=as.data.frame(ID)

# use na.omit function to diminish rows containing NA
ID <- ID[!is.na(ID$`Gene Symbol`), ]

# Transfer IDs to Gene Symbol
exp$ID=rownames(exp)
exp=merge(exp,ID,by="ID")
exp = as.data.frame(exp)
exp[1:3,130:137]
exp = exp[,-c(135,137)]
exp[1:5,1:5]
exp = exp[,-1]
rownames(exp) = exp$`Gene Symbol`

# duplicate gene symbols in gene expression dataframe
which(duplicated(exp$`Gene Symbol`))
k = !duplicated(exp$`Gene Symbol`);
table(k)
exp=exp[k,]

# copy gene symbol to rownames of the dataframe
rownames(exp)=exp$`Gene Symbol`
colnames(exp)[130:134]
ll=head(exp)
exp=exp[,-134]

# match wanted samples to gene expression dataframe
k=intersect(colnames(exp),rownames(ph))
ph=ph[k,]
exp=exp[,k]
identical(rownames(ph),colnames(exp))

####
## data processing
# delete contents after "///" in rownames
new_rownames = sub(" /// .*", "", rownames(exp))
# delete repeated rownames
new_rownames = make.unique(new_rownames)

# copy new_rownames to rownames of exp dataframe
rownames(exp) = new_rownames
rownames(exp)[1:10]

## save preprocessed data from GSE41568 as rdata
save(exp,ID,ph,file = "GSE41568_Raw_Data.Rdata")

###################################################################

## setting groups
# labels : metastasis & Primary CRC
# summary of labels 
table(ph$`metastatic tumor site:ch1`)
# setting group of samples to metastasis & primary
grouplist=ifelse(ph$`metastatic tumor site:ch1` =="NA","primary","metastasis")

# setting grouping labels as factors
grouplist=factor(grouplist,levels = c("metastasis","primary"))
table(grouplist)

# save grouping samples 
save(grouplist,file = "group_metastasis_primary.Rdata")

