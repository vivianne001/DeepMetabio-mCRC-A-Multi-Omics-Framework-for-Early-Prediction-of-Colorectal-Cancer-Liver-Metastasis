rm(list = ls())
#########################################################################
## GSE204805 dataset preprocessing
## set groups

## 1. GEO data download and analysis
rm(list=ls())
library(BiocManager)
library(GEOquery)
options(stringsAsFactors = F)
library(limma)
library(GEOquery) # load GEOquery R Package

# set working path to ".../raw data/Model Traning/GSE204805" if needed
# upload gene expression matrix from the downloaded "series_matrix.txt"
eSet=getGEO(filename = "GSE204805_series_matrix.txt.gz", getGPL = F) # open local GEO file

## 2. Read expression data (tsv.gz) and preprocess gene IDs
library(readr)
data = read_tsv("GSE204805_merged_hs_mm.tsv.gz")
data[1:10,1:10]
class(data)

data = as.data.frame(data)
colnames(data)

# remove "H_" prefix
data$Geneid = sub("H_","",data$Geneid)
# remove the suffix after "."
data$Geneid = sapply(strsplit(data$Geneid, "\\."), "[", 1)

# remove duplicated Geneid
data = data[!duplicated(data$Geneid),]

# set expression matrix
exp = data
range(exp) # check the range of exp

exp = as.data.frame(exp)
rownames(exp) = exp$Geneid
exp = exp[,-1]

# diagnostics for extreme values (optional)
# which(exp > 500, arr.ind = TRUE) # locate extreme values
# exp = exp[,-c(42,107)]           # remove samples with extreme values
# boxplot(exp)
# exp[is.na(exp)] = 0              # replace NA with 0 if needed
# range(exp)                       # check whether exp is log-transformed

## 3. Acquire clinical metadata
ph = pData(eSet) # clinical information
ph = as.data.frame(ph)
# split by comma and keep the first token
ph$title = sapply(strsplit(ph$title, ","), "[", 1)
ph$title[1:10]

rownames(ph) = ph$title
ph = ph[,-1]
# -----------------------------------------------------------------------------

## 4. Grouping info and annotation (probe ID -> Gene Symbol)
# NOTE: The following line clears the environment; keep only if intended.
rm(list=ls())
load(file="GSE41258_RAW DATA.Rdata") # ensure this file exists and contains required objects
exp = as.data.frame(exp)

# annotation: copy from GPL page to an Excel file
library(readxl)
ID = read_excel("GSE159216 ID.xlsx", sheet = 2, col_names = TRUE)
ID = as.data.frame(ID)
ID = ID[-c(1:14),]
colnames(ID) = ID[1,]
ID = ID[-1,]

ID = ID[,c(1,8)]
ID = as.data.frame(ID)
ID[1:3,1:2]

# parse gene_assignment column
split_strings <- strsplit(ID$gene_assignment, " // ")
second_words <- sapply(split_strings, function(x) x[2])
ID$gene_assignment <- second_words

# remove rows starting with "LOC" or equal to "---"
bad_rows <- grepl("^LOC", ID$gene_assignment) | ID$gene_assignment == "---"
ID <- ID[!bad_rows, ]

# drop NA in gene_assignment
ID <- ID[!is.na(ID$gene_assignment), ]

# map probe ID to Gene Symbol
exp = as.data.frame(exp)
exp$ID = rownames(exp)
exp = merge(exp, ID, by = "ID")
exp = as.data.frame(exp)
exp[1:3,280:285]

# strip " /// .*" from Gene Symbol
exp$`Gene Symbol` <- sub(" /// .*", "", exp$`Gene Symbol`)

# remove duplicated gene symbols (by gene_assignment)
duplicates <- duplicated(exp$gene_assignment)
exp <- exp[!duplicates, ]
rownames(exp) = exp$gene_assignment

# drop helper columns
exp = exp[,-285]
exp = exp[,-1]

## 5. Sample filtering and alignment
# filter treated and non-human samples
ph = ph[ph$`treatment:ch1` %in% "NA",]
ph = ph[ph$characteristics_ch1.1 %in% c("cell line: H"),]

# align samples between exp and ph
k = intersect(colnames(exp), rownames(ph))
ph = ph[k,]
exp = exp[,k]
identical(colnames(exp), rownames(ph))
print(setequal(colnames(exp), ph$title))

## 6. Gene name sanitization (valid identifiers and de-duplication)
# find invalid rownames (not matching make.names)
invalid_colnames <- rownames(exp)[!make.names(rownames(exp), unique=TRUE) %in% rownames(exp)]
print(invalid_colnames)

# generate valid, unique rownames
new_rownames <- make.names(rownames(exp), unique=TRUE)
rownames(exp) <- new_rownames

# standardize (remove version suffix like ".1")
rownames(exp) <- gsub("^(.*?)\\..*$", "\\1", rownames(exp))

# also keep a clean gene name column
exp$gene_name = rownames(exp)
exp$gene_name <- gsub("^(.*?)\\..*$", "\\1", exp$gene_name)
exp$gene_name[1:100]

# remove duplicated genes and reset rownames
exp = exp[!duplicated(exp$gene_name),]
rownames(exp) = exp$gene_name
exp = exp[,-964]

## 7. Save data and set groups
save(exp, ph, file = "GSE204805_Raw_Data.Rdata")

# define metastasis vs primary groups
grouplist = ifelse(ph$`cell type:ch1` == "LM", "metastasis", "primary")
# set factor levels (reference/control first)
grouplist = factor(grouplist, levels = c("metastasis","primary"))
table(grouplist)

################################
# build a sample-level data frame
dat = t(exp)
dat = as.data.frame(dat)
dat$status = grouplist

# save
save(grouplist, file = "group_metastasis_primary.Rdata")
