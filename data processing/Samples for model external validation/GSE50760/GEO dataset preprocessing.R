rm(list = ls())
#########################################################################
## GEO data download and storage
library(BiocManager)
library(GEOquery)

# Simple code to download GSE data
options(stringsAsFactors = F)
library(limma)
library(GEOquery) # Load GEOquery R package
gse = "GSE39582" # Select GSE series to download (example)
eSet = getGEO(gse,
              destdir = '.',
              getGPL = F) # Download GSE; eSet is a list containing expression matrices

# If the above method fails, download the series_matrix.txt file from GEO manually
eSet = getGEO(filename = "GSE50760_series_matrix.txt.gz", getGPL = F) # Load locally downloaded GEO file

########################
## Read expression matrix for GSE50760
library(readxl)
exp = read_excel("GSE50760-EXP.xlsx", col_names = T)

exp = exp[!duplicated(exp$genes), ]
exp = as.data.frame(exp)
rownames(exp) = exp$genes
exp = exp[, -1]

range(exp)

k = intersect(colnames(exp), rownames(ph))
exp = exp[, k]
ph = ph[k, ]

identical(colnames(exp), rownames(ph))

#############
## Import clinical metadata for GSE50760

# Note: some extreme (outlier) values were observed in `exp`
# which(exp > 500, arr.ind = T) # Identify outliers and their positions
# exp = exp[, -c(42,107)] # Remove patient samples with extreme values (if needed)
# boxplot(exp)
# exp[is.na(exp)] = 0 # Replace NA values in exp with 0 (optional)

# range(exp) # Use range() to inspect min/max and assess if data is log-transformed

# Obtain clinical metadata
ph = pData(eSet) # Retrieve phenotype/clinical information
ph = as.data.frame(ph)

# ----------------------------- Subset data to CRLM and primary CRC samples #################
table(ph$source_name_ch1)
ph = ph[ph$source_name_ch1 %in% c("metastasized cancer", "primary colorectal cancer"), ]

k = intersect(colnames(exp), rownames(ph))
exp = exp[, k]
ph = ph[k, ]

identical(colnames(exp), rownames(ph))

save(exp, ph, file = "GSE50760_Raw_Data.Rdata")

###################################################################

# 4. Group assignment
# 1) Group by presence of liver metastasis (metastasis / primary)
table(ph$source_name_ch1)

grouplist = ifelse(ph$source_name_ch1 == "metastasized cancer", "metastasis", "primary")

# Set reference level so control group appears first and treated/target group second
grouplist = factor(grouplist, levels = c("metastasis", "primary"))
table(grouplist)

# Save group labels
save(grouplist, file = "group_metastasis_primary.Rdata")