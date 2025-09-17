## GSE131418 dataset preprocessing
## Set groups
## PCA analysis

#########################################################################
## GEO data download and analysis
rm(list=ls())
library(BiocManager)
library(GEOquery)
options(stringsAsFactors = F)
library(limma)
library(GEOquery) # Load GEOquery R package

# Set working path to ".../raw data/Functional DEGs/GSE131418"
# Upload gene expression matrices
# Load from the downloaded file "series_matrix.txt"
eSet = getGEO(filename = "GSE131418_series_matrix.txt.gz", getGPL = F) # Load the locally downloaded GEO file

library(data.table)
# Load probe-level expression tables from the downloaded ".txt.gz" files
data <- read.table("GSE131418_Consortium_prim_met_GE_probe_level.txt.gz", header = TRUE, sep = "\t")
class(data)
data[1:10,1:10]

data2 = read.table("GSE131418_MCC_prim_met_GE_probe_level.txt.gz", header = T, sep = "\t")
# Combine two data.tables into one gene expression matrix
exp = cbind(data, data2)

# Normalize between arrays
boxplot(exp) # Inspect expression intensities across arrays
range(exp)
exp = normalizeBetweenArrays(exp)
range(exp)

# Determine whether to log-transform data
qx <- as.numeric(quantile(exp, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm = T))
LogC <- (qx[5] > 100) ||
  (qx[6] - qx[1] > 50 && qx[2] > 0)
if (LogC) {
  exp[which(exp <= 0)] <- NaN
  exp <- log2(exp + 1)
}

range(exp)

exp = as.data.frame(exp) # Convert matrix to data.frame

# Acquire sample metadata (clinical information)
ph = pData(eSet)
ph = as.data.frame(ph)

# ----------------------------------------------------------------------------- 

## Gene symbol mapping
# Upload list of gene symbols from file "GSE 131418 symbol.xlsx"
library(readxl)
ID = read_excel("GSE 131418 symbol.xlsx", sheet = 1, col_names = T) # Read gene symbol list
ID = as.data.frame(ID)

## Remove rows containing NA values in gene symbol column
ID <- ID[!is.na(ID$GeneSymbol), ]

## Map probe IDs to gene symbols
exp$ID = rownames(exp)
exp = merge(exp, ID, by = "ID")
exp = as.data.frame(exp)
exp[1:3,1135:1140]
exp = exp[,-c(1137,1138,1140)]
exp[1:5,1:5]
exp = exp[,-1]
rownames(exp) = exp$GeneSymbol

## Remove duplicate gene symbols in gene expression dataframe
which(duplicated(exp$GeneSymbol))
k = !duplicated(exp$GeneSymbol)
table(k)
exp = exp[k,]

## Set gene symbol as rownames of the dataframe
rownames(exp) = exp$GeneSymbol
colnames(exp)[1130:1136]
ll = head(exp)
exp = exp[,-1136]

##################
## Data filtering
# Filter out samples with prior medical treatment
table(ph$`treatment status classification (see description):ch1`)
ph = ph[ph$`treatment status classification (see description):ch1` %in% "PRE", ]

# Match untreated samples to gene expression dataframe
k = intersect(colnames(exp), ph$title)
ph = ph[ph$title %in% k, ]
exp = exp[, k]
identical(ph$title, colnames(exp))

## Ensure sample order consistency in gene expression dataframe
index = match(ph$title, colnames(exp))
exp = exp[, index]
identical(ph$title, colnames(exp))

## Save preprocessed data from GSE131418 as an RData file
save(exp, ID, ph, file = "GSE131418_Raw_Data.Rdata")

###################################################################

## Setting groups
# Sample labels: metastasis & primary colorectal cancer
# Summarize label distribution
table(ph$characteristics_ch1.6)
# Assign group labels
grouplist = ifelse(ph$characteristics_ch1.6 == "tumor type: METASTASIS", "metastasis", "primary")

# Convert group labels to factors
grouplist = factor(grouplist, levels = c("metastasis", "primary"))
table(grouplist)

# Save grouping labels
save(grouplist, file = "group_metastasis_primary.Rdata")

###################################################################
## PCA Analysis
rm(list=ls())

# Load sample data and group labels
# NOTE: ensure filename matches the saved grouping file above
load(file = "group_metastasis_primary.Rdata")
load(file = "GSE131418_Raw_Data.Rdata")

library(limma)
library(ggplot2)

## Data processing for PCA
dat = t(exp)
dat[is.na(dat)] = 0
df_pca = prcomp(dat)
df_pcs <- data.frame(df_pca$x, Species = grouplist)
head(df_pcs, 3)

## Visualization
plot(df_pca$x[,1], df_pca$x[,2])
ggplot(df_pcs, aes(x = PC1, y = PC2, color = Species)) + geom_point()

## PCA plot customization
ggplot(df_pcs, aes(x = PC1, y = PC2, color = Species)) +
  geom_point() +
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"))

## Add explained variance to PCA axis labels
percentage <- round(df_pca$sdev / sum(df_pca$sdev) * 100, 2)
percentage <- paste(colnames(df_pcs), "(", paste(as.character(percentage), "%", ")", sep = ""))
ggplot(df_pcs, aes(x = PC1, y = PC2, color = Species)) +
  geom_point() +
  xlab(percentage[1]) +
  ylab(percentage[2])

## Final PCA plot customization
ggplot(df_pcs, aes(x = PC1, y = PC2, color = Species)) +
  geom_point() +
  xlab(percentage[1]) +
  ylab(percentage[2]) +
  theme(panel.background = element_blank(),  # Remove background color
        panel.grid = element_blank())        # Remove grid lines

## Save plot as PDF in the working directory
pdf("PCA Analysis.pdf")
# Example plot (replace with actual plotting code if saving ggplot objects)
plot(1:10, 1:10, main = "Example plot")
# Close PDF device
dev.off()