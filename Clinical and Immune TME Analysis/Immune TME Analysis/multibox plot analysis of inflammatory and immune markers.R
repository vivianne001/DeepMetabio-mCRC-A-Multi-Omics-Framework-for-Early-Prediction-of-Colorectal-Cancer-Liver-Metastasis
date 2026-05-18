
###################### Multibox Plot Analysis of Inflammatory and Immune Markers in GSE204805 #########################
# GSE204805

rm(list = ls())

library(pheatmap)
library(RColorBrewer)
library(dplyr)
library(ggsignif)
library(ggplot2)
library(reshape)
library(reshape2)
library(tidyverse)
library(ggpubr)
library(IOBR)

# Load data
# load("GSE204805_processed.Rdata")

# Log2 transform expression data
exp = log2(exp + 1)

# Ensure sample order consistency using tpm gene list (tpm contains human gene symbols of interest)
kk = rownames(tpm) 
exp = exp[kk, ]

# Filter for LM cell type samples
ph = ph[ph$`cell type:ch1` == "LM", ]
k = intersect(rownames(ph), colnames(exp))
exp = exp[, k]
identical(rownames(ph), colnames(exp))

# Define marker gene sets
chemokine = c("CCR2", "CCR5", "CCR7", "CCL2", "CCL5", "CCL17", "CXCR4", "CXCL9", "CXCL10")
immune = c("ADORA2A", "BTLA", "CD160", "CD200", "CD200R1", "CD244", "CD27", "CD274", "CD28", "CD40", "CD40LG", "CD44", "CD48", "CD70", "CD80",
           "CD88", "ICOS", "ICOSLG", "IDO1", "IDO2", "KIR3DL1", "LAIR1", "LGALS9", "TMIGD2", "TNFRSF18", "TNFRSF25", "TNFRSF8", "TNFSF15")

# Prepare data for analysis: include immune markers plus ACMSD
dat = exp[c(immune, "ACMSD"), ]
dat = na.omit(dat)
dat = as.data.frame(t(dat))
dat$group = ifelse(dat$ACMSD > median(dat$ACMSD), "High", "Low")
dat$group = as.factor(dat$group)
dat <- subset(dat, select = -ACMSD)

# Example statistical tests for individual genes
t_test_result <- t.test(dat$CXCL11[dat$group == "High"], dat$CXCL11[dat$group == "Low"])
wilcox_result <- wilcox.test(dat$HAVCR2[dat$group == "High"], dat$HAVCR2[dat$group == "Low"])

# Function mapping p-values to significance symbols
map_pvalue_to_signif_mark <- function(pvalue) {
  if (pvalue < 0.001) {
    return("***")
  } else if (pvalue < 0.01) {
    return("**")
  } else if (pvalue < 0.05) {
    return("*")
  } else {
    return("ns")
  }
}

# Prepare long-format data for multibox plotting
dat$group <- factor(dat$group, levels = c("Low", "High"))
TME_NEW = melt(dat)
colnames(TME_NEW) = c("Group", "Gene", "Express")
TME_NEW$Express <- as.numeric(TME_NEW$Express)
TME_NEW$Group <- factor(TME_NEW$Group, levels = c("Low", "High"))

# Multibox plot for immune markers comparing ACMSD High vs Low groups
p <- ggplot(TME_NEW, aes(x = Gene, y = Express, fill = Group)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values = c("Low" = "lightblue", "High" = "red"), labels = c("Low", "High")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        axis.line = element_line(color = "black"),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "right") +
  stat_compare_means(aes(group = Group), method = "t.test", label = "p.signif", label.y = max(TME_NEW$Express)) +
  ylab("Expression (log2)") +
  xlab("Immune Marker")

print(p)

# Save marker gene list for downstream analyses
genes = c("ADORA2A", "BTLA", "CD160", "CD244", "CD28", "CD40", "CD40LG", "IDO1", "IDO2", "LAIR1", "TNFRSF25", "TNFRSF1A", "TNFRSF1B", "CCR5",
          "CCR6", "CXCL9", "CXCL10", "CXCR3", "IL7", "NFKB1", "TLR2", "TLR7")
save(genes, file = "immune_genes_for_boxplot_analysis.rdata")

# Additional boxplots for a selected immune marker subset
imm = c("BTLA", "CD160", "CD200R1", "CD244", "CD28", "CD40", "CD40LG", "CD48", "CD80", "ICOS", "IDO1", "IDO2", "TNFRSF18", "TNFRSF8")
dat = exp[c(imm, "ACMSD"), ]
dat = na.omit(dat)
dat = as.data.frame(t(dat))
dat$group = ifelse(dat$ACMSD > median(dat$ACMSD), "High", "Low")
dat <- subset(dat, select = -ACMSD)
TME_NEW = melt(dat)
colnames(TME_NEW) = c("Group", "Gene", "Express")
TME_NEW$Express <- as.numeric(TME_NEW$Express)
TME_NEW$Group <- factor(TME_NEW$Group, levels = c("Low", "High"))

p <- ggplot(TME_NEW, aes(x = Gene, y = Express, fill = Group)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values = c("Low" = "lightblue", "High" = "red"), labels = c("Low", "High")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        axis.line = element_line(color = "black"),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "right") +
  stat_compare_means(aes(group = Group), method = "t.test", label = "p.signif", label.y = max(TME_NEW$Express)) +
  ylab("Expression (log2)") +
  xlab("Immune Marker")

print(p)