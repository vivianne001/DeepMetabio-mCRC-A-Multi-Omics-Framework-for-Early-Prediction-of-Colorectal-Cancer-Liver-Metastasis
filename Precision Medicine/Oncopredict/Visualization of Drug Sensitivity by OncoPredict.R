###################### Visualization of Drug Sensitivity Differences Predicted by OncoPredict #########################

# Load OncoPredict results
result = load("Oncopredict-GSE204805 drugs.Rdata")
testPtype = result
rownames(testPtype) = result$X
testPtype = testPtype[, -1]

# Prepare group information
k = rownames(testPtype)
new_data = new_data[new_data$id %in% k, ]
testPtype$group = new_data$new_column
new_data$group = new_data$new_column

# Example for merging different drug response matrices (if needed)
# testPtype1 = testPtype_GDSC1
# testPtype2 = testPtype_GDSC2
# testPtype = merge(testPtype1, testPtype2, by = "V1")
# testPtype = testPtype_CTRP2

# Extract drug response columns
# Pick specific drugs of interest(Compounds targeting VEGFR and EGFR TKIs)
dabrafenib = testPtype$Dabrafenib_1373.x
axitinib = testPtype$Axitinib_1021
brivanib = testPtype$`Brivanib, BMS-540215_376`
cabozantinib = testPtype$Cabozantinib_249
cetuximab = testPtype$Cetuximab_1114
linifanib = testPtype$Linifanib_277
ponatinib = testPtype$Ponatinib_155
sunitinib = testPtype$Sunitinib_5
tivozanib = testPtype$Tivozanib_312

rowname = rownames(testPtype)
oncopredict = data.frame(
  "dabrafenib" = dabrafenib,
  "axitinib" = axitinib,
  "brivanib" = brivanib,
  "cabozantinib" = cabozantinib,
  "cetuximab" = cetuximab,
  "linifanib" = linifanib,
  "ponatinib" = ponatinib,
  "sunitinib" = sunitinib,
  "tivozanib" = tivozanib
)
rownames(oncopredict) = rowname
oncopredict$group = ph$group

write.csv(oncopredict, file = "oncopredict_results.csv")

# Match group information
table(ph$group)
k = intersect(testPtype$V1, rownames(ph))
ph = ph[k, ]
testPtype = testPtype[testPtype$V1 %in% k, ]
testPtype$group = ph$group
rownames(testPtype) = testPtype$V1
testPtype = testPtype[, -1]

# Remove NA values for a specific drug
testPtype.pr <- testPtype[!is.na(testPtype$BRD_K30748066), ]

# Statistical tests
t_test_result <- t.test(testPtype$Tivozanib_312[testPtype$group == "high"], 
                        testPtype$Tivozanib_312[testPtype$group == "low"])
wilcox_result <- wilcox.test(testPtype$dabrafenib[testPtype$group == "high"], 
                             testPtype$dabrafenib[testPtype$group == "low"])

library(ggsignif)
library(ggplot2)

map_pvalue_to_signif_mark <- function(pvalue) {
  if (pvalue < 0.001) {
    return("***")
  } else if (pvalue < 0.01) {
    return("**")
  } else if (pvalue < 0.055) {
    return("*")
  } else {
    return("ns")
  }
}

significance_symbols <- map_pvalue_to_signif_mark(t_test_result$p.value)
significance_symbols <- map_pvalue_to_signif_mark(wilcox_result$p.value)

# Identify outliers for a specific drug
outliers <- boxplot.stats(testPtype.pr$BRD_K30748066)$out

# Single drug boxplot (example: Dabrafenib_1373)
ggplot(testPtype, aes(x = group, y = Dabrafenib_1373, fill = group)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values = c("high" = "red", "low" = "blue")) +
  geom_signif(comparisons = list(c("high", "low")), 
              test = "t.test", 
              map_signif_level = TRUE) +
  geom_text(aes(x = 1.5, y = max(Dabrafenib_1373), label = paste("p-value =", format(t_test_result$p.value, digits = 3, scientific = TRUE))), 
            hjust = 0, vjust = 1) +
  theme_minimal() +
  ylab("BRD_K30748066 IC50(μM)")



# Multi-drug boxplot
options(stringsAsFactors = FALSE)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(ggsignif)
library(reshape)
library(IOBR)

TME_NEW = melt(testPtype)
colnames(TME_NEW) = c("Group", "drug", "IC_50")

# Order drugs by median IC50 in high group
plot_order = TME_NEW[TME_NEW$Group == "high", ] %>%
  group_by(drug) %>%
  summarise(m = median(IC_50)) %>%
  arrange(desc(m)) %>%
  pull(drug)

TME_NEW$drug <- gsub("_\\d+", "", TME_NEW$drug)
colnames(TME_NEW) = c("Group", "drug", "IC_50")

# Ensure Group is a factor with correct order
TME_NEW$Group <- factor(TME_NEW$Group, levels = c("low", "high"))

# Multi-drug boxplot with significance
p <- ggplot(TME_NEW, aes(x = drug, y = IC_50, fill = Group)) +
  geom_boxplot(width = 0.5) +
  scale_fill_manual(values = c("lightblue", "red")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(colour = "black")) +
  stat_compare_means(aes(group = Group), method = "t.test", label = "p.signif", label.y = 1.1