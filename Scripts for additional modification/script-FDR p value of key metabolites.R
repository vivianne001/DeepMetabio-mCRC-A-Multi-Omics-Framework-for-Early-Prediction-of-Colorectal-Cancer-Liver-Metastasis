
rm(list=ls())
######################### Quantitative Analysis of Metabolite sequencing data #################
## 1. Load Data
library(readxl)

# Import dataset - Metabolite LC/MS sequencing matrix
exp = read_excel("exp.xlsx", col_names = T, sheet = 1)
exp = as.data.frame(exp)
exp = exp[,-2]
rownames(exp) = exp$Compounds
exp = exp[,-1]

ph = read_excel("samples.xlsx", col_names = T, sheet = 1)
ph = as.data.frame(ph)

# View dataset
class(exp$`CRC-1`)

# Convert all columns to numeric
exp[] <- lapply(exp, as.numeric)
class(exp)
ph$Group = as.factor(ph$Group)
############
dat = as.data.frame(t(exp))
dat$Group = ph$Group

######################## Visualization ####################################

##############################################################
#########################
## 2. Boxplot
# Set the grouping column "Group" for the two groups
# Comparison of median values between the two groups
# Check p-value between the two groups for Erythrose
# Check p-value of IHC between CRLM & CRC AIS groups

# 1) Calculate t.test
t_test_result <- t.test(dat$`Gluconic acid`[dat$Group == "CRLM"], 
                        dat$`Gluconic acid`[dat$Group == "CRC"])

# 2) Calculate wilcox.test
wilcox_result <- wilcox.test(dat$`L-Proline`[dat$Group == "CRLM"], 
                             dat$`L-Proline`[dat$Group == "CRC"])

# Load ggsignif package
library(ggsignif)
library(ggplot2)

# Define a custom function to map p-values to significance marks
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

# Get significance marks
significance_symbols <- map_pvalue_to_signif_mark(wilcox_result$p.value)

# Calculate outliers
outliers <- boxplot.stats(pdat$Tamoxifen)$out

# Display t.test p-values on the plot
##################
# Ensure 'Group' is a factor and levels are in the desired order
IHC$Group <- factor(IHC$Group, levels = c("CRC AIS","CRLM"))

## 2) Plot boxplot and violin plot using ggplot2
## geom_violin for violin plot
library(ggplot2)
library(ggsignif)

ggplot(dat, aes(x = Group, y = `Quinolinic acid`, fill = Group)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.1, fill = "white") +
  scale_fill_manual(values = c("CRLM" = "red", "CRC" = "lightblue")) +
  geom_signif(comparisons = list(c("CRLM", "CRC")), 
              test = "wilcox.test", 
              map_signif_level = F) +
  theme_minimal() +
  theme(axis.line = element_line(color = "black"),  # Add axis lines for x and y axes
        axis.text.x = element_text(size = 10),      # Adjust x-axis text size
        axis.text.y = element_text(size = 10),      # Adjust y-axis text size
        axis.title.y = element_text(size = 14)) +   # Adjust y-axis title size
  ylab("QPRT") +
  xlab("Group: CRLM & CRC AIS")

######
## 3) Plot boxplot using geom_boxplot
ggplot(dat, aes(x = Group, y = `1,4-Dihydro-1-Methyl-4-Oxo-3-Pyridinecarboxamide`, fill = Group)) +
  geom_boxplot(width = 0.75, outlier.shape = NA) +
  scale_fill_manual(values = c("CRLM" = "red", "CRC" = "lightblue")) +
  geom_signif(comparisons = list(c("CRLM", "CRC")), 
              test = "wilcox.test", 
              map_signif_level = F) + 
  theme_minimal() +
  theme(axis.line = element_line(color = "black"),  # Add axis lines for x and y axes
        axis.text.x = element_text(size = 10),      # Adjust x-axis text size
        axis.text.y = element_text(size = 10),      # Adjust y-axis text size
        axis.title.y = element_text(size = 14)) +   # Adjust y-axis title size
  ylab("1,4-Dihydro-1-Methyl-4-Oxo-3-Pyridinecarboxamide") +
  xlab("Group: CRLM & CRC AIS")

################
# Load packages for calculating power and effect size
library(effsize)
library(pwr)
dat$`1,4-Dihydro-1-Methyl-4-Oxo-3-Pyridinecarboxamide`

# 1. Calculate effect size (Cohen's d)
# Using L-Tryptophan as an example
d_result <- cohen.d(`1,4-Dihydro-1-Methyl-4-Oxo-3-Pyridinecarboxamide` ~ Group, data = dat)
effect_size <- d_result$estimate
ci_low <- d_result$conf.int[1]
ci_high <- d_result$conf.int[2]
cat("Effect Size (Cohen's d):", effect_size, "\n")
cat("95% CI:", "[", ci_low, ",", ci_high, "]", "\n")

# 2. Calculate Post-hoc Power Analysis
# Use pwr.t2n.test for approximate estimation (or corresponding power transformation for Wilcoxon)
# n1=15, n2=15, d=calculated effect size, sig.level=0.05
power_val <- pwr.t2n.test(n1 = 15, n2 = 15, d = abs(effect_size), sig.level = 0.05)$power
cat("Statistical Power:", power_val, "\n")

# ==============================================================================
# Batch calculate P-values for all metabolites and perform FDR (Benjamini-Hochberg) multiple testing correction
# ==============================================================================
library(dplyr)

# 1. Extract column names of all metabolites (excluding the final "Group" column)
metabolites <- setdiff(colnames(dat), "Group")

# 2. Create an empty data frame to store results
results_df <- data.frame(
  Metabolite = character(),
  Raw_P_value = numeric(),
  stringsAsFactors = FALSE
)

# 3. Loop to calculate raw p-values for each metabolite
# Based on the methods described in the manuscript, Student's t-test is used here
for (metabolite in metabolites) {
  # Extract data for the two groups (Note: the control group might be "CRC" or "CRC AIS", handled for compatibility)
  crlm_data <- dat[dat$Group == "CRLM", metabolite]
  crc_data <- dat[dat$Group %in% c("CRC", "CRC AIS"), metabolite] 
  
  # Ensure numeric type
  crlm_data <- as.numeric(crlm_data)
  crc_data <- as.numeric(crc_data)
  
  # Run t-test (use tryCatch to prevent errors and interruption if some metabolites have zero abundance in all samples)
  test_res <- tryCatch({
    t.test(crlm_data, crc_data)
  }, error = function(e) {
    return(list(p.value = NA))
  })
  
  # Store results in the data frame
  results_df <- rbind(results_df, data.frame(
    Metabolite = metabolite,
    Raw_P_value = test_res$p.value
  ))
}

# Remove NA values from failed calculations
results_df <- results_df[!is.na(results_df$Raw_P_value), ]

# 4. Key step: Calculate FDR adjusted p-values (q-values)
# Use the common Benjamini-Hochberg (BH) method
results_df$FDR_P_value <- p.adjust(results_df$Raw_P_value, method = "BH")

# Sort by FDR values in ascending order
results_df <- results_df[order(results_df$FDR_P_value), ]

# ==============================================================================
# Extract FDR results of key metabolites focused on by reviewers in Table 6
# ==============================================================================
key_metabolites <- c("L-Tryptophan", 
                     "Nicotinamide", 
                     "N'-Methyl-2-pyridone-5-carboxamide", 
                     "1,4-Dihydro-1-Methyl-4-Oxo-3-Pyridinecarboxamide")

key_results <- results_df[results_df$Metabolite %in% key_metabolites, ]

print("====== FDR Results for Key Metabolites ======")
print(key_results)

# (Optional) Export FDR results of all metabolites to a CSV file for reviewer submission or supplementary materials
write.csv(results_df, "Metabolites_FDR_Results.csv", row.names = FALSE)
print("Complete FDR calculation results have been saved to 'Metabolites_FDR_Results.csv'")