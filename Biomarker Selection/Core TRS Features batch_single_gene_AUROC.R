rm(list = ls())
######################## Compare AUC values of candidate features ####################

dat = as.data.frame(t(exp))
dat$group = ph$`tissue:ch1`
dat$group = ifelse(dat$group == "metastatic colorectal cancer to the liver", "metastasis", "primary")
table(dat$group)

#############################################
library(ROCR)
library(pROC)
library(ggplot2)
library(colorspace)

# Prepare binary labels: metastasis = 1, primary = 0
labels <- as.numeric(dat$group == "metastasis")

# Genes to evaluate (ordered list)
metabolic_genes = c("CYP2C9","CYP2C8","CYP2A6","ADH1A","CYP1A2","RDH16","HSD17B6","ALDH8A1","DHTKD1","ACMSD")
n_genes <- length(metabolic_genes)

# Generate a qualitative color palette (supports up to 20 genes)
my_colors <- qualitative_hcl(n_genes, palette = "Dark 3")
names(my_colors) <- metabolic_genes

# Swap colors for ACMSD and ALDH8A1 (visual preference)
tmp <- my_colors["ACMSD"]
my_colors["ACMSD"] <- my_colors["ALDH8A1"]
my_colors["ALDH8A1"] <- tmp

# Compute AUC and 95% CI for each gene, and collect ROC curve data
auc_info <- data.frame(
  gene = character(),
  auc = numeric(),
  ci_low = numeric(),
  ci_high = numeric(),
  stringsAsFactors = FALSE
)
roc_data_list <- list()
for (gene in metabolic_genes) {
  pred <- prediction(dat[[gene]], labels)
  perf <- performance(pred, "tpr", "fpr")
  auc <- performance(pred, "auc")@y.values[[1]]
  ci <- ci.auc(roc(labels, dat[[gene]]))
  roc_data_list[[gene]] <- data.frame(
    fpr = perf@x.values[[1]],
    tpr = perf@y.values[[1]],
    group = gene
  )
  auc_info <- rbind(auc_info, data.frame(
    gene = gene,
    auc = auc,
    ci_low = ci[1],
    ci_high = ci[3]
  ))
}
roc_data <- do.call(rbind, roc_data_list)

# Order genes by descending AUC
auc_info <- auc_info[order(-auc_info$auc), ]

# Prepare annotation text for AUC and 95% CI
auc_texts <- paste0(
  auc_info$gene, " AUC: ", round(auc_info$auc, 2), " [", 
  round(auc_info$ci_low, 3), "-", round(auc_info$ci_high, 3), "]"
)

# Plot ROC curves; legend order follows AUC descending
p <- ggplot(roc_data, aes(x = fpr, y = tpr, color = group)) +
  geom_line(size = 1.2) +
  geom_abline(linetype = "dashed") +
  theme_minimal() +
  theme(
    axis.line = element_line(color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size = 12)
  ) +
  labs(title = "ROC curves", x = "False Positive Rate (1 - Specificity)", y = "Sensitivity") +
  scale_color_manual(values = my_colors, breaks = auc_info$gene) # legend ordered by AUC

# Add AUC annotation text on the plot (ordered by AUC)
for (i in seq_along(auc_texts)) {
  p <- p + annotate("text", x = 0.65, y = 0.45 - 0.045 * (i - 1), label = auc_texts[i], hjust = 0, size = 4)
}

print(p)