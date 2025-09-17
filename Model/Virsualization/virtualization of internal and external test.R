
rm(list = ls())
###################################
# Install and load required packages
# If a package is not installed, uncomment the install.packages() line and run it once.
# install.packages("PRROC")
# install.packages("ggplot2")
library(PRROC)
library(ggplot2)

# ------------------------------------------------------------------
# Load prediction results (internal or external validation)
# ------------------------------------------------------------------
# Choose one of the lines below to load your dataset.
# The expected columns used in this script: 
#  - probability column(s): e.g. "Pred_Probability" or "Predicted_Probability"
#  - true label: e.g. "True_Label" or "True_Label_int"
#  - predicted label: e.g. "Predicted_Label_int"
# Adjust column names in the code if your CSV uses different names.
#
# Example usage (uncomment the appropriate line):
# dat <- read.csv("test_predictions_internal_with_labels.csv")  # internal validation result
# dat <- read.csv("test_predictions_external_with_labels.csv")  # external validation result

# For demonstration, load one file (edit the filename as needed)
dat <- read.csv("test_predictions_internal_with_labels.csv")
dat <- as.data.frame(dat)

# ------------------------------------------------------------------
# Precision-Recall (PR) curve and AUPRC
# ------------------------------------------------------------------
# NOTE: PRROC expects scores for one class and the corresponding weights/labels.
# Adjust the column names below to match your data.
pr <- pr.curve(scores.class0 = dat$Pred_Probability, weights.class0 = dat$True_Label, curve = TRUE)

# Convert PR curve to a data frame for ggplot
pr_df <- data.frame(Recall = pr$curve[, 1], Precision = pr$curve[, 2])

# Compute AUPRC
auprc <- pr$auc.integral

# Plot PR curve with customized theme
ggplot(pr_df, aes(x = Recall, y = Precision)) +
  geom_line(size = 1.2, color = "#0072B5") + # main PR curve
  geom_hline(yintercept = 0, color = "black", size = 1) +
  geom_vline(xintercept = 0, color = "black", size = 1) +
  theme_minimal(base_size = 16) +
  labs(
    title = "Precision-Recall Curve",
    subtitle = paste("AUPRC =", round(auprc, 3)),
    x = "Recall",
    y = "Precision"
  ) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.title = element_text(face = "bold"),
    panel.grid.major = element_line(color = "#e0e0e0"),
    panel.grid.minor = element_blank()
  ) +
  annotate("text", x = 0.7, y = 0.2, label = paste("AUPRC =", round(auprc, 3)), size = 6, color = "#E18727")

# ------------------------------------------------------------------
# Confusion matrix visualization
# ------------------------------------------------------------------
# Build confusion table.
# Ensure your data frame has integer-coded labels (e.g., 0/1 or 1/2) in the columns used below.
cm <- table(True = dat$True_Label_int, Pred = dat$Predicted_Label_int)

# Optional: set readable row/column names (modify if your labels differ)
rownames(cm) <- c("Primary", "Metastasis")
colnames(cm) <- c("Primary", "Metastasis")

# Convert table to data frame for plotting
cm_df <- as.data.frame(as.table(cm))
colnames(cm_df) <- c("True_Label", "Pred_Label", "Freq")

# Plot confusion matrix as a heat-tile
ggplot(cm_df, aes(x = Pred_Label, y = True_Label, fill = Freq)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", Freq)), size = 10) +
  scale_fill_gradient(low = "#FBDFE2", high = "#B83945") +
  labs(
    title = "Confusion Matrix",
    x = "Predicted Label",
    y = "True Label"
  ) +
  theme_minimal(base_size = 18) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold")
  )

# ------------------------------------------------------------------
# ROC curve and AUROC (with 95% CI)
# ------------------------------------------------------------------
# Load package for ROC analysis
# install.packages("pROC") # uncomment if not installed
library(pROC)

# Compute ROC object.
# Make sure to use the correct probability column name from your data.
roc_obj <- roc(dat$True_Label_int, dat$Predicted_Probability)

# Compute AUC and 95% confidence interval
auc_val <- auc(roc_obj)
ci_auc <- ci.auc(roc_obj)

# Create ROC plotting data frame
roc_df <- data.frame(
  FPR = 1 - roc_obj$specificities,
  TPR = roc_obj$sensitivities
)

# Plot ROC curve
ggplot(roc_df, aes(x = FPR, y = TPR)) +
  geom_line(size = 1.2, color = "#0072B5") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black", size = 1) +
  geom_hline(yintercept = 0, color = "black", size = 1) +
  geom_vline(xintercept = 0, color = "black", size = 1) +
  theme_minimal(base_size = 16) +
  labs(
    title = "ROC Curve",
    x = "False Positive Rate (FPR)",
    y = "True Positive Rate (TPR)"
  ) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold")
  ) +
  annotate(
    "text",
    x = 0.6, y = 0.2,
    label = paste0("AUROC = ", round(auc_val, 3),
                   "\n95% CI: [", round(ci_auc[1], 3), ", ", round(ci_auc[3], 3), "]"),
    size = 6, color = "#0072B5"
  )