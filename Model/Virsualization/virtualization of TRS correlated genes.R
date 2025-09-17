
rm(list = ls())
# Load required packages
# If a package is missing, uncomment and run the install.packages() line once.
# install.packages("ggplot2")
library(ggplot2)
library(scales)

# ------------------------------------------------------------------
# Load data
# ------------------------------------------------------------------
# Provide one of the lines below to load the dataset you intend to plot.
# Expected input columns:
#   - "Gene"        : gene symbol/name
#   - "Correlation" : correlation coefficient with TRS (numeric)
#   - "P_Value"     : p-value for the correlation
#
# Example usage (uncomment the appropriate line and set correct filepath):
# dat <- read.csv("trs_gene_correlations.csv", header = TRUE)                # training & internal validation
# dat <- read.csv("trs_gene_correlation_in_CRLM_samples.csv", header = TRUE) # external validation (CRLM)
dat <- read.csv("trs_gene_correlations.csv", header = TRUE)
dat <- as.data.frame(dat)

# ------------------------------------------------------------------
# Generate significance stars from p-values
# ------------------------------------------------------------------
# Stars convention: *** < 0.001, ** < 0.01, * < 0.05, otherwise blank
dat$Signif <- cut(dat$P_Value,
                  breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                  labels = c("***", "**", "*", ""))

# ------------------------------------------------------------------
# Select top genes by correlation:
# take top 20 positively correlated and top 20 negatively correlated genes
# ------------------------------------------------------------------
# Sort and select top positive and negative correlations
top_pos <- dat[order(-dat$Correlation), ][1:20, ]
top_neg <- dat[order(dat$Correlation), ][1:20, ]
top_genes <- rbind(top_pos, top_neg)

# Assign colors: red for positive, blue for negative
top_genes$Color <- ifelse(top_genes$Correlation > 0, "#c23616", "#192a56")

# Preserve gene order by correlation for plotting (from low to high)
top_genes$Gene <- factor(top_genes$Gene, levels = top_genes$Gene[order(top_genes$Correlation)])

# ------------------------------------------------------------------
# Plot horizontal bar chart of correlations with significance labels
# ------------------------------------------------------------------
# Note: adjust text size and y-limits if necessary for your data
y_min <- min(top_genes$Correlation, na.rm = TRUE) - 0.1
y_max <- max(top_genes$Correlation, na.rm = TRUE) + 0.1

ggplot(top_genes, aes(x = Gene, y = Correlation, fill = Color)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_text(aes(label = sprintf("%.3f %s", Correlation, Signif)),
            hjust = ifelse(top_genes$Correlation > 0, -0.1, 1.1),
            color = "black", size = 5) +
  scale_fill_identity() +
  # Draw a vertical reference line at 0 (no correlation)
  geom_vline(xintercept = 0, color = "black", size = 1) +
  coord_flip() +
  labs(
    title = "Top 20 Positive and Top 20 Negative Genes Correlated with TRS",
    x = "Gene",
    y = "Correlation Coefficient (with TRS)"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold"),
    axis.line = element_line(color = "black", size = 1)
  ) +
  ylim(y_min, y_max)