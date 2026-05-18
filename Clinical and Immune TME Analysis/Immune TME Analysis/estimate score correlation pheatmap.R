
rm(list = ls())
#############################
library(ggplot2)
library(dplyr)
library(scales)
library(readxl)

# Required packages
library(ggplot2)
library(dplyr)

# --- Load dataset (ensure `dat` has the expected structure) ---
dat = read_xlsx("estimate.results.xlsx", col_names = TRUE)
dat = as.data.frame(dat)

# 1. Data preprocessing: set factor levels and generate significance labels from p-values
dat$group <- factor(dat$group, levels = c("CRC-PR", "CRLM-L", "CRLM-H")) # set x-axis order
dat$Variable <- factor(dat$Variable,
                       levels = rev(c("Tumor Purity Score", "Immune Score", "Stromal Score", "Estimate Score"))) # set y-axis order (reversed for visual preference)

dat <- dat %>%
  mutate(
    signif_label = case_when(
      pvalue < 0.001 ~ "***",
      pvalue < 0.01  ~ "**",
      pvalue < 0.05  ~ "*",
      TRUE           ~ ""  # blank when not significant
    )
  )

# Set color limits for correlation values
cor_limits <- c(-0.55, 0.59)

# 2. Plotting: correlation dotplot with significance annotation
p_final <- ggplot(dat, aes(x = group, y = Variable)) +
  # Points: color represents correlation value
  geom_point(aes(color = Cor), size = 8) +
  
  # Add significance stars
  geom_text(aes(label = signif_label), vjust = -1, size = 5, color = "black") +
  
  # Color scale: blue-white-red with midpoint at 0
  scale_color_gradient2(low = "#2166ac", mid = "white", high = "#b2182b",
                        midpoint = 0, limits = cor_limits, oob = scales::squish, name = "Correlation") +
  
  # Theme and axis adjustments
  theme_minimal(base_size = 12) +
  theme(
    axis.title = element_blank(), # hide axis titles
    axis.text.x = element_text(angle = 45, hjust = 1), # tilt x-axis labels
    panel.grid = element_blank(), # remove grid lines
    axis.line = element_line(colour = "black") # add axis lines
  ) +
  labs(title = "Correlation dotplot (color = Cor)")

# 3. Display the final plot
print(p_final)