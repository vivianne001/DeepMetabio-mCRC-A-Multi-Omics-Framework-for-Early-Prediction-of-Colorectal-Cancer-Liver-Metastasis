rm(list = ls())

# --- 0. Install and load required packages ---
# If you have not installed these packages yet, run the install.packages() lines once.
# install.packages("ggplot2")
# install.packages("patchwork")
# install.packages("ggrepel")

library(ggplot2)
library(patchwork)
library(ggrepel)

# --- 1. Prepare data and colors ---
# (Data and color definitions are included here for convenience; adjust as needed)
# Internal test data
internal_data <- data.frame(
  row.names = c("Gradient Boosting", "XGBoost", "LightGBM", "Random Forest", "Support Vector Machine", "Logistic Regression", "AdaBoost", "K-Nearest Neighbors", "Decision Tree", "Gaussian Naive Bayes", "1D-CNN"),
  AUC = c(0.989872685, 0.988040123, 0.987750772, 0.970582562, 0.947145062, 0.929301698, 0.920814043, 0.830102238, 0.763888889, 0.612316744, 0.9059),
  Accuracy = c(0.930555556, 0.935185185, 0.930555556, 0.907407407, 0.847222222, 0.893518519, 0.828703704, 0.800925926, 0.800925926, 0.699074074, 0.8472),
  Precision = c(0.910828025, 0.933333333, 0.927152318, 0.902597403, 0.813559322, 0.906040268, 0.849673203, 0.798816568, 0.834437086, 0.725714286, 0.8725),
  Recall = c(0.993055556, 0.972222222, 0.972222222, 0.965277778, 1, 0.9375, 0.902777778, 0.9375, 0.875, 0.881944444, 0.9028)
)

# External test data
external_data <- data.frame(
  row.names = c("LightGBM", "Gradient Boosting", "AdaBoost", "XGBoost", "Random Forest", "Gaussian Naive Bayes", "K-Nearest Neighbors", "Decision Tree", "Support Vector Machine", "Logistic Regression", "1D-CNN"),
  AUC = c(0.947530864, 0.932098765, 0.865740741, 0.845679012, 0.557098765, 0.5, 0.472222222, 0.444444444, 0.361111111, 0.148148148, 0.9136),
  Accuracy = c(0.5, 0.5, 0.777777778, 0.666666667, 0.611111111, 0.5, 0.472222222, 0.444444444, 0.5, 0.25, 0.8611),
  Precision = c(0, 0, 1, 0.607142857, 0.642857143, 0, 0.485714286, 0.470588235, 0.5, 0.333333333, 1),
  Recall = c(0, 0, 0.555555556, 0.944444444, 0.5, 0, 0.944444444, 0.888888889, 1, 0.5, 0.7222)
)

# Ensure external data rows align with internal data order
external_data <- external_data[match(rownames(internal_data), rownames(external_data)), ]

# Ensure data frames contain a "Model" column rather than only row names for ggplot compatibility
internal_data$Model <- rownames(internal_data)
external_data$Model <- rownames(external_data)

# --- Define high-contrast colors for models ---
custom_colors <- c(
  "Gradient Boosting" = "#FFD700", "XGBoost" = "#377EB8", "LightGBM" = "#4DAF4A",
  "Random Forest" = "#984EA3", "Support Vector Machine" = "#FF7F00", "Logistic Regression" = "#F781BF",
  "AdaBoost" = "#A65628", "K-Nearest Neighbors" = "#000000", "Decision Tree" = "#999999",
  "Gaussian Naive Bayes" = "#00CED1", "1D-CNN" = "#E41A1C"
)

# Prepare colors and short labels
model_names <- rownames(internal_data)
plot_colors <- custom_colors[model_names]
short_names_vec <- c("GB", "XGB", "LGBM", "RF", "SVM", "LR", "Ada", "KNN", "DT", "GNB", "1D-CNN")
internal_data$ShortName <- short_names_vec
external_data$ShortName <- short_names_vec

# --- 2. Axis ranges ---
# Set x-axis ranges: internal tests narrower to emphasize differences; external tests full range
internal_x_range <- c(0.6, 1.0)
external_x_range <- c(0.1, 1.0)
y_range <- c(0.2, 1.0)

# --- 3. Internal test scatter plot ---
p_internal <- ggplot(internal_data, aes(x = AUC, y = Accuracy, color = Model)) +
  geom_point(size = 5, alpha = 0.9) +
  # Optional: add text labels with ggrepel (commented out for cleaner figure)
  # geom_text_repel(aes(label = ShortName), size = 3)
  scale_color_manual(values = custom_colors) +
  coord_cartesian(xlim = internal_x_range, ylim = y_range) +
  scale_x_continuous(breaks = seq(0.6, 1.0, 0.1)) +
  scale_y_continuous(breaks = seq(0.2, 1.0, 0.1)) +
  labs(title = "Internal Test Results",
       subtitle = "Better models appear toward the top-right",
       x = "Area Under Curve (AUC)",
       y = "Accuracy") +
  theme_bw() +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5, color = "gray50"),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 10),
    legend.position = "bottom",
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 9),
    panel.grid.major = element_line(alpha = 0.3),
    panel.grid.minor = element_line(alpha = 0.1)
  )

# --- 4. External test scatter plot ---
p_external <- ggplot(external_data, aes(x = AUC, y = Accuracy, color = Model)) +
  geom_point(size = 5, alpha = 0.9) +
  # geom_text_repel(aes(label = ShortName), size = 3)  # optional
  scale_color_manual(values = custom_colors) +
  coord_cartesian(xlim = external_x_range, ylim = y_range) +
  scale_x_continuous(breaks = seq(0.1, 1.0, 0.1)) +
  scale_y_continuous(breaks = seq(0.2, 1.0, 0.1)) +
  labs(title = "External Test Results",
       subtitle = "Performance degradation is indicated by a shift toward the bottom-left",
       x = "Area Under Curve (AUC)",
       y = "Accuracy") +
  theme_bw() +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5, color = "gray50"),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 10),
    legend.position = "bottom",
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 9),
    panel.grid.major = element_line(alpha = 0.3),
    panel.grid.minor = element_line(alpha = 0.1)
  )

# --- 5. Combine plots with patchwork ---
final_plot <- p_internal + p_external +
  plot_layout(guides = 'collect') &
  theme(
    legend.position = 'bottom',
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10),
    legend.key.size = unit(0.8, "cm")
  )

# --- 6. Display the final combined figure ---
print(final_plot)