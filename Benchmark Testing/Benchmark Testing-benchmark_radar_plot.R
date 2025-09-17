# If the package is not installed, run this line once:
# install.packages("fmsb")

# Load libraries
library(fmsb)
library(RColorBrewer)

# --- 1. Create data frames ---
# Internal test results
internal_data <- data.frame(
  row.names = c("Gradient Boosting", "XGBoost", "LightGBM", "Random Forest", "Support Vector Machine", "Logistic Regression", "AdaBoost", "K-Nearest Neighbors", "Decision Tree", "Gaussian Naive Bayes", "1D-CNN"),
  AUC = c(0.989872685, 0.988040123, 0.987750772, 0.970582562, 0.947145062, 0.929301698, 0.920814043, 0.830102238, 0.763888889, 0.612316744, 0.9174),
  Accuracy = c(0.930555556, 0.935185185, 0.930555556, 0.907407407, 0.847222222, 0.893518519, 0.828703704, 0.800925926, 0.800925926, 0.699074074, 0.8472),
  Precision = c(0.910828025, 0.933333333, 0.927152318, 0.902597403, 0.813559322, 0.906040268, 0.849673203, 0.798816568, 0.834437086, 0.725714286, 0.8581),
  Recall = c(0.993055556, 0.972222222, 0.972222222, 0.965277778, 1, 0.9375, 0.902777778, 0.9375, 0.875, 0.881944444, 0.9236)
)

# External test results
external_data <- data.frame(
  row.names = c("LightGBM", "Gradient Boosting", "AdaBoost", "XGBoost", "Random Forest", "Gaussian Naive Bayes", "K-Nearest Neighbors", "Decision Tree", "Support Vector Machine", "Logistic Regression", "1D-CNN"),
  AUC = c(0.947530864, 0.932098765, 0.865740741, 0.845679012, 0.557098765, 0.5, 0.472222222, 0.444444444, 0.361111111, 0.148148148, 0.9722),
  Accuracy = c(0.5, 0.5, 0.777777778, 0.666666667, 0.611111111, 0.5, 0.472222222, 0.444444444, 0.5, 0.25, 0.9167),
  Precision = c(0, 0, 1, 0.607142857, 0.642857143, 0, 0.485714286, 0.470588235, 0.5, 0.333333333, 0.9412),
  Recall = c(0, 0, 0.555555556, 0.944444444, 0.5, 0, 0.944444444, 0.888888889, 1, 0.5, 0.8889)
)

# Ensure the external test rows align with internal test order for plotting
external_data <- external_data[match(rownames(internal_data), rownames(external_data)), ]

# --- 2. Prepare data format required by radarchart ---
# radarchart expects the first row to be the maxima and the second row to be the minima.
# Since metrics are in [0,1], set max = 1 and min = 0.
max_min <- data.frame(
  AUC = c(1, 0),
  Accuracy = c(1, 0),
  Precision = c(1, 0),
  Recall = c(1, 0)
)

##############
# Option 1: Small multiple layout (recommended)
# Set plotting area to a 3x4 grid
par(mfrow = c(3, 4), mar = c(1, 1, 2, 1)) # adjust margins

# Loop through each model and draw a radar chart
for (i in 1:nrow(internal_data)) {
  model_name <- rownames(internal_data)[i]
  
  # Combine max/min, internal, and external rows for current model
  plot_data <- rbind(max_min, internal_data[i, ], external_data[i, ])
  
  # Draw radar chart
  radarchart(
    plot_data,
    axistype = 1, # axis style
    # polygons: internal test in blue, external in orange
    pcol = c("#0072B2", "#D55E00"),
    pfcol = scales::alpha(c("#0072B2", "#D55E00"), 0.3), # semi-transparent fills
    plwd = 2, # line width
    plty = 1, # line type
    # grid customization
    cglcol = "grey",
    cglty = 1,
    axislabcol = "grey",
    caxislabels = seq(0, 1, 0.25), # axis ticks
    cglwd = 0.8,
    # label customization
    vlcex = 0.8, # vertex label size
    title = model_name # chart title
  )
}

# Add a shared legend
plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
legend(
  x = "center",
  legend = c("Internal Test", "External Test"),
  bty = "n",
  pch = 20,
  col = c("#0072B2", "#D55E00"),
  cex = 1.5,
  pt.cex = 3
)

#########
# Option 2: Two combined plots for all 11 models (internal vs external)
# Prepare fmsb-compatible data frames
internal_fmsb <- rbind(max_min, internal_data)
external_fmsb <- rbind(max_min, external_data)

# Use 11 distinct colors
colors_11 <- brewer.pal(11, "Paired")

# Set plotting area to 1 row x 2 columns
par(mfrow = c(1, 2), mar = c(2, 2, 2, 2))

# --- Internal test radar plot ---
radarchart(
  internal_fmsb,
  axistype = 1,
  pcol = colors_11,
  plwd = 2,
  plty = 1,
  cglcol = "grey",
  cglty = 1,
  axislabcol = "grey",
  cglwd = 0.8,
  vlcex = 0.8,
  title = "Internal Test Results"
)

# --- External test radar plot ---
radarchart(
  external_fmsb,
  axistype = 1,
  pcol = colors_11,
  plwd = 2,
  plty = 1,
  cglcol = "grey",
  cglty = 1,
  axislabcol = "grey",
  cglwd = 0.8,
  vlcex = 0.8,
  title = "External Test Results"
)

# Add legend (positioning for multi-plot mode is non-trivial)
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("bottom", legend = rownames(internal_data), bty = "n", horiz = TRUE,
       fill = colors_11, cex = 0.7)