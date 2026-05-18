# Clear environment
rm(list=ls())

# ==============================================================================
# Step 0: Load necessary R packages and data
# ==============================================================================
library(survival)
library(survminer)
library(ggplot2)
library(ggpubr)
library(dplyr)

# 1. Load clinical data and expression matrix of GSE39582
# Note: Local directory path translated to English for SCI submission
load("D:/CRC_Liver_Metastasis_Biomarker_Diagnosis/Github_Public/data_preprocess/Other_datasets/GSE39582/GSE39582_Raw_Data.Rdata")

# 2. Load the TRS score data calculated via Python 1D-CNN
# Assuming the prediction was saved as CSV in Python, with the first column as sample_id (GSM...) and the second as TRS_Score
trs_data <- read.csv("predicted_TRS_scores_gse39582.csv") 
class(trs_data)
trs_data = trs_data[,c(1,2)]

# For demonstration, merge trs_data with ph. Assume the merged data frame is df
# df <- merge(ph, trs_data, by.x = "row.names", by.y = "sample_id")
identical(trs_data$Sample_ID, ph$geo_accession)

# (Here ph is used instead of df for continuous demonstration, assuming ph already contains the "TRS_Score" column)
# Note: Before running, ensure that ph already contains the continuous predicted probability values "TRS_Score"
df <- ph 
df$TRS_Score = trs_data$TRS_Score

# ==============================================================================
# Step 1: Strict data filtering (excluding M1 samples with distant metastasis)
# ==============================================================================
# Keep only M0 (no distant metastasis) patients with valid clinical information
# This meets the reviewer's requirement for "initially non-metastatic" patients
df_m0 <- df %>% filter(`tnm.m:ch1` == "M0")

# ==============================================================================
# Step 2: Extract follow-up endpoint indicators
# ==============================================================================
# ⚠️ Note: Please find the actual follow-up time and event column names for GSE39582 in colnames(df_m0)!
# Follow-up data for GSE39582 is usually hidden in characteristics_ch1.xx, for example:
# characteristics_ch1.xx: "relapse-free survival (mes): 45.2" 
# characteristics_ch1.xx: "relapse (1=true): 0"

# Variable placeholders are used here, please replace them with the actual column names you found:
time_col <- "RFS_Time_Months"  # Replace with the actual Recurrence-Free Survival time column name
event_col <- "RFS_Event"       # Replace with the actual recurrence event column name (usually 0=no recurrence, 1=recurrence)

# Remove samples with missing follow-up data
# df_m0 <- df_m0[!is.na(df_m0[[time_col]]) & !is.na(df_m0[[event_col]]), ]
# df_m0$Time <- as.numeric(df_m0[[time_col]])
# df_m0$Event <- as.numeric(df_m0[[event_col]])

# ==============================================================================
# Added Step: Trend and correlation analysis between TRS Score and TNM Stage
# ==============================================================================
print("====== Analyzing the association between TRS Score and TNM Stage ======")

# Extract Stage 1, 2, 3 patients (excluding Stage 0 and 4)
df_stage <- df_m0 %>% 
  filter(`tnm.t:ch1` %in% c("T1", "T2", "T3","T4"))

# Format factor variables
df_stage$T_stage <- factor(df_stage$`tnm.t:ch1`, levels = c("T1", "T2", "T3","T4"), labels = c("T1", "T2","T3","T4"))

# 1. Draw boxplot and perform Kruskal-Wallis test for inter-group differences
p_stage <- ggboxplot(df_stage, x = "T_stage", y = "TRS_Score",
                     color = "T_stage", palette = c("#0073C2FF", "#EFC000FF", "#CD534CFF","orange"),
                     add = "jitter", add.params = list(alpha=0.5)) +
  stat_compare_means(method = "kruskal.test", label.y = 0.29) +
  labs(title = "TRS Distribution Across TNM Stages (M0 Patients)",
       x = "TNM Stage", y = "Transcriptomic Risk Score (TRS)") +
  scale_y_continuous(limits = c(0, 0.3)) +
  theme_classic() +
  theme(legend.position = "none")

print(p_stage)
# Save plot option: ggsave("TRS_vs_Stage_Boxplot.pdf", p_stage, width = 5, height = 5)

# 2. Trend correlation analysis (Spearman Correlation)
df_stage$T_stage <- factor(df_stage$`tnm.t:ch1`, levels = c("T1", "T2", "T3","T4"), labels = c(1,2,3,4))
df_stage$T_stage <- as.numeric(df_stage$T_stage)

cor_res <- cor.test(df_stage$T_stage, df_stage$TRS_Score, method = "spearman")
print(paste("Spearman Correlation between Stage and TRS: rho =", 
            round(cor_res$estimate, 3), ", p-value =", signif(cor_res$p.value, 3)))

# ==============================================================================
# Step 3: Survival Analysis
# ==============================================================================
print("====== Performing RFS/DFS survival analysis ======")

# 1. Kaplan-Meier (K-M) survival curve and Log-rank test
surv_obj <- Surv(time = df_m0$`rfs.delay:ch1`, event = df_m0$`rfs.event:ch1`)
fit <- survfit(surv_obj ~ Risk_Group, data = df_m0)

p_surv <- ggsurvplot(fit,
                     data = df_m0,
                     pval = TRUE,             # Display Log-rank P-value
                     conf.int = TRUE,         # Display 95% confidence interval
                     risk.table = TRUE,       # Display risk table at the bottom (Number at risk)
                     palette = c("#0073C2FF", "#CD534CFF"), # Blue and red color palette
                     title = "Recurrence-Free Survival by TRS Stratification",
                     xlab = "Time (Months)",
                     ylab = "Recurrence-Free Survival Probability",
                     legend.title = "Risk Group")

print(p_surv)
# Save plot option: ggsave("KM_Survival_Curve.pdf", print(p_surv), width = 6, height = 6)

# 2. Univariate Cox proportional hazards regression
cox_uni <- coxph(surv_obj ~ TRS_Score, data = df_m0)
print("--- Univariate Cox Regression (Continuous TRS) ---")
print(summary(cox_uni))

cox_uni_cat <- coxph(surv_obj ~ Risk_Group, data = df_m0)
print("--- Univariate Cox Regression (Categorical High vs Low) ---")
print(summary(cox_uni_cat))

# 3. Multivariate Cox proportional hazards regression (verifying if TRS is an independent prognostic factor)
# Combining age, gender, Stage, chemotherapy, etc. as covariates
cox_multi <- coxph(surv_obj ~ Risk_Group + `tnm.stage:ch1` + `chemotherapy.adjuvant:ch1` + `age.at.diagnosis (year):ch1`, data = df_m0)

# print("--- Multivariate Cox Regression ---")
print(summary(cox_multi))

print("Analysis completed!")