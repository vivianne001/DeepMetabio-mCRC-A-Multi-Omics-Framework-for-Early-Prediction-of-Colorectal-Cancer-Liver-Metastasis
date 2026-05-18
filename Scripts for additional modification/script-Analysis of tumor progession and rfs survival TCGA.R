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
library(stringr)

# 1. Assuming clinical data 'clin' and calculated 'trs_data' are already loaded
# clin <- read.csv("TCGA_clinical.csv")
# trs_data <- read.csv("TCGA_TRS_Scores.csv")
# df <- merge(clin, trs_data, by = "Sample_ID")

# Assuming you saved the Python prediction as a CSV, where the first column is sample_id (GSM...) and the second is TRS_Score
trs_data <- read.csv("predicted_TRS_scores_TCGA.csv") 
class(trs_data)
trs_data = trs_data[,c(1,2)]
samples = trs_data$Sample_ID

ph = ph[match(samples,ph$SAMPLE_ID),]
clin = clin[match(ph$PATIENT_ID,clin$PATIENT_ID),]

df <- clin # Assuming 'clin' already contains the "TRS_Score" column here

# For demonstration, merge trs_data with ph. Assume the merged dataframe is df
# df <- merge(ph, trs_data, by.x = "row.names", by.y = "Sample_ID")
df <- clin # Assuming 'ph' already contains the "TRS_Score" column
df <- clin 
df$TRS_Score = trs_data$TRS_Score

# ==============================================================================
# Step 1: Data preprocessing (Keep only initially non-metastatic M0 patients)
# ==============================================================================
# 1. Extract M0 patients, excluding M1, M1A, M1B, MX, etc., to meet the definition of prospective prediction
df_m0 <- df %>% filter(PATH_M_STAGE == "M0")

# 2. Clean AJCC Stage (Merge sub-stages, e.g., IIA/IIB/IIC -> II)
# Remove A/B/C suffixes, standardize to STAGE I, STAGE II, STAGE III
df_m0$Stage_clean <- gsub("[A-C]$", "", df_m0$AJCC_PATHOLOGIC_TUMOR_STAGE)

# Keep only explicit Stage I-III (further exclude potential Stage IV and NA)
df_m0 <- df_m0 %>% filter(Stage_clean %in% c("STAGE I", "STAGE II", "STAGE III"))
df_m0$Stage_clean <- factor(df_m0$Stage_clean, levels = c("STAGE I", "STAGE II", "STAGE III"))
df_m0$Stage_num <- as.numeric(df_m0$Stage_clean)

# 3. Process Progression-Free Survival (PFS) events and time
df_m0 <- df_m0 %>% filter(!is.na(PFS_STATUS) & !is.na(PFS_MONTHS) & PFS_MONTHS > 0)

# Set "1:PROGRESSION" to 1 (event occurred), "0:CENSORED" to 0 (censored)
df_m0$PFS_Event <- ifelse(grepl("1:PROGRESSION", df_m0$PFS_STATUS), 1, 0)
df_m0$PFS_Time <- as.numeric(as.character(df_m0$PFS_MONTHS))
df_m0$PFS_Time = df_m0$PFS_Time/12

# ==============================================================================
# Step 2: Trend analysis of TRS and TNM Stage (reflecting disease progression gradient)
# ==============================================================================
print("====== Analyzing the association between TRS scores and TNM Stage ======")

p_stage <- ggboxplot(df_m0, x = "Stage_clean", y = "TRS_Score",
                     color = "Stage_clean", palette = "jco",
                     add = "jitter", add.params = list(alpha=0.5)) +
  stat_compare_means(method = "kruskal.test", label.y = 1.2) +
  labs(title = "TRS across AJCC Stages (TCGA M0 Patients)",
       x = "AJCC Stage", y = "TRS Score") +
  scale_y_continuous(limits = c(0.5,1.2)) +
  theme_classic() + theme(legend.position = "none")
p_stage

# Spearman trend correlation analysis
cor_res <- cor.test(df_m0$Stage_num, df_m0$TRS_Score, method = "spearman", exact = FALSE)
print(paste("Stage vs TRS Correlation: rho =", round(cor_res$estimate, 3), "P =", signif(cor_res$p.value, 3)))

# ==============================================================================
# Step 3: Survival Analysis (PFS KM Curve)
# ==============================================================================
print("====== Performing Progression-Free Survival (PFS) analysis ======")

# Use survminer to find the optimal survival-related cutpoint
# 3) Calculate best cutoff using surv_cutpoint (prevent overwriting by median)
res.cut <- surv_cutpoint(
  data      = df_m0,
  time      = "PFS_Time",
  event     = "PFS_Event",
  variables = "TRS_Score",
  minprop   = 0.2
)
cutoff <- res.cut$cutpoint["TRS_Score", "cutpoint"]
print(paste("Best TRS cutoff =", round(cutoff, 4)))

# 4) Group by best cutoff
df_m0$Risk_Group <- ifelse(df_m0$TRS_Score >= cutoff, "High TRS", "Low TRS")
df_m0$Risk_Group <- factor(df_m0$Risk_Group, levels = c("Low TRS", "High TRS"))
table(df_m0$Risk_Group)

# 1. Kaplan-Meier (K-M) survival curve and Log-rank test
surv_obj <- Surv(time = df_m0$PFS_Time, event = df_m0$PFS_Event)
fit <- survfit(surv_obj ~ Risk_Group, data = df_m0)

p_surv <- ggsurvplot(fit,
                     data = df_m0,
                     pval = TRUE,             # Display Log-rank P-value
                     conf.int = TRUE,         # Display 95% confidence interval
                     risk.table = TRUE,       # Display risk table at the bottom (Number at risk)
                     palette = c("#0073C2FF", "#CD534CFF"), # Blue and red color palette
                     title = "Progression-Free Survival by TRS Stratification",
                     xlab = "Time (Months)",
                     ylab = "Progression-Free Survival Probability",
                     legend.title = "Risk Group")
print(p_surv)

##################--alternative
# (Optional) Perform univariate Cox regression to quantify Hazard Ratio (HR)
cox_res <- coxph(Surv(PFS_Time, PFS_Event) ~ Risk_Optimal, data = df_m0)
print(summary(cox_res))

df_tnm <- df_m0 %>%
  mutate(
    PFS_Time  = as.numeric(PFS_Time),
    PFS_Event = as.numeric(PFS_Event),
    TRS_Score = as.numeric(TRS_Score),
    AGE       = as.numeric(AGE),
    SEX       = factor(trimws(SEX)),
    PATH_T_STAGE = factor(trimws(PATH_T_STAGE)),
    PATH_N_STAGE = factor(trimws(PATH_N_STAGE))
  ) %>%
  filter(
    is.finite(PFS_Time), PFS_Time > 0,
    PFS_Event %in% c(0, 1),
    is.finite(TRS_Score), is.finite(AGE),
    !is.na(SEX), !is.na(PATH_T_STAGE), !is.na(PATH_N_STAGE), !is.na(PATH_M_STAGE)
  )

cox_multi_B <- coxph(Surv(PFS_Time, PFS_Event) ~ TRS_Score + AGE + SEX + PATH_T_STAGE + PATH_N_STAGE, data = df_tnm)
summary(cox_multi_B)