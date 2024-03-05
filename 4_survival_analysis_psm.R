# | ---------------------------------------
# | Author: Simplezzz
# | Date: 2024-02-22 16:25:26
# | LastEditTime: 2024-03-01 20:21:16
# | FilePath: \R_script\4_survival_analysis_psm.R
# | Description:
# | ---------------------------------------

library(tidyverse)
library(data.table)
library(survival)
library(survminer)
library(forestmodel)

setwd(mimic_path)

load("sepsis_P2Y12/output/sepsis_P2Y12.Rdata")

# 1---------------------------------------- survival analysis
# 2---------------------------------------- 30 days

survfit_psm_30 <- survfit(Surv(surv_30, status_30) ~ group, data = data_psm)

surv_pvalue(survfit_psm_30, method = "survdiff")$pval

surv_curve_psm_30 <- ggsurvplot(
    survfit_psm_30,
    conf.int = TRUE,
    pval = "Log-rank p = 0.094",
    pval.size = 6,
    pval.coord = c(0, 0.5),
    pval.method = TRUE,
    pval.method.size = 6,
    pval.method.coord = c(0, 0.5),
    legend.title = "Group",
    legend.labs = c("Aspirin", "Aspirin + P2Y12"),
    xlab = "Time/days",
    font.x = c(18, "bold"),
    font.y = c(18, "bold"),
    font.tickslab = c(14, "plain"),
    font.legend = c(16),
    risk.table = "absolute",
    risk.table.pos = "in",
    risk.table.title = "Number in risk",
    fontsize = 6,
    ggtheme = theme_classic()
)

table_y_limit <- 0.18

surv_curve_psm_30$plot <- surv_curve_psm_30$plot +
    annotate("text", x = 0, y = table_y_limit - 0.03, label = "Percentage in risk", color = "black", size = 6, hjust = "inward") +
    geom_hline(aes(yintercept = table_y_limit), linetype = 2) +
    annotate("text", x = 25, y = 0.77, label = "Survival probability = 84.6%", color = "#00bfc4", size = 6) +
    annotate("text", x = 25, y = 0.95, label = "Survival probability = 86.5%", color = "#f8766d", size = 6)

tiff(filename = "plot/psm/surv_curve_psm_30.tiff", width = 10, height = 6, res = 300, units = "in", compression = "lzw")

surv_curve_psm_30

dev.off()

# 1---------------------------------------- cox
# 2---------------------------------------- psm

cox_fit_psm <- data_psm %>%
    mutate(
        gender = factor(gender, levels = c("F", "M"), labels = c("Female", "Male")),
        co_diabetes = factor(co_diabetes, levels = c(0, 1), labels = c("No", "Yes")),
        co_hypertension = factor(co_hypertension, levels = c(0, 1), labels = c("No", "Yes")),
        co_neoplasm = factor(co_neoplasm, levels = c(0, 1), labels = c("No", "Yes")),
        co_COPD = factor(co_COPD, levels = c(0, 1), labels = c("No", "Yes")),
        group = factor(group, levels = c(0, 1), labels = c("No", "Yes")),
    ) %>%
    rename(
        "Age" = "admission_age",
        "Gender" = "gender",
        "Weight" = "weight",
        "Height" = "height_imputed",
        "SOFA score" = "sofa_score",
        "INR" = "inr",
        "PT" = "pt",
        "Diabetes" = "co_diabetes",
        "Hypertension" = "co_hypertension",
        "Neoplasm" = "co_neoplasm",
        "COPD" = "co_COPD",
        "Duration of heparin" = "duration_pres_heparin",
        "P2Y12 Inhibitor" = "group"
    ) %>%
    coxph(Surv(surv_30, status_30) ~ Age + Gender + Height + Weight + `SOFA score` + Diabetes + Hypertension + Neoplasm + COPD + `Duration of heparin` + INR + PT + `P2Y12 Inhibitor`, data = .)

panels <- list(
    list(width = 0.03),
    list(width = 0.1, display = ~variable, fontface = "bold", heading = "Variable"),
    list(width = 0.1, display = ~level),
    list(width = 0.05, display = ~n, hjust = 1, heading = "N"),
    list(width = 0.03, item = "vline", hjust = 0.5),
    list(
        width = 0.50, item = "forest", hjust = 0.5, heading = "Hazard ratio", linetype = "dashed",
        line_x = 0
    ),
    list(width = 0.03, item = "vline", hjust = 0.5),
    list(width = 0.15, heading = "95% CI", display = ~ ifelse(reference, "Reference", sprintf(
        "%0.2f (%0.2f, %0.2f)",
        trans(estimate), trans(conf.low), trans(conf.high)
    )), display_na = NA),
    list(
        width = 0.05,
        display = ~ ifelse(reference, "", format.pval(p.value, digits = 1, eps = 0.001)),
        display_na = NA, hjust = 1, heading = "P value"
    ),
    list(width = 0.03)
)

coxplot_psm <- forest_model(
    cox_fit_psm,
    format_options = forest_model_format_options(
        point_size = 4
    ),
    panels = panels,
    theme = theme_void()
)

tiff(filename = "plot/coxplot_psm.tiff", width = 10, height = 6, res = 300, units = "in", compression = "lzw")

coxplot_psm

dev.off()

coxplot_psm$data %>%
    write_csv(file = "sepsis_P2Y12/output/cox_result_psm.csv")

# 2---------------------------------------- total

cox_fit_total <- final_data %>%
    mutate(
        gender = factor(gender, levels = c("F", "M"), labels = c("Female", "Male")),
        co_diabetes = factor(co_diabetes, levels = c(0, 1), labels = c("No", "Yes")),
        co_hypertension = factor(co_hypertension, levels = c(0, 1), labels = c("No", "Yes")),
        co_neoplasm = factor(co_neoplasm, levels = c(0, 1), labels = c("No", "Yes")),
        co_COPD = factor(co_COPD, levels = c(0, 1), labels = c("No", "Yes")),
        group = factor(group, levels = c(0, 1), labels = c("No", "Yes")),
    ) %>%
    rename(
        "Age" = "admission_age",
        "Gender" = "gender",
        "Weight" = "weight",
        "Height" = "height_imputed",
        "SOFA score" = "sofa_score",
        "INR" = "inr",
        "PT" = "pt",
        "Diabetes" = "co_diabetes",
        "Hypertension" = "co_hypertension",
        "Neoplasm" = "co_neoplasm",
        "COPD" = "co_COPD",
        "Duration of heparin" = "duration_pres_heparin",
        "P2Y12 Inhibitor" = "group"
    ) %>%
    coxph(Surv(surv_30, status_30) ~ Age + Gender + Height + Weight + `SOFA score` + Diabetes + Hypertension + Neoplasm + COPD + `Duration of heparin` + INR + PT + `P2Y12 Inhibitor`, data = .)

coxplot_total <- forest_model(
    cox_fit_total,
    format_options = forest_model_format_options(
        point_size = 4
    ),
    panels = panels,
    theme = theme_void()
)

tiff(filename = "plot/coxplot_total.tiff", width = 10, height = 6, res = 300, units = "in", compression = "lzw")

coxplot_total

dev.off()

coxplot_total$data %>%
    write_csv(file = "sepsis_P2Y12/output/cox_result_total.csv")

save.image(file = "sepsis_P2Y12/output/sepsis_P2Y12.Rdata", compress = TRUE)

# ! end
