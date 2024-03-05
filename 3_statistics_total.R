# | ---------------------------------------
# | Author: Simplezzz
# | Date: 2024-02-06 15:48:59
# | LastEditTime: 2024-03-05 13:25:00
# | FilePath: d:\OneDrive\After_work\2024\sic_antiplatelet\github\3_statistics_total.R
# | Description:
# | ---------------------------------------

library(data.table)
library(tidyverse)
library(rstatix)
library(nortest)
library(gtsummary)

setwd(mimic_path)

load("sepsis_P2Y12/output/sepsis_P2Y12.Rdata")

# 1----------------------------------------

set.seed(2024)

data_tidy <- final_data %>%
    select(hadm_id, admission_age, gender, height_imputed, weight, starts_with("sofa"), race_imputed, starts_with("co_"), duration_pres_heparin, inr, pt, los_icu, los_hosp = los_hospital, surv_30, status_30, surv_90, status_90, group, -sofa_time) %>%
    as_tibble()

# 1---------------------------------------- numeric variable

numeric_var <- data_tidy %>%
    select_if(is.numeric) %>%
    select(-c(hadm_id, status_30, status_90, surv_30, surv_90)) %>%
    names()

# 2---------------------------------------- normality

ad.test.multi <- function(x) {
    ad.test(x) %>%
        broom::tidy()
}

normality <- map_df(
    data_tidy[numeric_var],
    ad.test.multi
) %>%
    bind_cols(variable = numeric_var) %>%
    select(-method) %>%
    rename("p_norm" = "p.value")

# 2---------------------------------------- variance test

variance <- map_df(
    data_tidy[numeric_var],
    function(x) {
        levene_test(data_tidy, x ~ data_tidy$group)
    }
) %>%
    select(p_vari = p) %>%
    cbind(variable = numeric_var)

# 2---------------------------------------- t test

t_test_res <- data_tidy %>%
    select(all_of(numeric_var), group) %>%
    pivot_longer(
        cols = all_of(numeric_var),
        names_to = "variable"
    ) %>%
    group_by(variable) %>%
    t_test(value ~ group, var.equal = TRUE) %>%
    select(variable, p_t_test = p)

# 2---------------------------------------- wilcox test

wilcox_test_res <- data_tidy %>%
    select(all_of(numeric_var), group) %>%
    pivot_longer(
        cols = all_of(numeric_var),
        names_to = "variable"
    ) %>%
    group_by(variable) %>%
    wilcox_test(value ~ group) %>%
    select(variable, p_wilcox_test = p)

p_numb <- normality %>%
    left_join(variance, by = "variable") %>%
    left_join(t_test_res, by = "variable") %>%
    left_join(wilcox_test_res, by = "variable") %>%
    group_by(variable) %>%
    mutate(p_final = if_else(p_norm >= 0.05 & p_vari >= 0.05, p_t_test, p_wilcox_test)) %>%
    add_significance(p.col = "p_final")

# 1----------------------------------------

nominal_var <- data_tidy %>%
    select_if(is.factor) %>%
    select(-c(group)) %>%
    names()

p_chisq <- data_tidy %>%
    select(all_of(nominal_var), group) %>%
    pivot_longer(
        cols = -group,
        names_to = "variable"
    ) %>%
    group_by(variable) %>%
    do(chisq_test(.$group, .$value)) %>%
    select(variable, p_chisq = p)

# 2---------------------------------------- fisher

freq <- freq_table(data_tidy, group, all_of(nominal_var))

fisher_map <- function(x) {
    temp <- freq_table(data_tidy, group, x) %>%
        select(-prop) %>%
        pivot_wider(
            names_from = "group",
            values_from = "n"
        ) %>%
        replace(is.na(.), 0)
    p_fisher <- temp %>%
        select(all_of(unique(data_tidy$group))) %>%
        fisher_test(simulate.p.value = TRUE)
    min_n <- temp %>%
        pivot_longer(
            cols = all_of(unique(data_tidy$group))
        ) %>%
        arrange(value) %>%
        slice(1) %>%
        select(value) %>%
        cbind(variable = x)
    cbind(p_fisher, min_n)
}

p_fisher <- map_df(
    all_of(nominal_var),
    fisher_map
) %>%
    select(p_fisher = p, min_n = value, variable)

p_norm <- tibble(variable = all_of(nominal_var)) %>%
    left_join(p_chisq, by = "variable") %>%
    left_join(p_fisher, by = "variable") %>%
    group_by(variable) %>%
    mutate(p_final = ifelse(nrow(data_tidy) > 40 & min_n >= 5, p_chisq, p_fisher))

p_numb
p_norm

p_value <- bind_rows(
    p_numb %>%
        select(variable, p_final),
    p_norm %>%
        select(variable, p_final)
)

# 1---------------------------------------- get summary

library(gtsummary)

theme_gtsummary_language("en", big.mark = "")

psm_summary <- data_tidy %>%
    tbl_summary(
        by = group,
        statistic = list(
            all_continuous2() ~ c("{mean} \u00B1 {sd}", "{median} ({p25}, {p75})", "{p_miss}"),
            all_categorical() ~ "{n} ({p}%)"
        ),
        digits = list(
            all_continuous() ~ 1,
            all_categorical() ~ 1
        )
    )

psm_summary %>%
    as_tibble() %>%
    rename("variable" = `**Characteristic**`) %>%
    left_join(p_value) %>%
    mutate(p_final = ifelse(p_final >= 0.001, round(p_final, 3), "< 0.001")) %>%
    write_csv(file = "sepsis_P2Y12/output/total_summary.csv")

# ! end