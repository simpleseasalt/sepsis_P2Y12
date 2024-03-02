# | ---------------------------------------
# | Author: Simplezzz
# | Date: 2024-01-19 19:30:25
# | LastEditTime: 2024-03-01 20:41:28
# | FilePath: \R_script\1_data_extraction.R
# | Description:
# | ---------------------------------------

# following packages are needed, install them before run this script.
# install.packages("data.table", "tidyverse", "magrittr", "tidymodels", "tidylog", "VIM", "naniar", "survival", "survminer", "forestmodel", "rstatix", "nortest", "MatchIt", "gtsummary")

library(data.table)
library(tidyverse)
library(magrittr)

tidymodels::tidymodels_prefer()

# 1----------------------------------------

setwd("D:/Physionet") # the directory of "mimic-iv-2.2" and "mimic_derived"

# The csv files under the "mimic_derived" folder were created by the sql scripts obtained from "https://github.com/MIT-LCP/mimic-code".

# 1---------------------------------------- comfirm patient with sepsis and antipltelet agents
# 2---------------------------------------- sepsis
sepsis <- fread("mimic-iv-2.2/mimic_derived/sepsis3.csv") %>%
    left_join(fread("mimic-iv-2.2/mimic_derived/icustay_detail.csv")) %>%
    left_join(select(fread("mimic-iv-2.2/hosp/admissions.csv"), hadm_id, deathtime)) %>%
    distinct(subject_id, hadm_id, stay_id, .keep_all = TRUE) %>%
    mutate(race = case_when(
        grepl("BLACK", race) ~ "Black",
        grepl("WHITE", race) ~ "White",
        grepl("ASIAN", race) ~ "Asian",
        grepl("HISPANIC", race) ~ "Hispanic",
        grepl("OTHER", race) ~ "Other"
    ))

# 2----------------------------------------

prescriptions <- fread("mimic-iv-2.2/hosp/prescriptions.csv")

# 3---------------------------------------- asipirin

aspirin_sorts <- prescriptions %>%
    distinct(drug) %>%
    filter(grepl("Aspirin", drug)) %>%
    mutate(drug_extracted = "aspirin")

aspirin <- prescriptions %>%
    filter(hadm_id %in% sepsis$hadm_id) %>%
    filter(drug %in% aspirin_sorts$drug) %>%
    left_join(select(sepsis, hadm_id, icu_intime)) %>%
    group_by(hadm_id) %>%
    filter(stoptime >= icu_intime & stoptime > starttime) %>% # some records had stoptime earlier than starttime
    ungroup() %>%
    mutate(dose_val_rx = as.numeric(dose_val_rx)) %>%
    group_by(hadm_id) %>%
    summarise(
        dose_min = min(dose_val_rx, na.rm = TRUE),
        dose_max = max(dose_val_rx, na.rm = TRUE),
        starttime = min(starttime),
        stoptime = max(stoptime)
    ) %>%
    mutate(
        duration_pres = ifelse(
            lubridate::hour(stoptime) >= 8,
            ceiling_date(stoptime, unit = "days") - floor_date(starttime, unit = "days"),
            ceiling_date(stoptime, unit = "days") - floor_date(starttime, unit = "days") - 1
        )
    )

# 3---------------------------------------- clopidogrel

clopidogrel_sorts <- prescriptions %>%
    distinct(drug) %>%
    filter(grepl("Clopidogrel", drug)) %>%
    mutate(drug_extracted = "clopidogrel")

clopidogrel <- prescriptions %>%
    filter(hadm_id %in% sepsis$hadm_id) %>%
    filter(drug %in% clopidogrel_sorts$drug) %>%
    left_join(select(sepsis, hadm_id, icu_intime)) %>%
    group_by(hadm_id) %>%
    filter(stoptime >= icu_intime & stoptime > starttime) %>% # some records had stoptime earlier than starttime
    ungroup() %>%
    mutate(dose_val_rx = as.numeric(dose_val_rx)) %>%
    group_by(hadm_id) %>%
    summarise(
        dose_min = min(dose_val_rx, na.rm = TRUE),
        dose_max = max(dose_val_rx, na.rm = TRUE),
        starttime = min(starttime),
        stoptime = max(stoptime)
    ) %>%
    mutate(
        duration_pres = ifelse(
            lubridate::hour(stoptime) >= 8,
            ceiling_date(stoptime, unit = "days") - floor_date(starttime, unit = "days"),
            ceiling_date(stoptime, unit = "days") - floor_date(starttime, unit = "days") - 1
        )
    )

# 3---------------------------------------- ticagrelor

ticagrelor_sorts <- prescriptions %>%
    distinct(drug) %>%
    filter(grepl("TiCAGRELOR", drug)) %>%
    mutate(drug_extracted = "ticagrelor")

ticagrelor <- prescriptions %>%
    filter(hadm_id %in% sepsis$hadm_id) %>%
    filter(drug %in% ticagrelor_sorts$drug) %>%
    left_join(select(sepsis, hadm_id, icu_intime)) %>%
    group_by(hadm_id) %>%
    filter(stoptime >= icu_intime & stoptime > starttime) %>% # some records had stoptime earlier than starttime
    ungroup() %>%
    mutate(dose_val_rx = as.numeric(dose_val_rx)) %>%
    group_by(hadm_id) %>%
    summarise(
        dose_min = min(dose_val_rx, na.rm = TRUE),
        dose_max = max(dose_val_rx, na.rm = TRUE),
        starttime = min(starttime),
        stoptime = max(stoptime)
    ) %>%
    mutate(
        duration_pres = ifelse(
            lubridate::hour(stoptime) >= 8,
            ceiling_date(stoptime, unit = "days") - floor_date(starttime, unit = "days"),
            ceiling_date(stoptime, unit = "days") - floor_date(starttime, unit = "days") - 1
        )
    ) %>%
    rename_at(vars(-c("hadm_id")), function(x) paste0(x, "_ticagrelor"))

# 2----------------------------------------  combine antiplatet agents

antiplatelet <- full_join(aspirin, clopidogrel, by = "hadm_id", suffix = c("_aspirin", "_clopidogrel")) %>%
    full_join(ticagrelor, by = "hadm_id")

# 2---------------------------------------- combine sic and antiplatet agents

target_patient <- sepsis %>%
    left_join(antiplatelet) %>%
    mutate(
        aspirin = ifelse(hadm_id %in% aspirin$hadm_id, 1, 0),
        clopidogrel = ifelse(hadm_id %in% clopidogrel$hadm_id, 1, 0),
        ticagrelor = ifelse(hadm_id %in% ticagrelor$hadm_id, 1, 0),
        antiplatelet_sorts = aspirin + clopidogrel + ticagrelor
    ) %>%
    tidylog::filter(aspirin == 1) %>%
    tidylog::filter(antiplatelet_sorts != 3) %>%
    mutate(group = ifelse(antiplatelet_sorts == 1, 0, 1))

# 1----------------------------------------

weight <- fread("mimic-iv-2.2/mimic_derived/first_day_weight.csv") %>%
    filter(stay_id %in% target_patient$stay_id) %>%
    select(stay_id, weight)

height <- fread("mimic-iv-2.2/mimic_derived/first_day_height.csv") %>%
    filter(stay_id %in% target_patient$stay_id)

# 1---------------------------------------- comorbidities and outcome

# Needs to read from uncompressed file to prevent automatically deleting the zeros before icd_code by Microsoft Excel.

diabetes_icd <- fread("mimic-iv-2.2/hosp/d_icd_diagnoses.csv.gz") %>%
    slice(2313:2372, 16339:16692)

hypertension_icd <- fread("mimic-iv-2.2/hosp/d_icd_diagnoses.csv.gz") %>%
    slice(4650:4682, 24499:24521)

neoplasm_icd <- fread("mimic-iv-2.2/hosp/d_icd_diagnoses.csv.gz") %>%
    slice(1275:2332, 13534:15627)

COPD_icd <- fread("mimic-iv-2.2/hosp/d_icd_diagnoses.csv.gz") %>%
    slice(26410:26414)

CA_surgery_icd <- fread("mimic-iv-2.2/hosp/d_icd_diagnoses.csv.gz") %>%
    slice(c(99666, 99667, 109654, 109658, 109749) - 1)

GI_icd <- fread("mimic-iv-2.2/hosp/d_icd_diagnoses.csv.gz") %>%
    slice(c(5063, 5065, 5592, 5607, 5608, 53120, 53121, 53200, 53201, 53220, 53221, 53300, 53301, 53320, 53321, 53400, 53401, 53420, 53421, 27628, 26044, 26047, 26899, 26917, 26925, 26927, 26935, 26937, 26945, 26947, 26978, 26981, 26994) - 1)

ICH_icd <- fread("mimic-iv-2.2/hosp/d_icd_diagnoses.csv.gz") %>%
    slice(4851:4855, 24823:24866)

VTE_icd <- fread("mimic-iv-2.2/hosp/d_icd_diagnoses.csv.gz") %>%
    slice(c(4740:4743, 4746, 7333:7337, 24627:24637, 5008:5025, 5038:5046, 25726:25841, 25925:25936) - 1)

CI_icd <- fread("mimic-iv-2.2/hosp/d_icd_diagnoses.csv.gz") %>%
    slice(c(4870, 4872, 24868:24983) - 1)

comorbidities <- fread("mimic-iv-2.2/hosp/diagnoses_icd.csv.gz") %>%
    filter(hadm_id %in% target_patient$hadm_id) %>%
    mutate(
        co_diabetes = ifelse(icd_code %in% diabetes_icd$icd_code, 1, 0),
        co_hypertension = ifelse(icd_code %in% hypertension_icd$icd_code, 1, 0),
        co_neoplasm = ifelse(icd_code %in% neoplasm_icd$icd_code, 1, 0),
        co_COPD = ifelse(icd_code %in% COPD_icd$icd_code, 1, 0),
        co_CA_surgery = ifelse(icd_code %in% CA_surgery_icd$icd_code, 1, 0),
        co_GI = ifelse(icd_code %in% GI_icd$icd_code, 1, 0),
        co_ICH = ifelse(icd_code %in% ICH_icd$icd_code, 1, 0),
        co_bleeding = if_else(co_GI == 1 | co_ICH == 1, 1, 0),
        co_VTE = ifelse(icd_code %in% VTE_icd$icd_code, 1, 0),
        co_CI = ifelse(icd_code %in% CI_icd$icd_code, 1, 0)
    ) %>%
    tidylog::filter(!(icd_code %in% GI_icd$icd_code & seq_num < 3)) %>%
    tidylog::filter(!(icd_code %in% ICH_icd$icd_code & seq_num < 3)) %>%
    tidylog::filter(!(icd_code %in% VTE_icd$icd_code & seq_num < 3)) %>%
    tidylog::filter(!(icd_code %in% CI_icd$icd_code & seq_num < 3)) %>%
    group_by(hadm_id) %>%
    summarise(
        co_diabetes = max(co_diabetes),
        co_hypertension = max(co_hypertension),
        co_neoplasm = max(co_neoplasm),
        co_COPD = max(co_COPD),
        co_CA_surgery = max(co_CA_surgery),
        co_GI = max(co_GI),
        co_ICH = max(co_ICH),
        co_bleeding = max(co_bleeding),
        co_VTE = max(co_VTE),
        co_CI = max(co_CI)
    )

# 1---------------------------------------- drugs
# 2---------------------------------------- heparin

heparin_sorts <- prescriptions %>%
    distinct(drug) %>%
    filter(grepl("Heparin", drug)) %>%
    filter(!grepl("Flush", drug))

heparin <- prescriptions %>%
    filter(hadm_id %in% target_patient$hadm_id) %>%
    filter(drug %in% heparin_sorts$drug) %>%
    left_join(select(target_patient, hadm_id, icu_intime)) %>%
    group_by(hadm_id) %>%
    filter(stoptime >= icu_intime & stoptime > starttime) %>% # some records had stoptime earlier than starttime
    ungroup() %>%
    mutate(dose_val_rx = as.numeric(dose_val_rx)) %>%
    group_by(hadm_id) %>%
    summarise(
        dose_min = min(dose_val_rx, na.rm = TRUE),
        dose_max = max(dose_val_rx, na.rm = TRUE),
        starttime = min(starttime),
        stoptime = max(stoptime)
    ) %>%
    mutate(
        duration_pres = ifelse(
            lubridate::hour(stoptime) >= 8,
            ceiling_date(stoptime, unit = "days") - floor_date(starttime, unit = "days"),
            ceiling_date(stoptime, unit = "days") - floor_date(starttime, unit = "days") - 1
        )
    ) %>%
    rename_at(vars(-c("hadm_id")), function(x) paste0(x, "_heparin"))

rm(prescriptions)

# 1---------------------------------------- baseline of lab result

coagulation_baseline <- fread("mimic-iv-2.2/mimic_derived/coagulation.csv") %>%
    mutate(
        starttime = floor_date(charttime, unit = "hour"),
        endtime = ceiling_date(charttime, unit = "hour")
    ) %>%
    select(-specimen_id) %>%
    rename("charttime_coagulation" = "charttime") %>%
    filter(hadm_id %in% target_patient$hadm_id) %>%
    group_by(hadm_id) %>%
    arrange(charttime_coagulation) %>%
    slice(1) %>%
    select(ends_with("id"), d_dimer:ptt)

blood_count_baseline <- fread("mimic-iv-2.2/mimic_derived/complete_blood_count.csv") %>%
    mutate(
        starttime = floor_date(charttime, unit = "hour"),
        endtime = ceiling_date(charttime, unit = "hour")
    ) %>%
    select(-specimen_id) %>%
    rename("charttime_blood" = "charttime") %>%
    filter(hadm_id %in% target_patient$hadm_id) %>%
    group_by(hadm_id) %>%
    arrange(charttime_blood) %>%
    slice(1) %>%
    select(ends_with("id"), hematocrit:wbc)

# 1---------------------------------------- dataset construction

data_raw <- target_patient %>%
    left_join(heparin) %>%
    left_join(height) %>%
    left_join(weight) %>%
    left_join(comorbidities) %>%
    left_join(coagulation_baseline) %>%
    left_join(blood_count_baseline) %>%
    tidylog::filter(admission_age >= 18) %>%
    tidylog::filter(los_icu > 2) %>%
    tidylog::filter(first_icu_stay == "t") %>%
    tidylog::filter(duration_pres_aspirin >= 3 | is.na(duration_pres_aspirin) | duration_pres_clopidogrel >= 3 | is.na(duration_pres_clopidogrel) | duration_pres_ticagrelor >= 3 | is.na(duration_pres_ticagrelor)) %>%
    ungroup() %>%
    tidylog::drop_na(co_diabetes:co_CA_surgery) %>%
    mutate(
        surv_time = as.double(as.period(deathtime - ymd_hm(sofa_time)), units = "days"),
        surv_30 = case_when(
            surv_time >= 30 ~ 30,
            surv_time == NA ~ 30,
            TRUE ~ surv_time
        ),
        surv_90 = case_when(
            surv_time >= 90 ~ 90,
            surv_time == NA ~ 90,
            TRUE ~ surv_time
        ),
        status_30 = ifelse(surv_30 == 30, 0, 1), # 0 is alive
        status_90 = ifelse(surv_90 == 90, 0, 1)
    ) %>%
    replace_na(
        list(
            surv_30 = 30,
            surv_90 = 90,
            status_30 = 0,
            status_90 = 0
        )
    )

# 1---------------------------------------- multiple imputation

library(mice)

baseline <- data_raw %>%
    select(
        Age = admission_age,
        Weight = weight,
        Height = height,
        Gender = gender,
        Race = race,
        Hypertension = co_hypertension,
        Diabetes = co_diabetes,
        Neoplasm = co_neoplasm,
        COPD = co_COPD,
        `History of CAS` = co_CA_surgery,
        INR = inr,
        PT = pt
    )

baseline %>%
    VIM::aggr(cex.label = 1.2, cex.aixs = 2, gap = 3, labels = names(.), oma = c(8, 6, 2, 3))

baseline %>%
    naniar::mcar_test()

data_simple_impute <- data_raw %>%
    replace_na(
        list(
            weight = median(.$weight, na.rm = T),
            inr = median(.$inr, na.rm = T),
            pt = median(.$pt, na.rm = T)
        )
    )

imp <- data_simple_impute %>%
    select(admission_age, weight, height, gender, race, co_diabetes, co_hypertension, co_neoplasm, co_COPD, co_CA_surgery, inr, pt) %>%
    mutate(race = as.factor(race)) %>%
    mice(method = "rf", m = 1, seed = 2024)

data_imputed <- complete(imp) %>%
    select(height, race) %>%
    rename_at(c(names(.)), function(x) paste0(x, "_imputed")) %>%
    as_tibble()

final_data <- data_simple_impute %>%
    bind_cols(data_imputed) %>%
    mutate(
        race = as.factor(race),
        gender = as.factor(gender),
        race = factor(race, levels = c("White", "Black", "Asian", "Hispanic", "Others")),
        co_diabetes = as.factor(co_diabetes),
        co_hypertension = as.factor(co_hypertension),
        co_neoplasm = as.factor(co_neoplasm),
        co_COPD = as.factor(co_COPD),
        co_CA_surgery = as.factor(co_CA_surgery),
        co_VTE = as.factor(co_VTE),
        co_CI = as.factor(co_CI),
        co_GI = as.factor(co_GI),
        co_ICH = as.factor(co_ICH),
        co_bleeding = as.factor(co_bleeding),
        group = as.factor(group)
    ) %>%
    replace_na(
        list(
            duration_pres_aspirin = 0,
            duration_pres_clopidogrel = 0,
            duration_pres_ticagrelor = 0,
            duration_pres_heparin = 0
        )
    )

save.image(file = "D:/OneDrive/After_work/2024/SIC_antiplatelet/output/sic_antiplatet.Rdata", compress = TRUE)

# ! end
