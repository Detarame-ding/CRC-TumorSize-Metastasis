##################################################################################################################
# Survival analysis was performed to compare long-term outcomes across tumor-size groups using Kaplan–Meier      #
# curves and log-rank tests.                                                                                     #
##################################################################################################################
rm(list=ls())
library(readxl)
library(survival)
library(survminer)
library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)

setwd('~/BS_Tumor')
SHCC_data <- readRDS("data/SHCC_data.RDS")

clean_ordered <- function(x, levels_vec) {
  x <- as.character(x)
  factor(x, levels = levels_vec, ordered = TRUE)
}
SHCC_data$Gender <- clean_ordered(SHCC_data$Gender, c("0","1"))
SHCC_data$pTstage <- clean_ordered(SHCC_data$pTstage, c("1", "2", "3", "4"))
SHCC_data$Location <- clean_ordered(SHCC_data$Location, c("1", "2", "3")) 
SHCC_data$Pathology_class <- clean_ordered(SHCC_data$Pathology_class, c("1", "2", "3")) 
SHCC_data$`Lymphovasscular invasion` <- clean_ordered(SHCC_data$`Lymphovasscular invasion`, c("0","1"))
SHCC_data$`Perineural invasion` <- clean_ordered(SHCC_data$`Perineural invasion`, c("0","1"))
SHCC_data$Differentiation <- clean_ordered(SHCC_data$Differentiation , c("1", "2", "3"))
SHCC_data$Metastasis <- clean_ordered(SHCC_data$Metastasis, c("0","1"))

cut_off_value <- 5

SHCC_data$TumorSizeGroup <- ifelse(
  !is.na(SHCC_data$max_diameter) & SHCC_data$max_diameter >= cut_off_value,
  "Large",
  ifelse(!is.na(SHCC_data$max_diameter) & SHCC_data$max_diameter < cut_off_value,
         "Small", NA)
)

SHCC_data$TumorSizeGroup <- factor(SHCC_data$TumorSizeGroup, levels = c("Small", "Large"))
df <- SHCC_data[!is.na(SHCC_data$TumorSizeGroup), ]

exclude_cols <- c(
  "ID", "pathID",
  "TumorSizeGroup",
  "max_diameter", 
  "TTPstatus", "TTPM",
  "OSstatus", "OSM",
  "DFSstatus", "DFSM",
  "TD",
  "Year",
  "EMR",
  "Location_minor",
  "Mucinous subtype",
  "CRM","pNstage","pMstage","TNMstage","livermeta","lungmeta","bonemeta","brainmeta",
  "ki67","PCNA","P21-1","P53-1","CD44-1","nm23","COX2","EGFR","KRAS_result_L1","NRAS_result_L1","BRAF_result_L1",
  "KRAS_result_L2","NRAS_result_L2","BRAF_result_L2","MLH1","MSH2","MSH6","PMS2","MMRstatus","MMRtest","MMRMSIstatus",
  "MMRMSItest","Tsize","Ttype","NAC"
)

vars_to_analyze <- setdiff(colnames(df), exclude_cols)

cols_needed <- c("TumorSizeGroup", vars_to_analyze)
df <- df[complete.cases(df[, cols_needed]), ]

df$Gender <- factor(df$Gender,
                    levels = c("0","1"),
                    labels = c("Female", "Male"))

df$pTstage <- factor(df$pTstage,
                     levels = c("1","2","3","4"),
                     labels = c("T1", "T2", "T3", "T4"))

df$Location <- factor(df$Location,
                      levels = c("1","2","3"),
                      labels = c("Rectum", "Distal colon", "Proximal colon"))

df$Pathology_class <- factor(df$Pathology_class,
                             levels = c("1","2","3"),
                             labels = c("Adenocarcinoma",
                                        "Adenocarcinoma with mucinous",
                                        "Adenocarcinoma with signet-ring cells"))

df$`Lymphovasscular invasion` <- factor(df$`Lymphovasscular invasion`,
                                        levels = c("0","1"),
                                        labels = c("No", "Yes"))

df$`Perineural invasion` <- factor(df$`Perineural invasion`,
                                   levels = c("0","1"),
                                   labels = c("No", "Yes"))

df$Differentiation <- factor(df$Differentiation,
                             levels = c("1","2","3"),
                             labels = c("Well", "Moderate", "Poor"))

df$Metastasis_num <- as.numeric(as.character(df$Metastasis))

surf_df <- df %>%
  filter(pTstage %in% c("T3", "T4"),
         max_diameter >= 0.5,
         max_diameter <= 12,
         Age < 50,
         Gender == 'Male'
  ) %>%
  mutate(
    size_group = case_when(
      max_diameter < 3 ~ "SmallT",
      max_diameter > 5 ~ "LargeT",
      TRUE ~ NA_character_
    ),
    TTPM = as.numeric(TTPM),
    TTPstatus = as.numeric(TTPstatus)
  ) %>%
  filter(!is.na(size_group), !is.na(TTPM), !is.na(TTPstatus)) %>%
  mutate(size_group = factor(size_group, levels = c("SmallT", "LargeT")))


fit <- survfit(Surv(TTPM, TTPstatus) ~ size_group, data = surf_df)

p <- ggsurvplot(
  fit,
  data = surf_df,
  pval = TRUE,
  risk.table = TRUE,
  xlim = c(0,120),
  xlab = "Time (months)",
  ylab = "Time to progress",
  legend.title = "",
)

ggsave(
  filename = "outs/TTP_KM_plot.pdf",
  plot = p$plot,
  width = 5,
  height = 5
)

library(dplyr)
library(survival)
library(survminer)
surf_df <- df %>%
  filter(pTstage %in% c("T3", "T4"),
         max_diameter >= 0.5,
         max_diameter <= 12,
         Age < 50,
  ) %>%
  mutate(
    size_group = case_when(
      max_diameter < 3 ~ "SmallT",
      max_diameter > 5 ~ "LargeT",
      TRUE ~ NA_character_
    ),
    TTPM = as.numeric(TTPM),
    TTPstatus = as.numeric(TTPstatus)
  ) %>%
  filter(!is.na(size_group), !is.na(TTPM), !is.na(TTPstatus), !is.na(Metastasis)) %>%
  mutate(size_group = factor(size_group, levels = c("SmallT", "LargeT")))

fit0 <- survfit(Surv(TTPM, TTPstatus) ~ size_group, data = filter(surf_df, Metastasis == 0))
fit1 <- survfit(Surv(TTPM, TTPstatus) ~ size_group, data = filter(surf_df, Metastasis == 1))

p0 <- ggsurvplot(
  fit0,
  data = filter(surf_df, Metastasis == 0),
  pval = TRUE,
  risk.table = TRUE,
  title = "Metastasis = 0",
  xlim = c(0,120),
  xlab = "Time (months)",
  ylab = "Time to progress",
  legend.title = ""
)

p1 <- ggsurvplot(
  fit1,
  data = filter(surf_df, Metastasis == 1),
  pval = TRUE,
  risk.table = TRUE,
  title = "Metastasis = 1",
  xlim = c(0,120),
  xlab = "Time (months)",
  ylab = "Time to progress",
  legend.title = ""
)

p0
p1
ggsave(
  filename = "outs/TTP_KM_plot-M0.pdf",
  plot = p0$plot,
  width = 5,
  height = 5
)
ggsave(
  filename = "outs/TTP_KM_plot-M1.pdf",
  plot = p1$plot,
  width = 5,
  height = 5
)

surf_df <- df %>%
  filter(pTstage %in% c("T3", "T4"),
         max_diameter >= 0.5,
         max_diameter <= 12,
         Age < 50,
         Gender == 'Male'
  ) %>%
  mutate(
    size_group = case_when(
      max_diameter < 3 ~ "SmallT",
      max_diameter > 5 ~ "LargeT",
      TRUE ~ NA_character_
    ),
    OSM = as.numeric(OSM),
    OSstatus = as.numeric(OSstatus)
  ) %>%
  filter(!is.na(size_group), !is.na(OSM), !is.na(OSstatus)) %>%
  mutate(size_group = factor(size_group, levels = c("SmallT", "LargeT")))


fit <- survfit(Surv(OSM, OSstatus) ~ size_group, data = surf_df)

p <- ggsurvplot(
  fit,
  data = surf_df,
  pval = TRUE,
  risk.table = TRUE,
  xlim = c(0,120),
  xlab = "Time (months)",
  ylab = "Overall Survival",
  legend.title = "",
)

print(p)
ggsave(
  filename = "outs/OS_KM_plot.pdf",
  plot = p$plot,
  width = 5,
  height = 5
)

library(dplyr)
library(survival)
library(survminer)

surf_df <- df %>%
  filter(pTstage %in% c("T3", "T4"),
         max_diameter >= 0.5,
         max_diameter <= 12,
         Age < 50
  ) %>%
  mutate(
    size_group = case_when(
      max_diameter < 3 ~ "SmallT",
      max_diameter > 5 ~ "LargeT",
      TRUE ~ NA_character_
    ),
    OSM = as.numeric(OSM),
    OSstatus = as.numeric(OSstatus)
  ) %>%
  filter(!is.na(size_group), !is.na(OSM), !is.na(OSstatus), !is.na(Metastasis)) %>%
  mutate(size_group = factor(size_group, levels = c("SmallT", "LargeT")))

fit0 <- survfit(Surv(OSM, OSstatus) ~ size_group, data = filter(surf_df, Metastasis == 0))
fit1 <- survfit(Surv(OSM, OSstatus) ~ size_group, data = filter(surf_df, Metastasis == 1))

p0 <- ggsurvplot(
  fit0,
  data = filter(surf_df, Metastasis == 0),
  pval = TRUE,
  risk.table = TRUE,
  title = "Metastasis = 0",
  xlim = c(0,120),
  xlab = "Time (months)",
  ylab = "Overall Survival",
  legend.title = ""
)

p1 <- ggsurvplot(
  fit1,
  data = filter(surf_df, Metastasis == 1),
  pval = TRUE,
  risk.table = TRUE,
  title = "Metastasis = 1",
  xlim = c(0,120),
  xlab = "Time (months)",
  ylab = "Overall Survival",
  legend.title = ""
)

print(p0)
print(p1)
ggsave(
  filename = "outs/OS_KM_plot-M0.pdf",
  plot = p0$plot,
  width = 5,
  height = 5
)
ggsave(
  filename = "outs/OS_KM_plot-M1.pdf",
  plot = p1$plot,
  width = 5,
  height = 5
)

surf_df <- df %>%
  filter(pTstage %in% c("T3", "T4"),
         max_diameter >= 0.5,
         max_diameter <= 12,
         Age < 50,
         Gender == 'Male'
  ) %>%
  mutate(
    size_group = case_when(
      max_diameter < 3 ~ "SmallT",
      max_diameter > 5 ~ "LargeT",
      TRUE ~ NA_character_
    ),
    DFSM = as.numeric(DFSM),
    DFSstatus = as.numeric(DFSstatus)
  ) %>%
  filter(!is.na(size_group), !is.na(DFSM), !is.na(DFSstatus)) %>%
  mutate(size_group = factor(size_group, levels = c("SmallT", "LargeT")))


fit <- survfit(Surv(DFSM, DFSstatus) ~ size_group, data = surf_df)

p <- ggsurvplot(
  fit,
  data = surf_df,
  pval = TRUE,
  risk.table = TRUE,
  xlim = c(0,120),
  xlab = "Time (months)",
  ylab = "Overall Survival",
  legend.title = "",
)

print(p)
ggsave(
  filename = "outs/DFS_KM_plot.pdf",
  plot = p$plot,
  width = 5,
  height = 5
)

library(dplyr)
library(survival)
library(survminer)

surf_df <- df %>%
  filter(pTstage %in% c("T3", "T4"),
         max_diameter >= 0.5,
         max_diameter <= 12,
         Age < 50
  ) %>%
  mutate(
    size_group = case_when(
      max_diameter < 3 ~ "SmallT",
      max_diameter > 5 ~ "LargeT",
      TRUE ~ NA_character_
    ),
    DFSM = as.numeric(DFSM),
    DFSstatus = as.numeric(DFSstatus)
  ) %>%
  filter(!is.na(size_group), !is.na(DFSM), !is.na(DFSstatus), !is.na(Metastasis)) %>%
  mutate(size_group = factor(size_group, levels = c("SmallT", "LargeT")))

fit0 <- survfit(Surv(DFSM, DFSstatus) ~ size_group, data = filter(surf_df, Metastasis == 0))
fit1 <- survfit(Surv(DFSM, DFSstatus) ~ size_group, data = filter(surf_df, Metastasis == 1))

p0 <- ggsurvplot(
  fit0,
  data = filter(surf_df, Metastasis == 0),
  pval = TRUE,
  risk.table = TRUE,
  title = "Metastasis = 0",
  xlim = c(0,120),
  xlab = "Time (months)",
  ylab = "Disease Free Survival",
  legend.title = ""
)

p1 <- ggsurvplot(
  fit1,
  data = filter(surf_df, Metastasis == 1),
  pval = TRUE,
  risk.table = TRUE,
  title = "Metastasis = 1",
  xlim = c(0,120),
  xlab = "Time (months)",
  ylab = "Disease Free Survival",
  legend.title = ""
)

print(p0)
print(p1)
ggsave(
  filename = "outs/DFS_KM_plot-M0.pdf",
  plot = p0$plot,
  width = 5,
  height = 5
)
ggsave(
  filename = "outs/DFS_KM_plot-M1.pdf",
  plot = p1$plot,
  width = 5,
  height = 5
)




