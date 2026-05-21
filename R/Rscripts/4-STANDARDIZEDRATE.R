##################################################################################################################
# Standardized metastasis rates were calculated to adjust for differences in clinicopathological characteristics #
# across tumor-size groups.                                                                                      #
##################################################################################################################

rm(list = ls())

library(readxl)
library(survival)
library(survminer)
library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)
library(patchwork)

setwd('~/BS_Tumor')
SHCC_data <- readRDS("data/SHCC_data.RDS")

clean_ordered <- function(x, levels_vec) {
  x <- as.character(x)
  factor(x, levels = levels_vec, ordered = TRUE)
}

SHCC_data$Gender <- clean_ordered(SHCC_data$Gender, c("0","1"))
SHCC_data$pTstage <- clean_ordered(SHCC_data$pTstage, c("1", "2", "3", "4"))
SHCC_data$Location <- clean_ordered(SHCC_data$Location, c("1", "2", "3")) # 1:直肠 2:远端 3:近端
SHCC_data$Pathology_class <- clean_ordered(SHCC_data$Pathology_class, c("1", "2", "3")) # 1:腺癌 2:腺癌伴粘液 3:腺癌伴印戒细胞
SHCC_data$`Lymphovasscular invasion` <- clean_ordered(SHCC_data$`Lymphovasscular invasion`, c("0","1"))
SHCC_data$`Perineural invasion` <- clean_ordered(SHCC_data$`Perineural invasion`, c("0","1"))
SHCC_data$Differentiation <- clean_ordered(SHCC_data$Differentiation, c("1", "2", "3"))
SHCC_data$Metastasis <- clean_ordered(SHCC_data$Metastasis, c("0","1"))

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

exclude_cols <- intersect(exclude_cols, colnames(SHCC_data))
vars_to_analyze <- setdiff(colnames(SHCC_data), exclude_cols)

cols_needed <- c("max_diameter", "Metastasis", vars_to_analyze)
cols_needed <- intersect(cols_needed, colnames(SHCC_data))

df0 <- SHCC_data[complete.cases(SHCC_data[, cols_needed]), ]

df0$Metastasis_num <- as.numeric(as.character(df0$Metastasis))

df0 <- df0 %>%
  filter(!is.na(max_diameter), !is.na(Metastasis_num))


std_vars <- c(
  "Gender",
  "pTstage",
  "Location",
  "Differentiation",
  "Pathology_class"
)

std_vars <- intersect(std_vars, colnames(df0))

cutoff_seq <- seq(1, 10, by = 0.1)


p_to_star <- function(p) {
  dplyr::case_when(
    is.na(p) ~ "NA",
    p < 0.001 ~ "***",
    p < 0.01  ~ "**",
    p < 0.05  ~ "*",
    TRUE ~ "ns"
  )
}

format_p <- function(p) {
  dplyr::case_when(
    is.na(p) ~ "P = NA",
    p < 0.001 ~ "P < 0.001",
    TRUE ~ paste0("P = ", sprintf("%.3f", p))
  )
}


calc_rate_by_cutoff <- function(data, cutoff_value, std_vars) {
  
  group_levels <- tibble(
    TumorSizeGroup = factor(c("Small", "Large"), levels = c("Small", "Large"))
  )
  
  df <- data %>%
    mutate(
      TumorSizeGroup = ifelse(max_diameter >= cutoff_value, "Large", "Small"),
      TumorSizeGroup = factor(TumorSizeGroup, levels = c("Small", "Large"))
    ) %>%
    select(TumorSizeGroup, Metastasis_num, all_of(std_vars)) %>%
    filter(complete.cases(.))
  
  crude_rate <- df %>%
    group_by(TumorSizeGroup) %>%
    summarise(
      N = n(),
      Metastasis_n = sum(Metastasis_num == 1),
      Crude_rate = mean(Metastasis_num == 1),
      .groups = "drop"
    ) %>%
    right_join(group_levels, by = "TumorSizeGroup") %>%
    mutate(
      N = replace_na(N, 0L),
      Metastasis_n = replace_na(Metastasis_n, 0L)
    ) %>%
    arrange(TumorSizeGroup)
  
  std_pop <- df %>%
    group_by(across(all_of(std_vars))) %>%
    summarise(std_n = n(), .groups = "drop") %>%
    mutate(std_weight = std_n / sum(std_n))
  
  stratum_rate <- df %>%
    group_by(TumorSizeGroup, across(all_of(std_vars))) %>%
    summarise(
      n = n(),
      events = sum(Metastasis_num == 1),
      rate = mean(Metastasis_num == 1),
      .groups = "drop"
    )
  
  std_rate <- crossing(
    TumorSizeGroup = factor(c("Small", "Large"), levels = c("Small", "Large")),
    std_pop
  ) %>%
    left_join(stratum_rate, by = c("TumorSizeGroup", std_vars)) %>%
    left_join(
      crude_rate %>% select(TumorSizeGroup, Crude_rate),
      by = "TumorSizeGroup"
    ) %>%
    mutate(
      rate_filled = ifelse(is.na(rate), Crude_rate, rate)
    ) %>%
    group_by(TumorSizeGroup) %>%
    summarise(
      Standardized_rate = ifelse(
        all(is.na(rate_filled)),
        NA_real_,
        sum(rate_filled * std_weight, na.rm = TRUE)
      ),
      .groups = "drop"
    )
  
  crude_rate %>%
    left_join(std_rate, by = "TumorSizeGroup") %>%
    mutate(cutoff = cutoff_value) %>%
    select(
      cutoff,
      TumorSizeGroup,
      N,
      Metastasis_n,
      Crude_rate,
      Standardized_rate
    )
}


rate_5cm <- calc_rate_by_cutoff(
  data = df0,
  cutoff_value = 5,
  std_vars = std_vars
)

print(rate_5cm)

write.csv(
  rate_5cm,
  "outs/metastasis_rate_5cm_crude_standardized.csv",
  row.names = FALSE
)


df_5cm <- df0 %>%
  mutate(
    TumorSizeGroup = ifelse(max_diameter >= 5, "Large", "Small"),
    TumorSizeGroup = factor(TumorSizeGroup, levels = c("Small", "Large")),
    Metastasis_num = factor(Metastasis_num, levels = c(0, 1))
  ) %>%
  select(TumorSizeGroup, Metastasis_num, all_of(std_vars)) %>%
  filter(complete.cases(.))


tab_5cm <- table(df_5cm$TumorSizeGroup, df_5cm$Metastasis_num)

p_crude_5cm <- suppressWarnings(
  chisq.test(tab_5cm, correct = FALSE)$p.value
)

fit_full <- glm(
  as.numeric(as.character(Metastasis_num)) ~ TumorSizeGroup + Gender + pTstage + Location + Differentiation + Pathology_class,
  data = df_5cm,
  family = binomial()
)

fit_null <- glm(
  as.numeric(as.character(Metastasis_num)) ~ Gender + pTstage + Location + Differentiation + Pathology_class,
  data = df_5cm,
  family = binomial()
)

p_std_5cm <- anova(fit_null, fit_full, test = "Chisq")$`Pr(>Chi)`[2]

p_5cm_df <- tibble(
  RateType = factor(
    c("Crude rate", "Standardized rate"),
    levels = c("Crude rate", "Standardized rate")
  ),
  Method = c(
    "Pearson chi-square test",
    "Adjusted logistic regression likelihood-ratio test"
  ),
  P_value = c(p_crude_5cm, p_std_5cm),
  P_label = format_p(P_value),
  Star = p_to_star(P_value)
)

print(p_5cm_df)

write.csv(
  p_5cm_df,
  "outs/metastasis_rate_5cm_p_values.csv",
  row.names = FALSE
)


plot_5cm_df <- rate_5cm %>%
  pivot_longer(
    cols = c(Crude_rate, Standardized_rate),
    names_to = "RateType",
    values_to = "Rate"
  ) %>%
  mutate(
    RateType = factor(
      RateType,
      levels = c("Crude_rate", "Standardized_rate"),
      labels = c("Crude rate", "Standardized rate")
    )
  )

bracket_df <- plot_5cm_df %>%
  group_by(RateType) %>%
  summarise(
    y = max(Rate, na.rm = TRUE) * 1.15,
    .groups = "drop"
  ) %>%
  left_join(p_5cm_df, by = "RateType") %>%
  mutate(
    label = paste0(Star, "\n", P_label),
    x1 = 1,
    x2 = 2,
    xm = 1.5
  )

y_max_5cm <- max(bracket_df$y, na.rm = TRUE) * 1.18

p_5cm <- ggplot(
  plot_5cm_df,
  aes(x = TumorSizeGroup, y = Rate, fill = TumorSizeGroup)
) +
  geom_col(width = 0.62) +
  geom_text(
    aes(label = percent(Rate, accuracy = 0.1)),
    vjust = -0.45,
    size = 4
  ) +
  geom_segment(
    data = bracket_df,
    aes(x = x1, xend = x2, y = y, yend = y),
    inherit.aes = FALSE,
    linewidth = 0.6
  ) +
  geom_segment(
    data = bracket_df,
    aes(x = x1, xend = x1, y = y * 0.985, yend = y),
    inherit.aes = FALSE,
    linewidth = 0.6
  ) +
  geom_segment(
    data = bracket_df,
    aes(x = x2, xend = x2, y = y * 0.985, yend = y),
    inherit.aes = FALSE,
    linewidth = 0.6
  ) +
  geom_text(
    data = bracket_df,
    aes(x = xm, y = y * 1.04, label = label),
    inherit.aes = FALSE,
    size = 4.2,
    lineheight = 0.9
  ) +
  facet_wrap(~ RateType, nrow = 1) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    limits = c(0, y_max_5cm)
  ) +
  scale_fill_manual(values = c(
    "Small" = "#72bcd5",
    "Large" = "#f7aa58"
  )) +
  labs(
    title = "Metastasis rate by 5-cm tumor-size cut-off",
    subtitle = "Crude rate was compared by Pearson chi-square test; standardized rate was compared by adjusted logistic regression",
    x = "Tumor size group",
    y = "Metastasis rate",
    fill = "Tumor size group"
  ) +
  theme_bw(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "top",
    strip.background = element_rect(fill = "grey90", color = "grey70"),
    strip.text = element_text(face = "bold")
  )

print(p_5cm)

ggsave(
  filename = "outs/metastasis_rate_5cm_crude_standardized.pdf",
  plot = p_5cm,
  width = 6,
  height = 8,
  dpi = 300
)


rate_all <- bind_rows(
  lapply(
    cutoff_seq,
    function(x) calc_rate_by_cutoff(df0, x, std_vars)
  )
)

write.csv(
  rate_all,
  "outs/metastasis_rate_by_cutoff_crude_standardized.csv",
  row.names = FALSE
)


plot_df <- rate_all %>%
  pivot_longer(
    cols = c(Crude_rate, Standardized_rate),
    names_to = "RateType",
    values_to = "Rate"
  ) %>%
  mutate(
    RateType = factor(
      RateType,
      levels = c("Crude_rate", "Standardized_rate"),
      labels = c("Crude rate", "Standardized rate")
    )
  )


p_rate <- ggplot(
  plot_df,
  aes(
    x = cutoff,
    y = Rate,
    color = TumorSizeGroup,
    linetype = RateType,
    group = interaction(TumorSizeGroup, RateType)
  )
) +
  geom_smooth(
    se = FALSE,
    method = "loess",
    span = 0.25,
    linewidth = 1.1
  ) +
  scale_x_continuous(limits = c(1, 10), expand = c(0, 0)) +
  scale_y_continuous(labels = percent_format(accuracy = 0.1)) +
  scale_color_manual(values = c(
    "Small" = "#72bcd5",
    "Large" = "#f7aa58"
  )) +
  scale_linetype_manual(values = c(
    "Crude rate" = "solid",
    "Standardized rate" = "dashed"
  )) +
  labs(
    title = "Crude and standardized metastasis rate by tumor-size cut-off",
    subtitle = "Direct standardization by Gender, pTstage, Location, Differentiation, and Pathology",
    x = "Tumor size cut-off",
    y = "Metastasis rate",
    color = "Tumor size group",
    linetype = "Rate type"
  ) +
  theme_bw(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "top"
  )

print(p_rate)

ggsave(
  filename = "outs/metastasis_rate_by_cutoff_crude_standardized.pdf",
  plot = p_rate,
  width = 10,
  height = 6,
  dpi = 300
)


p_n <- ggplot(
  rate_all,
  aes(
    x = cutoff,
    y = N,
    color = TumorSizeGroup,
    group = TumorSizeGroup
  )
) +
  geom_line(linewidth = 1.1) +
  scale_x_continuous(limits = c(1, 10), expand = c(0, 0)) +
  scale_color_manual(values = c(
    "Small" = "#72bcd5",
    "Large" = "#f7aa58"
  )) +
  labs(
    x = "Tumor size cut-off",
    y = "Sample size",
    color = "Tumor size group"
  ) +
  theme_bw(base_size = 13) +
  theme(
    legend.position = "none"
  )

print(p_n)

ggsave(
  filename = "outs/sample_size_by_cutoff.pdf",
  plot = p_n,
  width = 10,
  height = 4,
  dpi = 300
)


p_rate_top <- p_rate +
  labs(x = NULL) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.margin = margin(t = 5.5, r = 5.5, b = 0, l = 5.5)
  )

p_n_bottom <- p_n +
  theme(
    plot.margin = margin(t = 0, r = 5.5, b = 5.5, l = 5.5)
  )

p_combined <- p_rate_top / p_n_bottom +
  plot_layout(heights = c(3, 1), guides = "collect") &
  theme(legend.position = "top")

print(p_combined)

ggsave(
  filename = "outs/metastasis_rate_by_cutoff_crude_standardized_with_sample_size.pdf",
  plot = p_combined,
  width = 6,
  height = 8,
  dpi = 300
)