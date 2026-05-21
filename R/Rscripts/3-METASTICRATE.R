##################################################################################################################
# Metastasis rates were calculated across tumor-size groups, including lymph node metastasis,                    #
# distant metastasis, and overall metastasis.                                                                    #
##################################################################################################################
rm(list = ls())

library(dplyr)
library(ggplot2)
library(scales)
library(patchwork)
library(splines)

SHCC_data <- readRDS("data/SHCC_data_cleaned.RDS")

clean_ordered <- function(x, levels_vec) {
  x <- as.character(x)
  factor(x, levels = levels_vec, ordered = TRUE)
}

SHCC_data$Gender <- clean_ordered(SHCC_data$Gender, c("0", "1"))
SHCC_data$pTstage <- clean_ordered(SHCC_data$pTstage, c("1", "2", "3", "4"))
SHCC_data$Location <- clean_ordered(SHCC_data$Location, c("1", "2", "3")) 
SHCC_data$Pathology_class <- clean_ordered(SHCC_data$Pathology_class, c("1", "2", "3"))
SHCC_data$`Lymphovasscular invasion` <- clean_ordered(SHCC_data$`Lymphovasscular invasion`, c("0", "1"))
SHCC_data$`Perineural invasion` <- clean_ordered(SHCC_data$`Perineural invasion`, c("0", "1"))
SHCC_data$Differentiation <- clean_ordered(SHCC_data$Differentiation, c("1", "2", "3"))

to_numeric_safe <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x[x %in% c("", "NA", "N/A", "na", "n/a", "x", "X", ".", "-", "unknown", "Unknown")] <- NA
  suppressWarnings(as.numeric(x))
}

round_to_half <- function(x) {
  round(x * 2) / 2
}

candidate_distant_cols <- c(
  "livermeta",
  "lungmeta",
  "bonemeta",
  "brainmeta"
)

distant_cols_present <- intersect(candidate_distant_cols, colnames(SHCC_data))


print(distant_cols_present)

pN_num <- to_numeric_safe(SHCC_data$pNstage)
pM_num <- to_numeric_safe(SHCC_data$pMstage)

if (length(distant_cols_present) > 0) {
  distant_mat <- sapply(SHCC_data[, distant_cols_present, drop = FALSE], to_numeric_safe)
  if (is.vector(distant_mat)) {
    distant_mat <- matrix(distant_mat, ncol = 1)
    colnames(distant_mat) <- distant_cols_present[1]
  }
  distant_any <- apply(distant_mat, 1, function(v) {
    if (all(is.na(v))) {
      NA
    } else {
      any(v != 0, na.rm = TRUE)
    }
  })
} else {
  distant_any <- rep(NA, nrow(SHCC_data))
}

LN_metastasis <- ifelse(is.na(pN_num), NA, ifelse(pN_num != 0, 1, 0))

Distant_metastasis <- ifelse(
  is.na(pM_num) & is.na(distant_any),
  NA,
  ifelse((!is.na(pM_num) & pM_num != 0) | (!is.na(distant_any) & distant_any), 1, 0)
)

Overall_metastasis <- ifelse(
  is.na(LN_metastasis) & is.na(Distant_metastasis),
  NA,
  ifelse((!is.na(LN_metastasis) & LN_metastasis == 1) |
           (!is.na(Distant_metastasis) & Distant_metastasis == 1), 1, 0)
)

SHCC_data$LN_metastasis <- factor(LN_metastasis, levels = c(0, 1), ordered = TRUE)
SHCC_data$Distant_metastasis <- factor(Distant_metastasis, levels = c(0, 1), ordered = TRUE)
SHCC_data$Overall_metastasis <- factor(Overall_metastasis, levels = c(0, 1), ordered = TRUE)

bin_width <- 0.5
SHCC_data$max_diameter <- to_numeric_safe(SHCC_data$max_diameter)
SHCC_data$max_diameter_round <- round_to_half(SHCC_data$max_diameter)

diameter_df <- SHCC_data %>%
  dplyr::select(max_diameter_round) %>%
  dplyr::filter(!is.na(max_diameter_round))

xmin <- floor(min(diameter_df$max_diameter_round, na.rm = TRUE) / bin_width) * bin_width
xmax <- ceiling(max(diameter_df$max_diameter_round, na.rm = TRUE) / bin_width) * bin_width
breaks_seq <- seq(xmin, xmax + bin_width, by = bin_width)

make_three_panel_plot <- function(data, outcome_var, outcome_title,
                                  bin_width = 0.5,
                                  breaks_seq,
                                  xmin, xmax,
                                  min_bin_n = 10,
                                  spline_df = 4) {
  
  plot_df <- data %>%
    dplyr::select(max_diameter, max_diameter_round, all_of(outcome_var)) %>%
    dplyr::filter(!is.na(max_diameter_round), !is.na(.data[[outcome_var]])) %>%
    dplyr::mutate(
      outcome_num = as.numeric(as.character(.data[[outcome_var]])),
      outcome_label = factor(
        outcome_num,
        levels = c(0, 1),
        labels = c("Non-metastasis", "Metastasis")
      )
    )
  
  
  rate_df <- plot_df %>%
    group_by(max_diameter_round) %>%
    summarise(
      n = n(),
      metastasis_n = sum(outcome_num, na.rm = TRUE),
      metastasis_rate = mean(outcome_num, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    rename(x_mid = max_diameter_round) %>%
    arrange(x_mid)
  
  rate_df_plot <- rate_df %>%
    filter(n >= min_bin_n)
  

  p_hist_all <- ggplot(plot_df, aes(x = max_diameter_round)) +
    geom_histogram(
      binwidth = bin_width,
      boundary = xmin - bin_width / 2,
      closed = "left",
      fill = "#72bcd5",
      color = "black",
      linewidth = 0.2
    ) +
    scale_x_continuous(
      breaks = seq(xmin, xmax, by = 1)
    ) +
    labs(
      title = paste0(outcome_title, ": overall distribution"),
      subtitle = "Diameter rounded to nearest 0.5 cm",
      x = "Tumor max diameter (cm)",
      y = "Count"
    ) +
    theme_bw(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(size = 10),
      panel.grid.minor = element_blank()
    )
  
  p_hist_meta <- ggplot(plot_df, aes(x = max_diameter_round, fill = outcome_label)) +
    geom_histogram(
      binwidth = bin_width,
      boundary = xmin - bin_width / 2,
      closed = "left",
      position = "identity",
      alpha = 0.5,
      color = "black",
      linewidth = 0.2
    ) +
    scale_fill_manual(values = c(
      "Non-metastasis" = "#72bcd5",
      "Metastasis" = "#f7aa58"
    )) +
    scale_x_continuous(
      breaks = seq(xmin, xmax, by = 1)
    ) +
    labs(
      title = paste0(outcome_title, ": distribution by status"),
      subtitle = "Diameter rounded to nearest 0.5 cm",
      x = "Tumor max diameter (cm)",
      y = "Count",
      fill = "Status"
    ) +
    theme_bw(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(size = 10),
      legend.position = "top",
      panel.grid.minor = element_blank()
    )
  
  p_rate <- ggplot(rate_df_plot, aes(x = x_mid, y = metastasis_rate)) +
    geom_point(
      aes(size = n),
      shape = 21,
      fill = "#f7aa58",
      color = "#d95f02",
      stroke = 0.6,
      alpha = 0.85
    ) +
    geom_smooth(
      method = "lm",
      formula = y ~ ns(x, df = spline_df),
      se = TRUE,
      color = "#d95f02",
      fill = "#f7aa58",
      linewidth = 1.2,
      alpha = 0.25
    ) +
    scale_size_area(max_size = 9) +
    scale_x_continuous(
      breaks = seq(xmin, xmax, by = 1)
    ) +
    scale_y_continuous(
      labels = percent_format(accuracy = 1),
      limits = c(0, 1),
      expand = expansion(mult = c(0.02, 0.02))
    ) +
    labs(
      title = paste0(outcome_title, ": crude rate by diameter"),
      subtitle = "Bubble size reflects sample size at each 0.5-cm diameter; cubic spline fit",
      x = "Tumor max diameter (cm)",
      y = "Crude metastasis rate",
      size = "n"
    ) +
    theme_bw(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(size = 10),
      panel.grid.minor = element_blank()
    )
  
  p_combined <- p_hist_all | p_hist_meta | p_rate
  
  list(
    plot_df = plot_df,
    rate_df = rate_df,
    rate_df_plot = rate_df_plot,
    p_hist_all = p_hist_all,
    p_hist_meta = p_hist_meta,
    p_rate = p_rate,
    p_combined = p_combined
  )
}

res_overall <- make_three_panel_plot(
  data = SHCC_data,
  outcome_var = "Overall_metastasis",
  outcome_title = "Overall metastasis",
  bin_width = bin_width,
  breaks_seq = breaks_seq,
  xmin = xmin,
  xmax = xmax,
  min_bin_n =5,
  spline_df = 4
)

res_ln <- make_three_panel_plot(
  data = SHCC_data,
  outcome_var = "LN_metastasis",
  outcome_title = "Lymph node metastasis",
  bin_width = bin_width,
  breaks_seq = breaks_seq,
  xmin = xmin,
  xmax = xmax,
  min_bin_n = 5,
  spline_df = 4
)

res_dm <- make_three_panel_plot(
  data = SHCC_data,
  outcome_var = "Distant_metastasis",
  outcome_title = "Distant metastasis",
  bin_width = bin_width,
  breaks_seq = breaks_seq,
  xmin = xmin,
  xmax = xmax,
  min_bin_n = 5,
  spline_df = 4
)


print(res_overall$p_combined)
print(res_ln$p_combined)
print(res_dm$p_combined)

p_all_combined <- (res_overall$p_hist_all | res_overall$p_hist_meta | res_overall$p_rate) /
  (res_ln$p_hist_all      | res_ln$p_hist_meta      | res_ln$p_rate) /
  (res_dm$p_hist_all      | res_dm$p_hist_meta      | res_dm$p_rate)

print(p_all_combined)

ggsave(
  filename = "outs/max_diameter_three_panel_plot_overall_metastasis.pdf",
  plot = res_overall$p_combined,
  width = 15,
  height = 5.8,
  dpi = 300
)

ggsave(
  filename = "outs/max_diameter_three_panel_plot_ln_metastasis.pdf",
  plot = res_ln$p_combined,
  width = 15,
  height = 5.8,
  dpi = 300
)

ggsave(
  filename = "outs/max_diameter_three_panel_plot_distant_metastasis.pdf",
  plot = res_dm$p_combined,
  width = 15,
  height = 5.8,
  dpi = 300
)

ggsave(
  filename = "outs/max_diameter_nodal_and_distant_metastasis_all_panels.pdf",
  plot = p_all_combined,
  width = 15,
  height = 16.5,
  dpi = 300
)

