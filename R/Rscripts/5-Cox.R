##################################################################################################################
# Cox proportional hazards regression models were used to evaluate the associations between tumor size,          #
# clinicopathological variables, and survival outcomes.                                                          #
##################################################################################################################
rm(list = ls())
gc()

pkg_needed <- c(
  "survival", "dplyr", "broom", "ggplot2", "rms", "Hmisc"
)

pkg_to_install <- pkg_needed[!pkg_needed %in% installed.packages()[, "Package"]]
if (length(pkg_to_install) > 0) install.packages(pkg_to_install)

library(survival)
library(dplyr)
library(broom)
library(ggplot2)
library(rms)
library(Hmisc)

setwd("~/BS_Tumor")
dir.create("outs", showWarnings = FALSE)

input_file <- "data/SHCC_data_merged.RDS"

data <- readRDS(input_file)
data <- as.data.frame(data)

clinical_vars <- c(
  "OSM", "OSstatus",
  "Age",
  "Gender",
  "Location",
  "Pathology_class",
  "Differentiation",
  "Lymphovasscular invasion",
  "Perineural invasion",
  "pTstage",
  "max_diameter"
)

missing_vars <- setdiff(clinical_vars, colnames(data))
if (length(missing_vars) > 0) {
  stop("Missing variables: ", paste(missing_vars, collapse = ", "))
}

df <- data[, clinical_vars, drop = FALSE]

convert_x_to_na <- function(x) {
  if (is.factor(x)) x <- as.character(x)
  if (is.character(x)) {
    x <- trimws(x)
    x[x %in% c("", "x", "X", "NA", "N/A", "Unknown", "unknown")] <- NA
  }
  x
}

df[] <- lapply(df, convert_x_to_na)

names(df)[names(df) == "Lymphovasscular invasion"] <- "LVI"
names(df)[names(df) == "Perineural invasion"] <- "PNI"

numeric_vars <- c("OSM", "OSstatus", "Age", "max_diameter")
for (v in numeric_vars) {
  df[[v]] <- suppressWarnings(as.numeric(as.character(df[[v]])))
}

df <- df %>%
  mutate(
    AgeGroup = case_when(
      Age < 50 ~ "<50",
      Age >= 50 ~ ">=50",
      TRUE ~ NA_character_
    ),
    AgeGroup = factor(AgeGroup, levels = c(">=50", "<50"))
  )

df$Gender <- factor(
  as.character(df$Gender),
  levels = c("0", "1"),
  labels = c("Female", "Male")
)

df$Location <- factor(
  as.character(df$Location),
  levels = c("1", "2", "3"),
  labels = c("Rectum", "Distal_colon", "Proximal_colon")
)

df$Pathology_class <- factor(
  as.character(df$Pathology_class),
  levels = c("1", "2", "3"),
  labels = c("Adenocarcinoma", "Mucinous_adenocarcinoma", "Signet_ring_cell")
)

df$Differentiation <- factor(
  as.character(df$Differentiation),
  levels = c("1", "2", "3"),
  labels = c("Well", "Moderate", "Poor")
)

df$LVI <- factor(
  as.character(df$LVI),
  levels = c("0", "1"),
  labels = c("No", "Yes")
)

df$PNI <- factor(
  as.character(df$PNI),
  levels = c("0", "1"),
  labels = c("No", "Yes")
)

df$pTstage <- factor(
  as.character(df$pTstage),
  levels = c("1", "2", "3", "4"),
  labels = c("T1", "T2", "T3", "T4")
)

analysis_vars <- c(
  "OSM", "OSstatus",
  "AgeGroup",
  "Gender",
  "Location",
  "Pathology_class",
  "Differentiation",
  "LVI",
  "PNI",
  "pTstage",
  "max_diameter"
)

df <- df[, analysis_vars, drop = FALSE] %>%
  filter(
    !is.na(OSM),
    !is.na(OSstatus),
    OSM > 0,
    OSstatus %in% c(0, 1),
    !is.na(max_diameter),
    max_diameter > 0
  ) %>%
  na.omit() %>%
  droplevels()

if (nrow(df) == 0) stop("No samples available after cleaning.")
if (sum(df$OSstatus == 1) == 0) stop("No event=1 samples available.")

cat("Final n:", nrow(df), "\n")
cat("Events:", sum(df$OSstatus == 1), "\n")
cat("Censored:", sum(df$OSstatus == 0), "\n")

write.csv(df, "outs/cox_analysis_dataset_OS.csv", row.names = FALSE)

format_p <- function(p) {
  ifelse(
    is.na(p), NA_character_,
    ifelse(p < 0.001, "<0.001", sprintf("%.3f", p))
  )
}

p_stars <- function(p) {
  ifelse(
    is.na(p), "",
    ifelse(
      p < 0.001, "***",
      ifelse(p < 0.01, "**",
             ifelse(p < 0.05, "*", ""))
    )
  )
}

get_lrt_p <- function(fit_small, fit_large) {
  lrt <- anova(fit_small, fit_large, test = "Chisq")
  lrt_df <- as.data.frame(lrt)
  
  p_col <- grep("Pr\\(|P\\(", colnames(lrt_df), value = TRUE)
  if (length(p_col) == 0) p_col <- grep("^P", colnames(lrt_df), value = TRUE)
  
  if (length(p_col) == 0) return(NA_real_)
  
  suppressWarnings(as.numeric(lrt_df[nrow(lrt_df), p_col[1]]))
}

make_subgroup_n <- function(data, var) {
  tab <- table(data[[var]])
  paste0(names(tab), "=", as.integer(tab), collapse = "; ")
}

uni_vars <- c(
  "AgeGroup",
  "Gender",
  "Location",
  "Pathology_class",
  "Differentiation",
  "LVI",
  "PNI",
  "pTstage",
  "max_diameter"
)

run_uni_cox <- function(var, data) {
  fm <- as.formula(paste0("Surv(OSM, OSstatus) ~ ", var))
  fit <- coxph(fm, data = data)
  
  broom::tidy(fit, exponentiate = TRUE, conf.int = TRUE) %>%
    mutate(
      Variable = var,
      N = nrow(data),
      Events = sum(data$OSstatus == 1)
    ) %>%
    select(Variable, term, N, Events, estimate, conf.low, conf.high, p.value) %>%
    rename(
      HR = estimate,
      CI_lower = conf.low,
      CI_upper = conf.high,
      P_value = p.value
    )
}

uni_results <- lapply(uni_vars, function(v) {
  tryCatch(
    run_uni_cox(v, df),
    error = function(e) {
      data.frame(
        Variable = v,
        term = NA_character_,
        N = nrow(df),
        Events = sum(df$OSstatus == 1),
        HR = NA_real_,
        CI_lower = NA_real_,
        CI_upper = NA_real_,
        P_value = NA_real_
      )
    }
  )
}) %>%
  bind_rows() %>%
  mutate(
    HR_95CI = ifelse(
      is.na(HR), NA_character_,
      sprintf("%.3f (%.3f-%.3f)", HR, CI_lower, CI_upper)
    ),
    P_value_format = paste0(format_p(P_value), p_stars(P_value))
  )

write.csv(
  uni_results,
  "outs/cox_univariate_OS_results.csv",
  row.names = FALSE
)

multi_formula <- as.formula(
  "Surv(OSM, OSstatus) ~
   AgeGroup + Gender + Location + Pathology_class +
   Differentiation + LVI + PNI + pTstage + max_diameter"
)

multi_fit <- coxph(multi_formula, data = df)

multi_results <- broom::tidy(
  multi_fit,
  exponentiate = TRUE,
  conf.int = TRUE
) %>%
  mutate(
    N = nrow(df),
    Events = sum(df$OSstatus == 1)
  ) %>%
  select(term, N, Events, estimate, conf.low, conf.high, p.value) %>%
  rename(
    HR = estimate,
    CI_lower = conf.low,
    CI_upper = conf.high,
    P_value = p.value
  ) %>%
  mutate(
    HR_95CI = sprintf("%.3f (%.3f-%.3f)", HR, CI_lower, CI_upper),
    P_value_format = paste0(format_p(P_value), p_stars(P_value))
  )

write.csv(
  multi_results,
  "outs/cox_multivariable_OS_results.csv",
  row.names = FALSE
)

sink("outs/cox_multivariable_OS_summary.txt")
print(summary(multi_fit))
sink()

ph_test <- cox.zph(multi_fit)

capture.output(
  ph_test,
  file = "outs/cox_multivariable_OS_PH_assumption_test.txt"
)

pdf("outs/cox_multivariable_OS_PH_assumption_plots.pdf", width = 10, height = 8)
plot(ph_test)
dev.off()

dd <- datadist(df)
options(datadist = "dd")

rcs_fit <- cph(
  Surv(OSM, OSstatus) ~
    AgeGroup + Gender + Location + Pathology_class +
    Differentiation + LVI + PNI + pTstage +
    rcs(max_diameter, 4),
  data = df,
  x = TRUE,
  y = TRUE,
  surv = TRUE
)

sink("outs/rcs_cox_max_diameter_OS_summary.txt")
print(summary(rcs_fit))
sink()

anova_rcs <- anova(rcs_fit)

capture.output(
  anova_rcs,
  file = "outs/rcs_cox_max_diameter_OS_anova.txt"
)

anova_rcs_df <- as.data.frame(anova_rcs)
anova_rcs_df$Effect <- rownames(anova_rcs_df)

write.csv(
  anova_rcs_df,
  "outs/rcs_cox_max_diameter_OS_anova_table.csv",
  row.names = FALSE
)

# 提取 max_diameter overall 和 nonlinear P
p_col <- grep("^P", colnames(anova_rcs_df), value = TRUE)[1]

p_overall_rcs <- NA_real_
p_nonlinear_rcs <- NA_real_

if (!is.na(p_col)) {
  idx_overall <- which(anova_rcs_df$Effect == "max_diameter")
  idx_nonlinear <- grep("Nonlinear", anova_rcs_df$Effect)
  
  if (length(idx_overall) > 0) {
    p_overall_rcs <- as.numeric(anova_rcs_df[idx_overall[1], p_col])
  }
  if (length(idx_nonlinear) > 0) {
    p_nonlinear_rcs <- as.numeric(anova_rcs_df[idx_nonlinear[1], p_col])
  }
}

pred_rcs <- Predict(
  rcs_fit,
  max_diameter,
  fun = exp,
  ref.zero = TRUE
)

p_rcs <- ggplot(pred_rcs, aes(x = max_diameter, y = yhat)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.20) +
  geom_line(linewidth = 1.1) +
  geom_hline(yintercept = 1, linetype = 2, color = "gray50") +
  labs(
    title = "Restricted cubic spline Cox model for max_diameter",
    subtitle = paste0(
      "Outcome: OS; ",
      "P overall = ", format_p(p_overall_rcs), p_stars(p_overall_rcs),
      "; P nonlinear = ", format_p(p_nonlinear_rcs), p_stars(p_nonlinear_rcs)
    ),
    x = "max_diameter",
    y = "Hazard ratio"
  ) +
  theme_bw(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(size = 11)
  )

ggsave(
  "outs/rcs_cox_max_diameter_OS.pdf",
  p_rcs,
  width = 6,
  height = 5
)

rcs_basis <- Hmisc::rcspline.eval(
  df$max_diameter,
  nk = 4,
  inclx = TRUE
)

rcs_knots <- attr(rcs_basis, "knots")

rcs_basis <- as.data.frame(rcs_basis)
rcs_terms <- paste0("md_rcs", seq_len(ncol(rcs_basis)))
names(rcs_basis) <- rcs_terms

df_int <- cbind(df, rcs_basis)

add_rcs_basis <- function(newdata, knots, rcs_terms) {
  basis <- Hmisc::rcspline.eval(
    newdata$max_diameter,
    knots = knots,
    inclx = TRUE
  )
  
  basis <- as.data.frame(basis)
  names(basis) <- rcs_terms
  
  cbind(newdata, basis)
}

interaction_vars <- c(
  "AgeGroup",
  "Gender",
  "Location",
  "Differentiation",
  "Pathology_class",
  "pTstage",
  "LVI",
  "PNI"
)

adjust_vars <- c(
  "AgeGroup",
  "Gender",
  "Location",
  "Pathology_class",
  "Differentiation",
  "LVI",
  "PNI",
  "pTstage"
)

get_interaction_result <- function(var, data) {
  other_vars <- setdiff(adjust_vars, var)
  
  base_formula <- as.formula(
    paste(
      "Surv(OSM, OSstatus) ~",
      paste(c(other_vars, var, rcs_terms), collapse = " + ")
    )
  )
  
  linear_int_formula <- as.formula(
    paste(
      "Surv(OSM, OSstatus) ~",
      paste(c(
        other_vars,
        var,
        rcs_terms,
        paste0(rcs_terms[1], ":", var)
      ), collapse = " + ")
    )
  )
  
  full_int_formula <- as.formula(
    paste(
      "Surv(OSM, OSstatus) ~",
      paste(c(
        other_vars,
        var,
        rcs_terms,
        paste0(rcs_terms, ":", var)
      ), collapse = " + ")
    )
  )
  
  fit_base <- try(coxph(base_formula, data = data, x = TRUE), silent = TRUE)
  fit_linear <- try(coxph(linear_int_formula, data = data, x = TRUE), silent = TRUE)
  fit_full <- try(coxph(full_int_formula, data = data, x = TRUE), silent = TRUE)
  
  if (
    inherits(fit_base, "try-error") ||
    inherits(fit_linear, "try-error") ||
    inherits(fit_full, "try-error")
  ) {
    return(data.frame(
      Interaction_term = paste0("max_diameter × ", var),
      Variable = var,
      Subgroup_N = make_subgroup_n(data, var),
      N = nrow(data),
      Events = sum(data$OSstatus == 1),
      P_for_interaction_overall = NA_real_,
      P_for_interaction_nonlinear = NA_real_
    ))
  }
  
  p_overall <- get_lrt_p(fit_base, fit_full)
  p_nonlinear <- get_lrt_p(fit_linear, fit_full)
  
  data.frame(
    Interaction_term = paste0("max_diameter × ", var),
    Variable = var,
    Subgroup_N = make_subgroup_n(data, var),
    N = nrow(data),
    Events = sum(data$OSstatus == 1),
    P_for_interaction_overall = p_overall,
    P_for_interaction_nonlinear = p_nonlinear
  )
}

interaction_results <- lapply(
  interaction_vars,
  get_interaction_result,
  data = df_int
) %>%
  bind_rows() %>%
  mutate(
    P_for_interaction_overall_format =
      paste0(format_p(P_for_interaction_overall), p_stars(P_for_interaction_overall)),
    P_for_interaction_nonlinear_format =
      paste0(format_p(P_for_interaction_nonlinear), p_stars(P_for_interaction_nonlinear))
  )

write.csv(
  interaction_results,
  "outs/cox_RCS_interaction_OS_results.csv",
  row.names = FALSE
)

make_interaction_plot <- function(var, data, knots, rcs_terms, interaction_results) {
  other_vars <- setdiff(adjust_vars, var)
  
  full_int_formula <- as.formula(
    paste(
      "Surv(OSM, OSstatus) ~",
      paste(c(
        other_vars,
        var,
        rcs_terms,
        paste0(rcs_terms, ":", var)
      ), collapse = " + ")
    )
  )
  
  fit <- coxph(full_int_formula, data = data, x = TRUE)
  
  x_seq <- seq(
    quantile(data$max_diameter, 0.01, na.rm = TRUE),
    quantile(data$max_diameter, 0.99, na.rm = TRUE),
    length.out = 120
  )
  
  ref_x <- median(data$max_diameter, na.rm = TRUE)
  
  newdata <- do.call(
    rbind,
    lapply(levels(data[[var]]), function(lv) {
      tmp <- data.frame(max_diameter = x_seq)
      
      for (v in adjust_vars) {
        tmp[[v]] <- factor(
          levels(data[[v]])[1],
          levels = levels(data[[v]])
        )
      }
      
      tmp[[var]] <- factor(lv, levels = levels(data[[var]]))
      tmp
    })
  )
  
  refdata <- newdata
  refdata$max_diameter <- ref_x
  
  newdata <- add_rcs_basis(newdata, knots, rcs_terms)
  refdata <- add_rcs_basis(refdata, knots, rcs_terms)
  
  tt <- delete.response(terms(fit))
  
  X <- model.matrix(tt, newdata, contrasts.arg = fit$contrasts)
  X_ref <- model.matrix(tt, refdata, contrasts.arg = fit$contrasts)
  
  beta <- coef(fit)
  keep <- !is.na(beta)
  
  beta <- beta[keep]
  V <- vcov(fit)[keep, keep, drop = FALSE]
  
  X <- X[, names(beta), drop = FALSE]
  X_ref <- X_ref[, names(beta), drop = FALSE]
  
  X_diff <- X - X_ref
  
  eta <- as.vector(X_diff %*% beta)
  se <- sqrt(rowSums((X_diff %*% V) * X_diff))
  
  plot_df <- newdata %>%
    mutate(
      HR = exp(eta),
      lower = exp(eta - 1.96 * se),
      upper = exp(eta + 1.96 * se),
      subgroup = .data[[var]]
    )
  
  p_row <- interaction_results %>%
    filter(Variable == var)
  
  subtitle_text <- paste0(
    "Reference max_diameter = ", round(ref_x, 2),
    "; P overall = ",
    p_row$P_for_interaction_overall_format[1],
    "; P nonlinear = ",
    p_row$P_for_interaction_nonlinear_format[1]
  )
  
  p <- ggplot(
    plot_df,
    aes(x = max_diameter, y = HR, color = subgroup, fill = subgroup)
  ) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.12, color = NA) +
    geom_line(linewidth = 1.0) +
    geom_hline(yintercept = 1, linetype = 2, color = "gray50") +
    scale_y_log10() +
    labs(
      title = paste0("Interaction plot: max_diameter × ", var),
      subtitle = subtitle_text,
      x = "max_diameter",
      y = "Hazard ratio for OS",
      color = var,
      fill = var
    ) +
    theme_bw(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(size = 10),
      legend.position = "right"
    )
  
  ggsave(
    filename = paste0("outs/interaction_plot_max_diameter_by_", var, "_OS.pdf"),
    plot = p,
    width = 7.5,
    height = 5.5
  )
  
  p
}

interaction_plots <- lapply(
  interaction_vars,
  make_interaction_plot,
  data = df_int,
  knots = rcs_knots,
  rcs_terms = rcs_terms,
  interaction_results = interaction_results
)

names(interaction_plots) <- interaction_vars

pdf("outs/interaction_plots_max_diameter_OS_all.pdf", width = 7.5, height = 5.5)
for (p in interaction_plots) print(p)
dev.off()

saveRDS(df, "outs/cox_analysis_dataset_OS.RDS")
saveRDS(multi_fit, "outs/cox_multivariable_OS_model.rds")
saveRDS(rcs_fit, "outs/rcs_cox_max_diameter_OS_model.rds")
saveRDS(interaction_results, "outs/cox_RCS_interaction_OS_results.rds")

format_p2 <- function(p) {
  ifelse(
    is.na(p), "",
    ifelse(p < 0.001, "<0.001", sprintf("%.3f", p))
  )
}

p_stars2 <- function(p) {
  ifelse(
    is.na(p), "",
    ifelse(
      p < 0.001, "***",
      ifelse(
        p < 0.01, "**",
        ifelse(p < 0.05, "*", "")
      )
    )
  )
}

format_p_with_star <- function(p) {
  paste0(format_p2(p), p_stars2(p))
}

format_hr <- function(hr, low, high) {
  ifelse(
    is.na(hr), "",
    sprintf("%.2f (%.2f,%.2f)", hr, low, high)
  )
}

anova_rcs_df <- as.data.frame(anova_rcs)
anova_rcs_df$Effect <- rownames(anova_rcs_df)

p_col <- grep("^P", colnames(anova_rcs_df), value = TRUE)[1]

get_anova_p <- function(var) {
  idx <- which(anova_rcs_df$Effect == var)
  if (length(idx) == 0 || is.na(p_col)) return(NA_real_)
  as.numeric(anova_rcs_df[idx[1], p_col])
}

get_nonlinear_p <- function(var) {
  idx_var <- which(anova_rcs_df$Effect == var)
  idx_non <- grep("Nonlinear", anova_rcs_df$Effect)
  
  if (length(idx_var) == 0 || length(idx_non) == 0 || is.na(p_col)) {
    return(NA_real_)
  }
  
  idx_non_after <- idx_non[idx_non > idx_var[1]]
  if (length(idx_non_after) == 0) return(NA_real_)
  
  as.numeric(anova_rcs_df[idx_non_after[1], p_col])
}


coef_vec <- coef(rcs_fit)
vc <- vcov(rcs_fit)
se_vec <- sqrt(diag(vc))

coef_table <- data.frame(
  term = names(coef_vec),
  beta = as.numeric(coef_vec),
  SE = as.numeric(se_vec),
  row.names = NULL
) %>%
  mutate(
    HR = exp(beta),
    CI_lower = exp(beta - 1.96 * SE),
    CI_upper = exp(beta + 1.96 * SE)
  )

get_coef_row <- function(var, level) {
  target1 <- paste0(var, "=", level)
  target2 <- paste0(var, level)
  
  idx <- which(coef_table$term == target1)
  if (length(idx) == 0) idx <- which(coef_table$term == target2)
  
  if (length(idx) == 0) {
    idx <- grep(
      paste0("^", var, ".*", level),
      coef_table$term
    )
  }
  
  if (length(idx) == 0) {
    return(data.frame(
      HR = NA_real_,
      CI_lower = NA_real_,
      CI_upper = NA_real_
    ))
  }
  
  coef_table[idx[1], c("HR", "CI_lower", "CI_upper")]
}

var_label <- c(
  max_diameter = "max_diameter",
  AgeGroup = "Age",
  Gender = "Gender",
  Location = "Location",
  Differentiation = "Differentiation",
  Pathology_class = "Pathology Class",
  pTstage = "pTstage",
  LVI = "LVI",
  PNI = "PNI"
)

level_label <- list(
  AgeGroup = c(
    "<50" = "<50",
    ">=50" = "≥50"
  ),
  Gender = c(
    "Female" = "Female",
    "Male" = "Male"
  ),
  Location = c(
    "Rectum" = "Rectum",
    "Distal_colon" = "Distal",
    "Proximal_colon" = "Proximal"
  ),
  Differentiation = c(
    "Well" = "Well",
    "Moderate" = "Medium",
    "Poor" = "Poor"
  ),
  Pathology_class = c(
    "Adenocarcinoma" = "Class I",
    "Mucinous_adenocarcinoma" = "Class II",
    "Signet_ring_cell" = "Class III"
  ),
  pTstage = c(
    "T1" = "T1",
    "T2" = "T2",
    "T3" = "T3",
    "T4" = "T4"
  ),
  LVI = c(
    "No" = "No",
    "Yes" = "Yes"
  ),
  PNI = c(
    "No" = "No",
    "Yes" = "Yes"
  )
)

make_factor_rows <- function(var, data, compact_binary = FALSE) {
  
  lvls <- levels(data[[var]])
  ref <- lvls[1]
  p_overall <- get_anova_p(var)
  
  # 二分类变量压缩为一行，例如 LVI、PNI
  if (compact_binary && length(lvls) == 2) {
    target_level <- lvls[2]
    hr_row <- get_coef_row(var, target_level)
    
    return(data.frame(
      Variable = var_label[[var]],
      N = as.character(sum(data[[var]] == target_level, na.rm = TRUE)),
      `HR (95%CI)` = format_hr(
        hr_row$HR,
        hr_row$CI_lower,
        hr_row$CI_upper
      ),
      `P value` = format_p_with_star(p_overall),
      `P for non-linear` = "",
      check.names = FALSE
    ))
  }
  
  out <- data.frame(
    Variable = var_label[[var]],
    N = "",
    `HR (95%CI)` = "",
    `P value` = format_p_with_star(p_overall),
    `P for non-linear` = "",
    check.names = FALSE
  )
  
  for (lv in lvls) {
    
    show_lv <- level_label[[var]][[lv]]
    if (is.null(show_lv)) show_lv <- lv
    
    if (lv == ref) {
      hr_text <- "Ref"
    } else {
      hr_row <- get_coef_row(var, lv)
      hr_text <- format_hr(
        hr_row$HR,
        hr_row$CI_lower,
        hr_row$CI_upper
      )
    }
    
    tmp <- data.frame(
      Variable = paste0("    ", show_lv),
      N = as.character(sum(data[[var]] == lv, na.rm = TRUE)),
      `HR (95%CI)` = hr_text,
      `P value` = "",
      `P for non-linear` = "",
      check.names = FALSE
    )
    
    out <- bind_rows(out, tmp)
  }
  
  out
}

rcs_main_row <- data.frame(
  Variable = "max_diameter",
  N = as.character(nrow(df)),
  `HR (95%CI)` = "NA",
  `P value` = format_p_with_star(get_anova_p("max_diameter")),
  `P for non-linear` = format_p_with_star(get_nonlinear_p("max_diameter")),
  check.names = FALSE
)

rcs_pub_table <- bind_rows(
  rcs_main_row,
  make_factor_rows("AgeGroup", df),
  make_factor_rows("Gender", df),
  make_factor_rows("Location", df),
  make_factor_rows("Differentiation", df),
  make_factor_rows("Pathology_class", df),
  make_factor_rows("pTstage", df),
  make_factor_rows("LVI", df, compact_binary = TRUE),
  make_factor_rows("PNI", df, compact_binary = TRUE)
)

write.csv(
  rcs_pub_table,
  "outs/rcs_cox_multivariable_OS_publication_table.csv",
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

rcs_pub_table