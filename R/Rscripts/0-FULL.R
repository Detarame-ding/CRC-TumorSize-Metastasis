##################################################################################################################
# The dataset was preprocessed by cleaning variables, recoding categorical features, handling missing values,    #
# and preparing analysis-ready variables for subsequent statistical and modeling analyses.                       #
##################################################################################################################

rm(list = ls())
library(readxl)
library(survival)
library(survminer)
library(ggplot2)
setwd('~/BS_Tumor')

SHCC_data <- read_excel("data/SHCC_data.xlsx")
SHCC_data <- SHCC_data[!is.na(SHCC_data$ID), ]
SHCC_data$Gender <- SHCC_data$`Gender(0女1男)`
SHCC_data$Location_minor <- factor(SHCC_data$`部位0未指明的 1回盲部2阑尾3升结肠（右半结肠）4肝曲5横结肠6脾曲7降结肠（左半结肠）8乙状结肠9直乙结肠10直肠11肛管`)
SHCC_data$Location <- factor(SHCC_data$`Location（0未制定1直肠2远端3近端4多发）`)
SHCC_data$max_diameter <- SHCC_data$`肿瘤最大径（cm）...33`
SHCC_data$Pathology_class <- factor(SHCC_data$`病理分类（0未鉴定1腺癌2腺癌伴黏液3腺癌伴印戒细胞4高级别局灶癌变）`)
SHCC_data$LymphNodePositive <- factor(SHCC_data$`LN检出数量`)
SHCC_data$LymphNodeTotal <- factor(SHCC_data$`LN检查总数`)
SHCC_data$recurrence <- factor(SHCC_data$`术后转移复发`)
SHCC_data$KRAS_result_L1 <- factor(SHCC_data$`KRAS-1级(0\tWT\r\n1\tMT/可疑突变\r\n2\t未检测\r\n3\t检测失败/无法判断)`)
SHCC_data$NRAS_result_L1 <- factor(SHCC_data$`NRAS-1级(0\tWT\r\n1\tMT/可疑突变\r\n2\t未检测\r\n3\t检测失败/无法判断)`)
SHCC_data$BRAF_result_L1 <- factor(SHCC_data$`BRAF-1级(0\tWT\r\n1\tMT/可疑突变\r\n2\t未检测\r\n3\t检测失败/无法判断)`)
SHCC_data$KRAS_result_L2 <- factor(SHCC_data$`KRAS-2级(0\tWT\r\n112\t12号外显子突变\r\n114\t14号外显子突变\r\n12\t2号外显子突变\r\n123\t2、3号都突变\r\n124\t2、4号都突变\r\n13\t3号突变\r\n134\t3、4号都突变\r\n14\t4号突变)`)
SHCC_data$NRAS_result_L2 <- factor(SHCC_data$`NRAS-2级(NRAS结果-2级\t\r\n0\tWT\r\n13\t3号外显子突变\r\n12\t2号外显子突变\r\n14\t4号突变)`)
SHCC_data$BRAF_result_L2 <- factor(SHCC_data$`BRAF-2级(0\tWT\r\n115\t15号外显子突变\r\n114\t14号外显子突变\r\n14\t4号突变\r\n117\t17号外显子突变\r\n118\t18号外显子突变\r\n1610\t6/10都突变)`)
SHCC_data$MSI_result_L1 <- factor(SHCC_data$`MSI-1级(0\tMSS\r\n1\tMSI\r\n2\t未检测\r\n3\t检测失败/无法判断)`)
SHCC_data$MSI_result_L2 <- factor(SHCC_data$`MSI-2级(0\tMSS\r\n11\tMSI-H\r\n12\tMSI-L)`)

# Metastasis was defined as the presence of either lymph node metastasis or distant metastasis.
SHCC_data$Metastasis <- as.numeric((SHCC_data$pNstage != 0) | (SHCC_data$pMstage != 0) | (SHCC_data$livermeta != 0) | (SHCC_data$bonemeta != 0) | (SHCC_data$lungmeta != 0) | (SHCC_data$brainmeta != 0))

cutoff <- 5
SHCC_data$Tsize <- ifelse(SHCC_data$max_diameter < cutoff, "Small", "Large")
SHCC_data$Ttype <- ifelse(SHCC_data$Tsize == 'Small',
                          ifelse(SHCC_data$Metastasis == 0, 'S+NM', 'S+M'),
                          ifelse(SHCC_data$Metastasis == 0, 'L+NM', 'L+M'))

if(TRUE){
  features <- c('ID','pathID','Year','Age','Gender','EMR','NAC', 'Location_minor', 'Location','Pathology_class','Mucinous subtype',
                'Differentiation','Lymphovasscular invasion','Perineural invasion','CRM','TD','LymphNodePositive','LymphNodeTotal',
                'pTstage','pNstage','pMstage','TNMstage','max_diameter','livermeta','lungmeta','bonemeta','brainmeta','recurrence',
                'TTPstatus','TTPM','OSstatus','OSM','DFSstatus','DFSM','ki67','PCNA','P21-1','P53-1','CD44-1','nm23','COX2','EGFR',
                'KRAS_result_L1','NRAS_result_L1','BRAF_result_L1','KRAS_result_L2','NRAS_result_L2','BRAF_result_L2','MLH1','MSH2',
                'MSH6','PMS2','MSI_result_L1','MSI_result_L2','MMRstatus','MMRtest','MMRMSIstatus','MMRMSItest','Metastasis','Tsize','Ttype')
  
  SHCC_data <- SHCC_data[, features, drop = FALSE]
  
  char_cols <- c("ID", "pathID")
  num_cols <- c(
    "Age",  "Year", "LymphNodePositive", "LymphNodeTotal", "max_diameter",
    "TTPM", "OSM", "DFSM",
    "ki67", "PCNA"
  )
  
  factor_cols <- c(
    "Gender", "EMR", "NAC",
    "Location", "Pathology_class", "Mucinous subtype",
    "Lymphovasscular invasion", "Perineural invasion", "CRM", "TD",
    "livermeta", "lungmeta", "bonemeta", "brainmeta", "recurrence",
    "TTPstatus", "OSstatus", "DFSstatus",
    "KRAS_result_L1", "NRAS_result_L1", "BRAF_result_L1",
    "KRAS_result_L2", "NRAS_result_L2", "BRAF_result_L2",
    "MLH1", "MSH2", "MSH6", "PMS2",
    "MMRstatus", "MMRtest", "MMRMSIstatus", "MMRMSItest",
    "Metastasis", "Tsize", "Ttype", "P21-1", "P53-1", "CD44-1", "nm23", "COX2", "EGFR"
  )
  
  ord_cols <- c("Differentiation", "pTstage", "pNstage", "pMstage", "TNMstage", "Location_minor")
  
  SHCC_data[char_cols] <- lapply(SHCC_data[char_cols], as.character)
  SHCC_data[num_cols] <- lapply(
    SHCC_data[num_cols],
    function(x) as.numeric(as.character(x))
  )
  SHCC_data[factor_cols] <- lapply(
    SHCC_data[factor_cols],
    function(x) {
      x <- as.character(x)
      x[is.na(x)] <- "x"
      as.factor(x)
    }
  )
  clean_ordered <- function(x, levels_vec) {
    x <- as.character(x)
    factor(x, levels = levels_vec, ordered = TRUE)
  }
  SHCC_data$Differentiation <- clean_ordered(SHCC_data$Differentiation, c("1", "2", "3")) 
  SHCC_data$pTstage <- clean_ordered(SHCC_data$pTstage, c("0", "1", "2", "3", "4"))
  SHCC_data$pNstage <- clean_ordered(SHCC_data$pNstage, c("0", "1", "2"))
  SHCC_data$pMstage <- clean_ordered(SHCC_data$pMstage, c("0", "1"))
  SHCC_data$TNMstage <- clean_ordered(SHCC_data$TNMstage, c("0", "1", "2", "3", "4"))
  SHCC_data$Location_minor <- clean_ordered(SHCC_data$Location_minor, c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11"))
  SHCC_data <- SHCC_data[complete.cases(SHCC_data[, ord_cols, drop = FALSE]), ]
}

saveRDS(SHCC_data, 'data/SHCC_data_cleaned.RDS')
write_xlsx(SHCC_data, "data/SHCC_data_cleaned.xlsx")
write.csv(SHCC_data, "data/SHCC_data_cleaned.csv", row.names = FALSE)

##################################################################################################################
# Descriptive statistical tables were generated to summarize baseline clinicopathological characteristics        #
# and compare variables across predefined patient groups.                                                        #
##################################################################################################################
rm(list = ls())
library(readxl)
library(survival)
library(survminer)
library(ggplot2)
setwd('~/BS_Tumor')

SHCC_data <- readRDS("data/SHCC_data_cleaned.RDS")
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
convert_if_needed <- function(x) {
  if (is.character(x)) as.factor(x) else x
}
is_binary_var <- function(x) {
  ux <- unique(x[!is.na(x)])
  length(ux) == 2
}
is_multiclass_var <- function(x) {
  ux <- unique(x[!is.na(x)])
  length(ux) > 2
}
fmt_mean_sd <- function(x) {
  if (all(is.na(x))) return("NA")
  sprintf("%.2f ± %.2f", mean(x, na.rm = TRUE), sd(x, na.rm = TRUE))
}
get_cat_pvalue <- function(x, group) {
  tab <- table(x, group)
  if (nrow(tab) < 2 || ncol(tab) < 2) return(NA)
  
  suppressWarnings({
    chisq_res <- try(chisq.test(tab), silent = TRUE)
  })
  
  if (!inherits(chisq_res, "try-error")) {
    expected <- chisq_res$expected
    if (all(expected >= 5)) {
      return(chisq_res$p.value)
    } else {
      fisher_res <- try(fisher.test(tab), silent = TRUE)
      if (!inherits(fisher_res, "try-error")) {
        return(fisher_res$p.value)
      } else {
        return(chisq_res$p.value)
      }
    }
  } else {
    fisher_res <- try(fisher.test(tab), silent = TRUE)
    if (!inherits(fisher_res, "try-error")) {
      return(fisher_res$p.value)
    } else {
      return(NA)
    }
  }
}
get_num_pvalue <- function(x, group) {
  dtmp <- data.frame(x = x, group = group)
  dtmp <- dtmp[!is.na(dtmp$x) & !is.na(dtmp$group), ]
  
  if (length(unique(dtmp$group)) != 2) return(NA)
  if (nrow(dtmp) == 0) return(NA)
  res <- try(t.test(x ~ group, data = dtmp), silent = TRUE)
  if (!inherits(res, "try-error")) {
    return(res$p.value)
  } else {
    return(NA)
  }
}
fmt_n_pct_level <- function(x, level) {
  x2 <- x[!is.na(x)]
  if (length(x2) == 0) return("NA")
  n <- sum(x2 == level, na.rm = TRUE)
  pct <- n / length(x2) * 100
  sprintf("%d (%.1f%%)", n, pct)
}
result_list <- list()
for (v in vars_to_analyze) {
  x <- df[[v]]
  x <- convert_if_needed(x)
  g <- df$TumorSizeGroup
  
  if (is.numeric(x)) {
    small_x <- x[g == "Small"]
    large_x <- x[g == "Large"]
    
    pval <- get_num_pvalue(x, g)
    
    one_row <- data.frame(
      Variable = v,
      Level = "",
      Type = "Numeric",
      Small = fmt_mean_sd(small_x),
      Large = fmt_mean_sd(large_x),
      P_value = pval,
      stringsAsFactors = FALSE
    )
    result_list[[length(result_list) + 1]] <- one_row
  }
  else {
    x <- as.factor(x)
    lvls <- levels(x)
    pval <- get_cat_pvalue(x, g)
    
    cat_rows <- lapply(seq_along(lvls), function(i) {
      lev <- lvls[i]
      
      small_x <- x[g == "Small"]
      large_x <- x[g == "Large"]
      
      data.frame(
        Variable = v,
        Level = as.character(lev),
        Type = ifelse(length(lvls) == 2, "Binary", "Multiclass"),
        Small = fmt_n_pct_level(small_x, lev),
        Large = fmt_n_pct_level(large_x, lev),
        P_value = ifelse(i == 1, pval, NA),
        stringsAsFactors = FALSE
      )
    })
    
    result_list[[length(result_list) + 1]] <- do.call(rbind, cat_rows)
  }
}
result_table <- do.call(rbind, result_list)
result_table$P_value_fmt <- ifelse(
  is.na(result_table$P_value),
  "",
  ifelse(result_table$P_value < 0.001, "<0.001", sprintf("%.3f", result_table$P_value))
)
print(result_table, row.names = FALSE)
write.csv(result_table, "TumorSizeGroup_descriptive_statistics.csv", row.names = FALSE)

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

##################################################################################################################
# SEER-based validation was performed to assess the external reproducibility of the associations between tumor   #
#size, metastasis risk, and survival outcomes.                                                                   #
##################################################################################################################
rm(list = ls())
library(readr)
library(dplyr)
library(stringr)
txt_file <- "data/seer.txt"
dic_file <- "data/seer.dic"

dic_lines <- readLines(dic_file, warn = FALSE, encoding = "UTF-8")
var_lines <- dic_lines[grepl("^Var[0-9]+Name=", dic_lines)]
var_info <- data.frame(
  var_order = as.integer(str_match(var_lines, "^Var([0-9]+)Name=")[, 2]),
  var_name  = sub("^Var[0-9]+Name=", "", var_lines),
  stringsAsFactors = FALSE
) %>%
  arrange(var_order)

col_names <- var_info$var_name
col_names <- make.unique(col_names)
print(col_names)

seer_data <- read_delim(
  file = txt_file,
  delim = "\t",
  col_names = col_names,
  col_types = cols(.default = col_character()),
  na = c("", " ", "NA", "N/A", "Unknown"),
  trim_ws = TRUE,
  quote = "\"",
  locale = locale(encoding = "UTF-8")
)

dim(seer_data)
head(seer_data)
colnames(seer_data)
saveRDS(seer_data, "data/seer.rds")
seer_data <- readRDS('data/seer.rds')

table(seer_data$`Age recode with <1 year olds and 90+`)
table(seer_data$`Year of diagnosis`)
table(seer_data$`CS tumor size (2004-2015)`)

clean_blank <- function(x) {
  x <- as.character(x)
  x <- str_squish(x)
  x[x %in% c("", " ", "Blank(s)", "Unknown", "Not available", "Not applicable")] <- NA
  return(x)
}

parse_tumor_size_mm <- function(x) {
  x <- clean_blank(x)
  val <- suppressWarnings(as.numeric(str_extract(x, "\\d+")))
  val[val %in% c(0, 990:999)] <- NA_real_
  return(val)
}

extract_T <- function(x) {
  x <- str_to_upper(clean_blank(x))
  str_extract(x, "TIS|TX|T0|T[1-4][A-D]?")
}

extract_N <- function(x) {
  x <- str_to_upper(clean_blank(x))
  str_extract(x, "NX|N0|N[1-3][A-C]?")
}

extract_M <- function(x) {
  x <- str_to_upper(clean_blank(x))
  str_extract(x, "MX|M0|M1[A-C]?|M1")
}

recode_yes_no <- function(x) {
  y <- str_to_lower(clean_blank(x))
  
  case_when(
    is.na(y) ~ NA_character_,
    str_detect(y, "unknown|not documented|not applicable|blank") ~ NA_character_,
    str_detect(y, "not identified|not present|negative|absent|none|^0$|no") ~ "Negative",
    str_detect(y, "identified|present|positive|^1$|yes") ~ "Positive",
    TRUE ~ NA_character_
  )
}

seer_clean <- seer_data %>%
  mutate(across(everything(), clean_blank)) %>%
  mutate(
    Year = as.integer(`Year of diagnosis`),
    tumor_size_cs_mm = parse_tumor_size_mm(`CS tumor size (2004-2015)`),
    tumor_size_summary_mm = parse_tumor_size_mm(`Tumor Size Summary (2016+)`),
    tumor_size_mm = coalesce(
      tumor_size_cs_mm,
      tumor_size_summary_mm
    ),
    
    tumor_size_cm = tumor_size_mm / 10,
    
    tumor_size_source = case_when(
      !is.na(tumor_size_cs_mm) ~ "CS tumor size (2004-2015)",
      !is.na(tumor_size_summary_mm) ~ "Tumor Size Summary (2016+)",
      TRUE ~ NA_character_
    ),
    T_raw = coalesce(
      `T value - based on AJCC 3rd (1988-2003)`,
      `Derived AJCC T, 6th ed (2004-2015)`,
      `Derived SEER Combined T (2016-2017)`,
      `Derived EOD 2018 T Recode (2018+)`
    ),
    T_stage = extract_T(T_raw),
    T_major = case_when(
      str_detect(T_stage, "^TIS$") ~ "Tis",
      str_detect(T_stage, "^T0$")  ~ "T0",
      str_detect(T_stage, "^T1")   ~ "T1",
      str_detect(T_stage, "^T2")   ~ "T2",
      str_detect(T_stage, "^T3")   ~ "T3",
      str_detect(T_stage, "^T4")   ~ "T4",
      str_detect(T_stage, "^TX$")  ~ NA_character_,
      TRUE ~ NA_character_
    ),
    T_source = case_when(
      !is.na(`T value - based on AJCC 3rd (1988-2003)`) ~ "AJCC 3rd",
      !is.na(`Derived AJCC T, 6th ed (2004-2015)`) ~ "AJCC 6th",
      !is.na(`Derived SEER Combined T (2016-2017)`) ~ "SEER Combined 2016-2017",
      !is.na(`Derived EOD 2018 T Recode (2018+)`) ~ "EOD 2018",
      TRUE ~ NA_character_
    ),
    N_raw = coalesce(
      `N value - based on AJCC 3rd (1988-2003)`,
      `Derived AJCC N, 6th ed (2004-2015)`,
      `Derived SEER Combined N (2016-2017)`,
      `Derived EOD 2018 N Recode (2018+)`
    ),
    N_stage = extract_N(N_raw),
    N_major = case_when(
      str_detect(N_stage, "^N0$") ~ "N0",
      str_detect(N_stage, "^N1")  ~ "N1",
      str_detect(N_stage, "^N2")  ~ "N2",
      str_detect(N_stage, "^N3")  ~ "N3",
      str_detect(N_stage, "^NX$") ~ NA_character_,
      TRUE ~ NA_character_
    ),
    N_positive = case_when(
      N_major %in% c("N1", "N2", "N3") ~ 1,
      N_major == "N0" ~ 0,
      TRUE ~ NA_real_
    ),
    N_source = case_when(
      !is.na(`N value - based on AJCC 3rd (1988-2003)`) ~ "AJCC 3rd",
      !is.na(`Derived AJCC N, 6th ed (2004-2015)`) ~ "AJCC 6th",
      !is.na(`Derived SEER Combined N (2016-2017)`) ~ "SEER Combined 2016-2017",
      !is.na(`Derived EOD 2018 N Recode (2018+)`) ~ "EOD 2018",
      TRUE ~ NA_character_
    ),
    M_raw = coalesce(
      `M value - based on AJCC 3rd (1988-2003)`,
      `Derived AJCC M, 6th ed (2004-2015)`,
      `Derived SEER Combined M (2016-2017)`,
      `Derived EOD 2018 M Recode (2018+)`
    ),
    M_stage = extract_M(M_raw),
    M_major = case_when(
      str_detect(M_stage, "^M0$") ~ "M0",
      str_detect(M_stage, "^M1")  ~ "M1",
      str_detect(M_stage, "^MX$") ~ NA_character_,
      TRUE ~ NA_character_
    ),
    Distant_metastasis = case_when(
      M_major == "M1" ~ 1,
      M_major == "M0" ~ 0,
      TRUE ~ NA_real_
    ),
    M_source = case_when(
      !is.na(`M value - based on AJCC 3rd (1988-2003)`) ~ "AJCC 3rd",
      !is.na(`Derived AJCC M, 6th ed (2004-2015)`) ~ "AJCC 6th",
      !is.na(`Derived SEER Combined M (2016-2017)`) ~ "SEER Combined 2016-2017",
      !is.na(`Derived EOD 2018 M Recode (2018+)`) ~ "EOD 2018",
      TRUE ~ NA_character_
    ),
    LVI = recode_yes_no(`Lymph-vascular Invasion (2004+ varying by schema)`),
    PNI = recode_yes_no(`Perineural Invasion Recode (2010+)`),
    Histology_raw = `Histologic Type ICD-O-3`,
    Histology_code = str_extract(Histology_raw, "\\d{4}"),
    Histology_group = case_when(
      Histology_code %in% c("8490") |
        str_detect(str_to_lower(Histology_raw), "signet") ~ "Signet ring cell carcinoma",
      Histology_code %in% c("8480", "8481") |
        str_detect(str_to_lower(Histology_raw), "mucinous") ~ "Mucinous adenocarcinoma",
      Histology_code %in% c(
        "8140", "8141", "8143", "8144",
        "8210", "8211",
        "8220", "8221",
        "8260", "8261", "8262", "8263",
        "8570", "8574"
      ) |
        str_detect(str_to_lower(Histology_raw), "adenocarcinoma") ~ "Adenocarcinoma",
      is.na(Histology_raw) ~ NA_character_,
      TRUE ~ "Other"
    ),
    Grade_raw = coalesce(
      `Grade Recode (thru 2017)`,
      `Grade Pathological (2018+)`
    ),
    Grade_detail = case_when(
      str_detect(str_to_lower(Grade_raw), "well|grade i\\b|g1") ~ "Well differentiated",
      str_detect(str_to_lower(Grade_raw), "moderately|moderate|grade ii\\b|g2") ~ "Moderately differentiated",
      str_detect(str_to_lower(Grade_raw), "poorly|poor|grade iii\\b|g3") ~ "Poorly differentiated",
      str_detect(str_to_lower(Grade_raw), "undifferentiated|anaplastic|grade iv\\b|g4") ~ "Undifferentiated",
      str_detect(str_to_lower(Grade_raw), "low grade") ~ "Low grade",
      str_detect(str_to_lower(Grade_raw), "high grade") ~ "High grade",
      TRUE ~ NA_character_
    ),
    
    Differentiation = case_when(
      Grade_detail %in% c("Well differentiated", "Moderately differentiated", "Low grade") ~ "Well/Moderate",
      Grade_detail %in% c("Poorly differentiated", "Undifferentiated", "High grade") ~ "Poor/Undifferentiated",
      TRUE ~ NA_character_
    ),
    survival_months = suppressWarnings(as.numeric(`Survival months`)),
    
    CSS_status = case_when(
      str_detect(
        str_to_lower(`SEER cause-specific death classification`),
        "dead.*attributable|dead.*this cancer"
      ) ~ 1,
      str_detect(
        str_to_lower(`SEER cause-specific death classification`),
        "alive|other cause"
      ) ~ 0,
      TRUE ~ NA_real_
    )
  )

saveRDS(seer_clean, "data/seer_clean.rds")

seer_clean2 <- seer_clean[!is.na(seer_clean$tumor_size_cm), ]
seer_clean2 <- seer_clean2[!is.na(seer_clean2$T_major), ]
seer_clean2 <- seer_clean2[!is.na(seer_clean2$N_major), ]
seer_clean2 <- seer_clean2[!is.na(seer_clean2$M_major), ]
seer_clean2 <- seer_clean2[!is.na(seer_clean2$PNI), ]

library(dplyr)
library(stringr)
seer_clean2 <- seer_clean2 %>%
  mutate(
    Grade_raw_clean = str_squish(as.character(Grade_raw)),
    
    Differentiation_3group = case_when(
      Grade_raw_clean %in% c("1", "Well differentiated; Grade I") ~ "well",
      
      Grade_raw_clean %in% c("2", "Moderately differentiated; Grade II") ~ "moderate",
      
      Grade_raw_clean %in% c(
        "3",
        "4",
        "Poorly differentiated; Grade III",
        "Undifferentiated; anaplastic; Grade IV"
      ) ~ "poor",
      
      Grade_raw_clean %in% c("9", "Unknown", "Blank(s)") ~ NA_character_,
      
      TRUE ~ NA_character_
    ),
    
    Differentiation_3group = factor(
      Differentiation_3group,
      levels = c("well", "moderate", "poor")
    )
  )
table(seer_clean2$Differentiation_3group, useNA = "ifany")
seer_clean2 <- seer_clean2[!is.na(seer_clean2$Differentiation_3group), ]
seer_final <- seer_clean2[,c("Sex"  ,"Age recode with <1 year olds and 90+"  ,
                             "tumor_size_cm"  , "T_major"  , "N_major"  , "M_major"   ,
                             "PNI","Histology_group" ,"Differentiation_3group" , "survival_months" , "CSS_status" ,"Site rec ICD-O-3/WHO 2008 (individual sites only)")]
seer_final <- seer_final[!is.na(seer_final$CSS_status),]
seer_final <- seer_final %>%
  mutate(
    Age = as.numeric(
      str_extract(`Age recode with <1 year olds and 90+`, "\\d+")
    )
  )
table(seer_final$Age, useNA = "ifany")
saveRDS(seer_final, 'data/seer_final.rds')

library(dplyr)
library(survival)
library(survminer)

df_surv <- seer_final %>%
  filter(
    !is.na(tumor_size_cm),
    !is.na(survival_months),
    !is.na(CSS_status),
    survival_months >= 0,
    tumor_size_cm < 4 | tumor_size_cm > 7
  ) %>%
  mutate(
    SizeGroup = case_when(
      tumor_size_cm < 3 ~ "Less than 3 cm",
      tumor_size_cm > 5 ~ "Greater than 5 cm"
    ),
    SizeGroup = factor(
      SizeGroup,
      levels = c("Less than 3 cm", "Greater than 5 cm")
    )
  )

df_surv <- df_surv %>%
  filter(T_major %in% c('T4')) %>%
  filter(tumor_size_cm <14) %>%
  filter(Age < 50)


fit <- survfit(
  Surv(survival_months, CSS_status) ~ SizeGroup,
  data = df_surv
)

p <- ggsurvplot(
  fit,
  data = df_surv,
  pval = TRUE,
  conf.int = FALSE,
  risk.table = TRUE,
  xlab = "Time (months)",
  ylab = "Cancer-specific survival",
  legend.title = "Tumor size",
  legend.labs = c("Less than 3 cm", "Greater than 5 cm"),
  ggtheme = theme_classic()
)

print(p)

pdf('outs/Seer-KM.pdf',width = 6,height=8)
print(p)
dev.off()

library(writexl)

write_xlsx(seer_final, "data/seer_final.xlsx")