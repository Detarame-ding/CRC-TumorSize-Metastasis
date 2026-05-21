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

