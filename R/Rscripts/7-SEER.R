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




