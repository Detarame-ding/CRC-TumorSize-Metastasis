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