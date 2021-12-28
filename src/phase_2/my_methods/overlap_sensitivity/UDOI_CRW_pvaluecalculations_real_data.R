# Calculating p-values from data julia provided
# CRW tracks from the REAL data
# Nov 23 2021
# Dallas Jordan

library(sjPlot)
# load in CRW UDOI values provided by Julia 
CRW_values_read <- file.choose()
CRW_values <- read.csv(CRW_values_read)
CRW_values <- CRW_values[,3:12]

# load in test stats for real data
real_test_stats_read <- file.choose()
test_stats <- load(real_test_stats_read)

ALAB95 <- alab_UDOI_95_test_stat
ALAB50 <- alab_UDOI_50_test_stat

MLTL95 <- mltl_UDOI_95_test_stat
MLTL50 <- mltl_UDOI_50_test_stat

MBTB95 <- mbtb_UDOI_95_test_stat
MBTB50 <- mbtb_UDOI_50_test_stat

MLMB95 <- mlmb_UDOI_95_test_stat
MLMB50 <- mlmb_UDOI_50_test_stat

TLTB95 <- tltb_UDOI_95_test_stat
TLTB50 <- tltb_UDOI_50_test_stat

ALAB95_pvalue <- (sum(CRW_values[,1]<ALAB95))/100
ALAB50_pvalue <- (sum(CRW_values[,2]<ALAB50))/100

MLTL95_pvalue <- (sum(CRW_values[,3]<MLTL95))/100
MLTL50_pvalue <- (sum(CRW_values[,4]<MLTL50))/100

MBTB95_pvalue <- (sum(CRW_values[,5]<MBTB95))/100
MBTB50_pvalue <- (sum(CRW_values[,6]<MBTB50))/100

MLMB95_pvalue <- (sum(CRW_values[,7]<MLMB95))/100
MLMB50_pvalue <- (sum(CRW_values[,8]<MLMB50))/100

TLTB95_pvalue <- (sum(CRW_values[,9]<TLTB95))/100
TLTB50_pvalue <- (sum(CRW_values[,10]<TLTB50))/100

# all commented values here accurate Nov 30 2021
ALAB95_pvalue # 0
ALAB50_pvalue # 0

MLTL95_pvalue # 0.83
MLTL50_pvalue # 0.83

MBTB95_pvalue # 0.5
MBTB50_pvalue # 0.02

MLMB95_pvalue # 0.91
MLMB50_pvalue # 0 

TLTB95_pvalue # 0 
TLTB50_pvalue # 0 
