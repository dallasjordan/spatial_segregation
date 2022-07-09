# Calculating p-values from CRW overlap data calculated temporally and in aggregate
# CRW tracks recalculated mid May 2022
# May 1 2022
# Dallas Jordan

# Run script 3 times, once for each season
# read in CRW overlap values by category (aggregate, or by season)
# read in test stats by season
# calculate

# load in CRW overlap values calculated in CRW_pvalue_processing
CRW_values_read <- file.choose() # sim_overlap_values_SEASON.RDS in your /data/overlap_sensitivity/ folder
CRW_values <- load(CRW_values_read)
CRW_values <- bind_rows(sim_overlap_values)

# calculate test stats for real data

# load in test stats for real data
real_test_stats_read <- file.choose() # all_test_stats_SEASON.Rdata, in your "pre_defense" folder
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

ALAB95_pvalue <- (sum(CRW_values[,1]<ALAB95))/1000
ALAB50_pvalue <- (sum(CRW_values[,2]<ALAB50))/1000

MLTL95_pvalue <- (sum(CRW_values[,3]<MLTL95))/1000
MLTL50_pvalue <- (sum(CRW_values[,4]<MLTL50))/1000

MBTB95_pvalue <- (sum(CRW_values[,5]<MBTB95))/1000
MBTB50_pvalue <- (sum(CRW_values[,6]<MBTB50))/1000

MLMB95_pvalue <- (sum(CRW_values[,7]<MLMB95))/1000
MLMB50_pvalue <- (sum(CRW_values[,8]<MLMB50))/1000

TLTB95_pvalue <- (sum(CRW_values[,9]<TLTB95))/1000
TLTB50_pvalue <- (sum(CRW_values[,10]<TLTB50))/1000

# Aggregate ---------------------------------------------------------------

# all commented UDOI values here accurate May 2 2022
ALAB95_pvalue # 0
ALAB50_pvalue # 0

MLTL95_pvalue # 0.979
MLTL50_pvalue # 0.984

MBTB95_pvalue # 0.552
MBTB50_pvalue # 0.103

MLMB95_pvalue # 0.643
MLMB50_pvalue # 0.013

TLTB95_pvalue # 0
TLTB50_pvalue # 0


# June July ---------------------------------------------------------------

# all commented values here accurate May 14 2022
ALAB95_pvalue # 0
ALAB50_pvalue # 0

MLTL95_pvalue # 0.984
MLTL50_pvalue # 0

MBTB95_pvalue # 0.977
MBTB50_pvalue # 0.21

MLMB95_pvalue # 0.001
MLMB50_pvalue # 0

TLTB95_pvalue # 0
TLTB50_pvalue # 0

# August September --------------------------------------------------------

# all commented values here accurate May 14 2022
ALAB95_pvalue # 0.002
ALAB50_pvalue # 0.084

MLTL95_pvalue # 0.204
MLTL50_pvalue # 0.236

MBTB95_pvalue # 0.999
MBTB50_pvalue # 0.114

MLMB95_pvalue # 0.017
MLMB50_pvalue # 0

TLTB95_pvalue # 0.032
TLTB50_pvalue # 0

# October November --------------------------------------------------------

# all commented values here accurate May 14 2022
ALAB95_pvalue # 1
ALAB50_pvalue # 1

MLTL95_pvalue # 0.507
MLTL50_pvalue # 0.426

MBTB95_pvalue # 1
MBTB50_pvalue # 1

MLMB95_pvalue # 1
MLMB50_pvalue # 0.983

TLTB95_pvalue # 0.972
TLTB50_pvalue # 0.576


# Spring ------------------------------------------------------------------

# all commented values here accurate May 2 2022
ALAB95_pvalue # 0.001
ALAB50_pvalue # 0

MLTL95_pvalue # 0.999
MLTL50_pvalue # 0.934

MBTB95_pvalue # NA
MBTB50_pvalue # NA

MLMB95_pvalue # NA
MLMB50_pvalue # NA

TLTB95_pvalue # 0.001
TLTB50_pvalue # 0

# Summer ------------------------------------------------------------------

# all commented values here accurate May 2 2022
ALAB95_pvalue # 0
ALAB50_pvalue # 0

MLTL95_pvalue # 0.839
MLTL50_pvalue # 0.433

MBTB95_pvalue # 0.815
MBTB50_pvalue # 0.059

MLMB95_pvalue # 0
MLMB50_pvalue # 0

TLTB95_pvalue # 0
TLTB50_pvalue # 0

# Fall --------------------------------------------------------------------

# all commented values here accurate May 2 2022
ALAB95_pvalue # 1
ALAB50_pvalue # 0.885

MLTL95_pvalue # 0.585
MLTL50_pvalue # 0.424

MBTB95_pvalue # 1
MBTB50_pvalue # 1

MLMB95_pvalue # 1
MLMB50_pvalue # 0.149

TLTB95_pvalue # 1
TLTB50_pvalue # 0.954

