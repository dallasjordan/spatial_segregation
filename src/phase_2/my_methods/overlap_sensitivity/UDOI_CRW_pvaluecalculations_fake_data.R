# Calculating p-values from the "fake" simulated data julia provided
# Nov 23 2021
# Dallas Jordan

library(sjPlot)
# load in CRW UDOI values provided by Julia 
CRW_values_read <- file.choose()

CRW_values <- read.csv(CRW_values_read)

load(file.choose()) # go choose the all_test_stats_xxxOL.Rdata file

tot_UDOI95 <- ab_UDOI_95_test_stat
tot_UDOI50 <- ab_UDOI_50_test_stat

some_UDOI95 <- ab_UDOI_95_test_stat
some_UDOI50 <- ab_UDOI_50_test_stat

no_UDOI95 <- ab_UDOI_95_test_stat
no_UDOI50 <- ab_UDOI_50_test_stat


tot_UDOI95_pvalue <- (sum(CRW_values[,3]<tot_UDOI95))/100
tot_UDOI50_pvalue <- (sum(CRW_values[,6]<tot_UDOI50))/100

some_UDOI95_pvalue <- (sum(CRW_values[,3]<some_UDOI95))/100
some_UDOI50_pvalue <- (sum(CRW_values[,6]<some_UDOI50))/100

no_UDOI95_pvalue <- (sum(CRW_values[,3]<no_UDOI95))/100
no_UDOI50_pvalue <- (sum(CRW_values[,6]<no_UDOI50))/100

# Values updated as of Nov 28
tot_UDOI95_pvalue # 0.88
tot_UDOI50_pvalue # 0.34

some_UDOI95_pvalue # 0.2  # so interesting point here - the test stat is super low (UDOI = 0.007) but still 20/100 simulated were lower
some_UDOI50_pvalue # 0

no_UDOI95_pvalue # 0.45 # interesting again - UDOI value for test stat is like 0.00004484, but many randomized are lower, resulting in 0.71
                  # PHR value here is 0.0017 for the test stat, maybe this is more meaningful? 
no_UDOI50_pvalue # 0

