library('ggplot2')
library('dplyr')
# library("data.table")
#devtools::install_github("lazappi/mcriPalettes")
library("mcriPalettes")
my.colours = c(1,4,5)

source('pooled_paper_functions.R')

### Functions

### GATK individual calling
data_dir = 'data/'
f = 'Parental-pool.compare.csv'

## Constant Depth
all_data_real = read.csv(paste0(data_dir, f), stringsAsFactors = F)
all_data_real$variantcaller = "gatk individual"
# Note, no false positives in this data

# Calculate recall
all_data_real %>% group_by(pool) %>% 
  summarise(n_recovered = sum(recovered_proband, na.rm = T), total = length(recovered_proband)) -> poolstats_real
poolstats_real$percent_recovered = poolstats_real$n_recovered/poolstats_real$total*100

# Cacluate recall by allele count
all_data_real %>% group_by(pool, nonref_alleles_probands) %>% 
  summarise(n_recovered = sum(recovered_proband, na.rm = T), total = length(recovered_proband)) -> poolstats_alleles_real
poolstats_alleles_real$percent_recovered = poolstats_alleles_real$n_recovered/poolstats_alleles_real$total*100

poolstats_real
poolstats_alleles_real

###### Calculate some numbers after filtering
table(all_data_real$proband)

# Starting number of variants that are in any of the probands (not just the pool) and have QD>=2
proband_vars = all_data_real[all_data_real$QD_proband >= 2,]
table(proband_vars[,"proband"])
# Also remove variants where more than 2 nonref alleles found in the probands
proband_vars = all_data_real[all_data_real$nonref_alleles_probands <=2,]
table(proband_vars[,"proband"])
# Limit to gnomad very rare: 0.0005 or NA
filter = proband_vars$AF_EXOMESgnomad < 0.0005 | is.na(proband_vars$AF_EXOMESgnomad)
table(proband_vars[filter,"proband"])
# Limit to not recovered in pool
filter = !proband_vars$recovered_proband
table(proband_vars[filter,"proband"])
# Limit to not recovered in pool and gnomad AF_EXOMESgnomad very rare: 0.0005 or NA
pool_gnomad_filter = !proband_vars$recovered_proband & (proband_vars$AF_EXOMESgnomad < 0.0005 | is.na(proband_vars$AF_EXOMESgnomad))
table(proband_vars[pool_gnomad_filter,"proband"])

