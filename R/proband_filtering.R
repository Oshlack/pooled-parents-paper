library('ggplot2')
library('dplyr')
# library("data.table")
#devtools::install_github("lazappi/mcriPalettes")
library("mcriPalettes")
my.colours = c(1,4,5)

source('R/simulation_plots.R') # Potentially takes ~1 hour if doing simulation stats from scratch!
source('R/pooled_paper_functions.R')

### Functions

### GATK individual calling
data_dir = 'data/real_analysis/'
f = 'Parental-pool.compare.csv'
#f = 'Parental-pool.compare_notsplit.csv'

## Parse data
all_data_real = read.csv(paste0(data_dir, f), stringsAsFactors = F)
all_data_real$variantcaller = "gatk individual"
all_data_real$pool = 4
# Calculate nonref_allele_count_truth over all probands
group_by(all_data_real, position) %>% 
  summarise(nonref_alleles_probands = sum(nonref_alleles_proband, na.rm = T),
            total_alleles_probands = sum(total_alleles_proband, na.rm = T),
            nonref_reads_probands = sum(nonref_reads_proband, na.rm = T)
  ) -> summed_alleles
all_data_real = merge(all_data_real, summed_alleles, all.x = T)
all_data_real$recovered = all_data_real$recovered_proband
# Notes:
# No false positives in this data 
# This is unfiltered vcfs, so doesn't take into account different capture regions

# Calculate recall
all_data_real %>% group_by(pool) %>% 
  summarise(n_recovered = sum(recovered_proband, na.rm = T), total = length(recovered_proband)) -> poolstats_real
poolstats_real$percent_recovered = poolstats_real$n_recovered/poolstats_real$total*100

# Cacluate recall by allele count
all_data_real %>% group_by(pool, nonref_alleles_probands) %>% 
  summarise(n_recovered = sum(recovered_proband, na.rm = T), total = length(recovered_proband)) -> poolstats_alleles_real
poolstats_alleles_real$percent_recovered = poolstats_alleles_real$n_recovered/poolstats_alleles_real$total*100

print(poolstats_real)
print(poolstats_alleles_real[1,])

real_recall = round(poolstats_real[[1,'percent_recovered']], 1)
real_recall_singleton = round(poolstats_alleles_real[[1,'percent_recovered']], 1)

# Plot real data recall over simulations
data_subset = subset(poolstats, sim_type == 'constant depth')
ggplot(data=data_subset, aes(x=pool, y=percent_recovered, colour = variantcaller)) + 
  stat_summary(fun.y=mean,geom="line",aes(group=interaction(sim_type,variantcaller), 
                                          colour=variantcaller), size=1.5) +
  geom_point(size=2.5) +
  geom_point(x=8, y=real_recall, colour='black', size=10, shape='*') +
  annotate("text", x=8.5, y=real_recall, label=paste0(real_recall,"%"), size=5) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 20)) +
  labs(x="Number of samples in pool", y="Recall %") +
  scale_color_manual(values = mcriPalette("symbol")[my.colours]) +
  theme_classic(base_size = 18) +
  theme(legend.position=c(0.2, 0.15), legend.background=element_blank())
ggsave('plots/recall_constant_plus_real.jpg', height=8, width=8)

data_subset = subset(poolstats_alleles, sim_type == 'constant depth' 
                     & nonref_allele_count_truth == 1)
ggplot(data=data_subset, 
       aes(x=pool, y=percent_recovered, colour = variantcaller)) + 
  geom_point(size=2.5) + 
  geom_point(x=8, y=real_recall_singleton, colour='black', size=10, shape='*') +
  annotate("text", x=8.5, y=real_recall_singleton, label=paste0(real_recall_singleton,"%"), size=5) +
  stat_summary(fun.y=mean,geom="line",aes(group=interaction(sim_type,variantcaller), 
                                          colour=variantcaller), size=1.5) + 
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 20)) +
  labs(x="Number of samples in pool", y="Recall %") +
  scale_color_manual(values = mcriPalette("symbol")[my.colours]) +
  theme_classic(base_size = 18) +
  theme(legend.position='none')
ggsave('plots/recall_singles_constant_plus_real.jpg', height=8, width=8)




# Try to match up variants between the pool recall data and the vep data
all_data_real$location_str = sapply(all_data_real$variant, variant_to_location)
all_data_real$proband = formatC(all_data_real$proband, width = 9, format = "d", flag = "0")
all_data_real$in_pool = all_data_real$recovered_proband

# get variant annotations
all_vep = extract_vep_tsvs(data_dir)
# Reduce the number of consequences by categorising them
all_vep$Consequence = sapply(all_vep$Consequence, simplify_consequence)
all_vep$Consequence[all_vep$Consequence %in% c('5_prime_UTR', '3_prime_UTR')] = "UTR"
all_vep$Consequence[all_vep$Consequence %in% c('mature_miRNA', 'non_coding_transcript_exon')] = "non_coding_transcript"
all_vep$Consequence[all_vep$Consequence %in% c('downstream_gene', 'upstream_gene', 'intergenic')] = "intergenic"

all_vep$Consequence = factor_by_freq(all_vep$Consequence)
all_vep$IMPACT = factor(all_vep$IMPACT, levels = c('HIGH', 'MODERATE', 'LOW', 'MODIFIER'))

all_merged = merge(all_data_real, all_vep, all.x=T) # Should I do all=T or all.x=T instead?
all_merged$gnomAD_AF_joint = pmax(all_merged$gnomAD_AF, all_merged$AF_EXOMESgnomad, na.rm = T)
#all_merged$gnomAD_AF_joint = all_merged$gnomAD_AF
#all_merged$gnomAD_AF_joint = all_merged$AF_EXOMESgnomad

all_merged_intersect = subset(all_merged, !is.na(Consequence))

gnomad_filtered = subset(all_merged_intersect, gnomAD_AF_joint < 0.0005 | is.na(gnomAD_AF_joint))
#gnomad_filtered = subset(all_merged_intersect, gnomAD_AF < 0.0005 | is.na(gnomAD_AF))
gnomad_filtered$Consequence = factor_by_freq(gnomad_filtered$Consequence)

# Calculate recall
gnomad_filtered %>% group_by(pool) %>% 
  summarise(n_recovered = sum(recovered_proband, na.rm = T), total = length(recovered_proband)) -> gnomad_poolstats
gnomad_poolstats$percent_recovered = gnomad_poolstats$n_recovered/gnomad_poolstats$total*100

# Cacluate recall by allele count
gnomad_filtered %>% group_by(pool, nonref_alleles_probands) %>% 
  summarise(n_recovered = sum(recovered_proband, na.rm = T), total = length(recovered_proband)) -> gnomad_poolstats_alleles
gnomad_poolstats_alleles$percent_recovered = gnomad_poolstats_alleles$n_recovered/gnomad_poolstats_alleles$total*100

# Recall results for genomAD filtered variants
gnomad_poolstats
gnomad_poolstats_alleles[1,]

ggplot(data=all_merged_intersect,
       aes(x=Consequence, fill = in_pool)) +
  geom_bar() + coord_flip() +
  facet_wrap(~proband) +
  theme_classic(base_size = 16)
#ggsave('plots/consequence_barchart.jpg', height=8, width=10)









###### Calculate some numbers after filtering
table(all_merged$proband)

# Starting number of variants that are in any of the probands (not just the pool) and have QD>=2
proband_vars = all_merged[all_merged$QD >= 2,]
table(proband_vars[,"proband"])
# Also remove variants where more than 2 nonref alleles found in the probands
proband_vars = proband_vars[proband_vars$nonref_alleles_probands <=2,]
table(proband_vars[,"proband"])
# Limit to gnomad very rare: 0.0005 or NA
gnomad_filter = (proband_vars$gnomAD_AF_joint < 0.0005 | is.na(proband_vars$gnomAD_AF_joint))
table(proband_vars[gnomad_filter,"proband"])
# Limit to not recovered in pool
pool_filter = !proband_vars$recovered_proband
table(proband_vars[pool_filter,"proband"])
# Limit to not recovered in pool and gnomad AF_EXOMESgnomad very rare: 0.0005 or NA
pool_gnomad_filter = gnomad_filter & pool_filter
table(proband_vars[pool_gnomad_filter,"proband"])

# pooled parent + gnomAD filters / gnomAD only filter (as percentage)
table(proband_vars[pool_gnomad_filter,"proband"])/table(proband_vars[gnomad_filter,"proband"])*100

# What percentage of all proband variants are filtered out by gnomAD filter only
gnomad_filter_all = (1-table(all_merged[(proband_vars$gnomAD_AF_joint < 0.0005 | is.na(proband_vars$gnomAD_AF_joint)),
                 "proband"])/table(all_merged[,"proband"]))*100
round(gnomad_filter_all, 1)
# Europeans
round(mean(gnomad_filter_all[c(1,2,4)]), 1)
# Pacific Islander
round(gnomad_filter_all[c(3)], 1)

# Plot consequences
filtered_qual_multi_gnomad = proband_vars[gnomad_filter & !is.na(proband_vars$Consequence),]
filtered_qual_multi_gnomad$Consequence = factor_by_freq(filtered_qual_multi_gnomad$Consequence)

# Rename samples to anonymous IDs
old_names = sort(unique(filtered_qual_multi_gnomad$proband))
new_names = c('Proband 1', 'Proband 2', 'Proband 3', 'Proband 4')
filtered_qual_multi_gnomad$proband_anon = c(new_names, filtered_qual_multi_gnomad$proband)[match(filtered_qual_multi_gnomad$proband, c(old_names, filtered_qual_multi_gnomad$proband))]

ggplot(data=filtered_qual_multi_gnomad,
       aes(x=Consequence, fill = in_pool)) +
  geom_bar() + coord_flip() +
  facet_wrap(~proband_anon) +
  theme_classic(base_size = 16)
ggsave('plots/consequence_barchart_gnomad_filtered.jpg', height=8, width=10)

ggplot(data=subset(filtered_qual_multi_gnomad, 
                   !Consequence %in% c('synonymous', 'missense', 'non_coding_transcript')),
       aes(x=Consequence, fill = in_pool)) +
  geom_bar() + coord_flip() +
  facet_wrap(~proband_anon) +
  theme_classic(base_size = 16)
ggsave('plots/consequence_barchart_gnomad_filtered_zoom.jpg', height=8, width=10)

ggplot(data=gnomad_filtered,
       aes(x=IMPACT, fill = in_pool)) +
  geom_bar() +
  facet_wrap(~proband) +
  theme_classic(base_size = 16)
#ggsave('plots/impact_barchart_gnomad_filtered.jpg', height=8, width=10)



