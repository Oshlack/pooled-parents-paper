library('ggplot2')
#devtools::install_github("lazappi/mcriPalettes")
library("mcriPalettes")
my.colours = c(1,4,5)

# Either source to re-generate stats, or just load previous data
source('R/simulation_stats.R') # Uses about 80GB memory!
#load('data/simulation_summary_data.Robj')

plot_prefix = 'plots/compare_callers.'

### Plots functions with subsets

# Plot recall as a percentage
plot_recall = function(variantcallers, sim_types, label){
  data_subset = subset(poolstats, variantcaller %in% variantcallers & sim_type %in% sim_types)
  ggplot(data=data_subset, aes(x=pool, y=percent_recovered, colour = variantcaller, shape = sim_type)) + 
    stat_summary(fun.y=mean,geom="line",aes(group=interaction(sim_type,variantcaller), 
                                            colour=variantcaller,
                                            linetype = sim_type), size=1.5) +
    geom_point(size=2.5) +
    scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 20)) +
    labs(x="Number of samples in pool", y="Recall %") +
    scale_color_manual(values = mcriPalette("symbol")[my.colours]) +
    theme_classic(base_size = 18) +
    theme(legend.position=c(0.2, 0.3), legend.background=element_blank())
  ggsave(paste0(plot_prefix, label, 'recall_percent.jpg'), height=8, width=8)
}

plot_recall_single = function(variantcallers, sim_types, label){
  data_subset = subset(poolstats_alleles, variantcaller %in% variantcallers & sim_type %in% sim_types 
                       & nonref_allele_count_truth == 1)
  ggplot(data=data_subset, 
         aes(x=pool, y=percent_recovered, colour = variantcaller, shape = sim_type)) + 
    geom_point(size=2.5) + 
    stat_summary(fun.y=mean,geom="line",aes(group=interaction(sim_type,variantcaller), 
                                            colour=variantcaller,
                                            linetype = sim_type), size=1.5) + 
    scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 20)) +
    labs(x="Number of samples in pool", y="Recall %") +
    scale_color_manual(values = mcriPalette("symbol")[my.colours]) + 
    theme_classic(base_size = 18) +
    theme(legend.position='none')
  ggsave(paste0(plot_prefix, label, 'recall_percent_singles.jpg'), height=8, width=8)
}

plot_fp = function(variantcallers, sim_types, label){
  data_subset = subset(poolstats_fp, variantcaller %in% variantcallers & sim_type %in% sim_types)
  ggplot(data=data_subset, aes(x=pool, y=percent_fp, colour = variantcaller, shape = sim_type)) +
    stat_summary(fun.y=mean,geom="line",aes(group=interaction(sim_type,variantcaller), 
                                            colour=variantcaller,
                                            linetype = sim_type), size=1.5) +
    scale_y_continuous(limits = c(0, 39), breaks = seq(0, 100, by = 10)) +
    labs(x="Number of samples in pool", y="False positives %") +
    scale_color_manual(values = mcriPalette("symbol")[my.colours]) +
    theme_classic(base_size = 18) +
    theme(legend.position='none')
  ggsave(paste0(plot_prefix, label, 'false_pos_percent.jpg'), height=8, width=8)
}


### Plots for paper:
plot_recall(c("gatk individual", "gatk joint", "freebayes highqual"),
                   c("constant depth", "additive depth"),
                   label = 'all.')
plot_recall_single(c("gatk individual", "gatk joint", "freebayes highqual"),
                   c("constant depth", "additive depth"),
                   label = 'all.')

plot_fp(c("gatk individual", "gatk joint", "freebayes highqual"),
        c("constant depth", "additive depth"),
        label = 'all.')


ggplot(data=subset(poolstats, sim_type %in% 'constant depth'),
       aes(y=total, x=pool, colour=variantcaller)) +
  stat_summary(fun.y=mean,geom="line",aes(group=variantcaller), size=1.5) +
  geom_point(position='dodge', size=3, aes(shape=variantcaller)) +
  labs(x="Number of samples in pool", y="Mean variants called per pool") +
  scale_color_manual(values = mcriPalette("symbol")[my.colours]) +
  theme_classic(base_size = 20) + theme(legend.position='none') +
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE))
ggsave(paste0(plot_prefix, 'total_variants.jpg'), height=8, width=8)


# Summary stats for paper
library('dplyr')
library(tidyr)

poolstats %>% group_by(pool, sim_type, variantcaller) %>% 
  summarise(mean_recall = mean(percent_recovered)) -> poolstats_mean
poolstats_mean_wide = spread(poolstats_mean, key = pool, value = mean_recall)
poolstats_mean_wide[,c('2','4','6','8','10')] = round(poolstats_mean_wide[,c('2','4','6','8','10')],1)
write.csv(poolstats_mean_wide, 'simulation_mean_recall_all.csv')

subset(poolstats_alleles, nonref_allele_count_truth == 1) %>% 
  group_by(pool, sim_type, variantcaller) %>% 
  summarise(mean_recall = mean(percent_recovered)) -> poolstats_alleles_mean
poolstats_alleles_mean_wide = spread(poolstats_alleles_mean, key = pool, value = mean_recall)
poolstats_alleles_mean_wide[,c('2','4','6','8','10')] = round(poolstats_alleles_mean_wide[,c('2','4','6','8','10')],1)
write.csv(poolstats_alleles_mean_wide, 'simulation_mean_recall_single.csv')

poolstats_fp %>% group_by(pool, sim_type, variantcaller) %>% 
  summarise(mean_fp = mean(percent_fp)) -> poolstats_fp_mean
poolstats_fp_mean_wide = spread(poolstats_fp_mean, key = pool, value = mean_fp)
poolstats_fp_mean_wide[,c('2','4','6','8','10')] = signif(poolstats_fp_mean_wide[,c('2','4','6','8','10')],3)
write.csv(poolstats_fp_mean_wide, 'simulation_mean_fp_all.csv')

poolstats_fp %>% group_by(pool, sim_type, variantcaller) %>% 
  summarise(mean_total = mean(total)) -> pool_totals_mean
pool_totals_mean_wide = spread(pool_totals_mean, key = pool, value = mean_total)
pool_totals_mean_wide[,c('2','4','6','8','10')] = round(pool_totals_mean_wide[,c('2','4','6','8','10')],0)
write.csv(pool_totals_mean_wide, 'simulation_total_per_pool.csv')

# Use the pools of 10 as the largest, so capturing most individuals
poolstats$total_per_ind = poolstats$total/poolstats$pool
subset(poolstats, pool == '10') %>% group_by(sim_type, variantcaller) %>% 
  summarise(mean_total_per_ind = mean(total_per_ind)) -> poolstats_total_per_ind_mean
poolstats_total_per_ind_mean_wide = spread(poolstats_total_per_ind_mean, key = variantcaller, value = mean_total_per_ind)
poolstats_total_per_ind_mean_wide[,c('gatk individual','gatk joint','freebayes highqual')] = round(poolstats_total_per_ind_mean_wide[,c('gatk individual','gatk joint','freebayes highqual')],0)
write.csv(poolstats_total_per_ind_mean_wide, 'simulation_total_per_ind_all.csv')

ggplot(data= poolstats, aes(x = pool, y = total_per_ind,
    colour = variantcaller, shape = sim_type)) +
  geom_point()



