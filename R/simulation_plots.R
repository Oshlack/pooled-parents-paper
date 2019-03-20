library('ggplot2')
#devtools::install_github("lazappi/mcriPalettes")
library("mcriPalettes")
my.colours = c(1,4,5)

# Either source to re-generate stats, or just load previous data
source('R/simulation_stats.R') # Uses about 80GB memory!
load('data/simulation_summary_data.Robj')

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
  labs(x="Number of samples in pool", y="Mean variants called per pool") +
  scale_color_manual(values = mcriPalette("symbol")[my.colours]) +
  theme_classic(base_size = 20) + theme(legend.position='none')
ggsave(paste0(plot_prefix, 'total_variants.jpg'), height=8, width=8)

# Plot real data recall over simulations
data_subset = subset(poolstats, sim_type == 'constant depth')
ggplot(data=data_subset, aes(x=pool, y=percent_recovered, colour = variantcaller)) + 
  stat_summary(fun.y=mean,geom="line",aes(group=interaction(sim_type,variantcaller), 
                                          colour=variantcaller), size=1.5) +
  geom_point(size=2.5) +
  geom_point(x=8, y=78, colour='black', size=10, shape='*') +
  annotate("text", x=8.5, y=78, label= "78%", size=5) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 20)) +
  labs(x="Number of samples in pool", y="Recall %") +
  scale_color_manual(values = mcriPalette("symbol")[my.colours]) +
  theme_classic(base_size = 18) +
  theme(legend.position='none')
ggsave('plots/recall_constant_plus_real.jpg', height=8, width=8)

data_subset = subset(poolstats_alleles, sim_type == 'constant depth' 
                     & nonref_allele_count_truth == 1)
ggplot(data=data_subset, 
       aes(x=pool, y=percent_recovered, colour = variantcaller)) + 
  geom_point(size=2.5) + 
  geom_point(x=8, y=72, colour='black', size=10, shape='*') +
  annotate("text", x=8.5, y=72, label= "72%", size=5) +
  stat_summary(fun.y=mean,geom="line",aes(group=interaction(sim_type,variantcaller), 
                                          colour=variantcaller), size=1.5) + 
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 20)) +
  labs(x="Number of samples in pool", y="Recall %") +
  scale_color_manual(values = mcriPalette("symbol")[my.colours]) +
  theme_classic(base_size = 18) +
  theme(legend.position='none')
ggsave('plots/recall_singles_constant_plus_real.jpg', height=8, width=8)