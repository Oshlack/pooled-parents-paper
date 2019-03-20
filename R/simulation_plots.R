start.time <- Sys.time()

library('ggplot2')
library('dplyr')
library("data.table")
#devtools::install_github("lazappi/mcriPalettes")
library("mcriPalettes")
my.colours = c(1,4,5)

source('R/pooled_paper_functions.R')


plot_prefix = 'plots/compare_callers.'
variantcaller_order = c("gatk individual", "gatk joint", "freebayes highqual")
sim_type_order = c("constant depth", "additive depth")

## Parse variant comparisons from individually called data

### GATK individual calling (unfiltered)
data_dir = 'data/gatk_individual/'

## Constant Depth
replicates = c('random_pools', 'random_pools2', 'random_pools3')
all_data_ind_const = extract_ind_csvs(data_dir, replicates, 
                                      sim_type = 'constant depth')
#all_data_ind_const = all_data_ind_const[all_data_ind_const$pool < 20,] # remove large pools
all_data_ind_const$variantcaller = "gatk individual"

# Separate into variants called in probands calculate false positives, then delete original
all_data_ind_const_probands = all_data_ind_const[!is.na(all_data_ind_const$proband),]
poolstats_ind_const_probands_fp = falsepos_rate(all_data_ind_const)
poolstats_ind_const_probands_fp$variantcaller = "gatk individual"
poolstats_ind_const_probands_fp$sim_type = 'constant depth'
rm(all_data_ind_const)

## additive Depth
replicates = c('depth_pools', 'depth_pools2', 'depth_pools3')
all_data_ind_prop = extract_ind_csvs(data_dir, replicates, 
                                     sim_type = 'additive depth')
all_data_ind_prop$variantcaller = "gatk individual"
# Separate into variants called in probands calculate false positives, then delete original
all_data_ind_prop_probands = all_data_ind_prop[!is.na(all_data_ind_prop$proband),]
poolstats_ind_prop_probands_fp = falsepos_rate(all_data_ind_prop)
poolstats_ind_prop_probands_fp$variantcaller = "gatk individual"
poolstats_ind_prop_probands_fp$sim_type = 'additive depth'
rm(all_data_ind_prop)

### GATK joint calling (unfiltered)
data_dir = 'data/gatk_joint/'

# Constant Depth
replicates = c('random_pools_joint', 'random_pools2_joint', 'random_pools3_joint')#, 
#'random_pools4_joint', 'random_pools5_joint', 'random_pools6_joint')
all_data_joint_const = extract_multi_csvs(data_dir, replicates, 
                                          sim_type = 'constant depth')
all_data_joint_const$variantcaller = "gatk joint"
#all_data1 = all_data1 [all_data1$total_alleles_probands == all_data1$pool*2,] # Only keep loci called in all probands
all_data_joint_const_probands = all_data_joint_const[!is.na(all_data_joint_const$proband),]
poolstats_joint_const_probands_fp = falsepos_rate(all_data_joint_const)
poolstats_joint_const_probands_fp$variantcaller = "gatk joint"
poolstats_joint_const_probands_fp$sim_type = 'constant depth'
rm(all_data_joint_const)

# Proportinal Depth
replicates = c('depth_pools_joint', 'depth_pools2_joint', 'depth_pools3_joint')#,
#'depth_pools4_joint', 'depth_pools5_joint', 'depth_pools6_joint')
all_data_joint_prop = extract_multi_csvs(data_dir, replicates,
                                         sim_type = 'additive depth')
all_data_joint_prop$variantcaller = "gatk joint"
#all_data1 = all_data1 [all_data1$total_alleles_probands == all_data1$pool*2,] # Only keep loci called in all probands
all_data_joint_prop_probands = all_data_joint_prop[!is.na(all_data_joint_prop$proband),]
poolstats_joint_prop_probands_fp = falsepos_rate(all_data_joint_prop)
poolstats_joint_prop_probands_fp$variantcaller = "gatk joint"
poolstats_joint_prop_probands_fp$sim_type = 'additive depth'
rm(all_data_joint_prop)

### freebayes individual calling (QUAL > 20 filtered)
data_dir = 'data/freebayes/'

## Constant Depth freebayes
replicates = c('random_pools_joint', 'random_pools2_joint', 'random_pools3_joint')
# replicates = c('random_pools_joint', 'random_pools2_joint', 'random_pools3_joint',
#                'random_pools4_joint', 'random_pools5_joint', 'random_pools6_joint')
all_data_fb_const = extract_ind_csvs(data_dir, replicates,
                               sim_type = 'constant depth')
all_data_fb_const$variantcaller = "freebayes highqual"

# Separate into variants called in probands calculate false positives, then delete original
all_data_fb_const_probands = all_data_fb_const[!is.na(all_data_fb_const$proband),]
poolstats_fb_const_probands_fp = falsepos_rate(all_data_fb_const)
poolstats_fb_const_probands_fp$variantcaller = "freebayes highqual"
poolstats_fb_const_probands_fp$sim_type = 'constant depth'
#rm(all_data_fb_const)

### Calculate stats

# Calculate recall over all variants in all probands 
# from a file containing many different simulation types and variant callers
recall_all_probands_multi = function(variants) {
  variants %>% group_by(pool, simulation, sim_type, variantcaller) %>% 
    summarise(n_recovered = sum(recovered_proband, na.rm = T), total = length(recovered_proband)) -> vars_alleles
  vars_alleles$percent_recovered = vars_alleles$n_recovered/vars_alleles$total*100
  return(vars_alleles)
}

# Group by number of non-reference alleles at locus for recall 
# from a file containing many different simulation types and variant callers
group_nonref_alleles_multi = function(variants) {
  variants %>% group_by(pool, simulation, sim_type, variantcaller, nonref_allele_count_truth) %>% 
    summarise(n_recovered = sum(recovered_proband), total = length(recovered_proband)) -> vars_alleles
  vars_alleles$percent_recovered = vars_alleles$n_recovered/vars_alleles$total*100
  return(vars_alleles)
}

# Make row names consistent then commbine data
setnames(all_data_ind_const_probands, old = c('recovered'), new = c('recovered_all'))
setnames(all_data_ind_prop_probands, old = c('recovered'), new = c('recovered_all'))
setnames(all_data_joint_const_probands, old = c('recovered'), new = c('recovered_proband'))
setnames(all_data_joint_prop_probands, old = c('recovered'), new = c('recovered_proband'))
setnames(all_data_fb_const_probands, old = c('recovered'), new = c('recovered_all'))

shared_names = intersect(names(all_data_ind_const_probands), names(all_data_joint_const_probands))

# all_data_probands = rbind(all_data_ind_const_probands[,shared_names], all_data_ind_prop_probands[,shared_names], 
#                           all_data_joint_const_probands[,shared_names], all_data_joint_prop_probands[,shared_names])
all_data_probands = rbind(all_data_ind_const_probands[,shared_names], all_data_ind_prop_probands[,shared_names],
                          all_data_joint_const_probands[,shared_names], all_data_joint_prop_probands[,shared_names],
                          all_data_fb_const_probands[,shared_names])
# Set factors to help with plotting XXX doesn't work when done here?
all_data_probands$variantcaller = factor(as.character(all_data_probands$variantcaller), levels = variantcaller_order)
poolstats = recall_all_probands_multi(all_data_probands)
poolstats_alleles = group_nonref_alleles_multi(all_data_probands)

poolstats_fp = rbind(poolstats_ind_const_probands_fp, poolstats_ind_prop_probands_fp,
                     poolstats_joint_const_probands_fp, poolstats_joint_prop_probands_fp,
                     poolstats_fb_const_probands_fp)

# Set plotting order
poolstats$variantcaller = factor(as.character(poolstats$variantcaller), levels = variantcaller_order)
poolstats_alleles$variantcaller = factor(as.character(poolstats_alleles$variantcaller), levels = variantcaller_order)
poolstats_fp$variantcaller = factor(as.character(poolstats_fp$variantcaller), levels = variantcaller_order)
poolstats$sim_type = factor(as.character(poolstats$sim_type), levels = sim_type_order)
poolstats_alleles$sim_type = factor(as.character(poolstats_alleles$sim_type), levels = sim_type_order)
poolstats_fp$sim_type = factor(as.character(poolstats_fp$sim_type), levels = sim_type_order)

#Save objects
# save(file = 'pool_sim_data.Robj',
#      all_data_ind_const_probands, all_data_ind_prop_probands,
#      all_data_joint_const_probands, all_data_joint_prop_probands,
#      poolstats_ind_const_probands_fp, poolstats_ind_prop_probands_fp,
#      poolstats_joint_const_probands_fp, poolstats_joint_prop_probands_fp,
#      all_data_fb_const_probands,poolstats_fb_const_probands_fp,
#      all_data_probands
#      )
# save(file = 'pool_sim_poolstats_data.Robj', poolstats, poolstats_alleles, all_data_probands)
# load('pool_sim_data.Robj')
# load('pool_sim_poolstats_data.Robj')
#for (thing in ls()) { message(thing); print(object.size(get(thing)), units='auto') }


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

end.time <- Sys.time()
time.taken <- end.time - start.time
print(time.taken)
