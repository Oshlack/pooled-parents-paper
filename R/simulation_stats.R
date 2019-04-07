start.time <- Sys.time()


library('dplyr')
library("data.table")

source('R/pooled_paper_functions.R')

variantcaller_order = c("gatk individual", "gatk joint", "freebayes highqual")
sim_type_order = c("constant depth", "additive depth")

## Parse variant comparisons from individually called data

### GATK individual calling (unfiltered)
data_dir = 'data/gatk_individual/'

## Constant Depth
replicates = c('const_rep1', 'const_rep2', 'const_rep3')
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
replicates = c('depth_rep1', 'depth_rep2', 'depth_rep3')
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
replicates = c('random_pools_joint', 'random_pools2_joint', 'random_pools3_joint')
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
replicates = c('depth_pools_joint', 'depth_pools2_joint', 'depth_pools3_joint')
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
replicates = c('const_rep1', 'const_rep2', 'const_rep3')
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

# Save R objects
save(file = 'data/simulation_summary_data.Robj', poolstats, poolstats_alleles, poolstats_fp)

end.time <- Sys.time()
time.taken <- end.time - start.time
print('Total time for simulations_stats.R')
print(time.taken)