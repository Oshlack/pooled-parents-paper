library('dplyr')

# Annotate with multi-alleleic and indel
check.allele = function(variant.id) {
  variant.id = as.character(variant.id)
  split.alleles = unlist(strsplit(variant.id, '_'))
  ref = split.alleles[length(split.alleles ) - 1]
  alt = split.alleles[length(split.alleles )]
  split.alt = unlist(strsplit(alt, "/"))
  # Check if multialleleic
  multi = length(split.alt) > 1
  # Check if indel by checking if any of the alleles are different sizes
  indel = diff(range(nchar(c(ref,split.alt)))) > 0
  
  result = c(multi, indel)
  names(result) = c('multi', 'indel')
  return(result)
}

# Group by number of non-reference alleles at locus for recall
group_nonref_alleles = function(variants, variantcaller) {
  variants %>% group_by(pool, nonref_allele_count_truth) %>% 
    summarise(n_recovered = sum(recovered), total = length(recovered)) -> vars_alleles
  vars_alleles$percent_recovered = vars_alleles$n_recovered/vars_alleles$total*100
  vars_alleles$pool_id = paste0('Pool_', vars_alleles$pool)
  vars_alleles$variantcaller = variantcaller
  return(vars_alleles)
  
  # Group by number of non-reference alleles at locus for recall
  group_nonref_alleles_named = function(variants, variantcaller) {
    variants %>% group_by(pool, simulation, nonref_allele_count_truth) %>% 
      summarise(n_recovered = sum(recovered), total = length(recovered)) -> vars_alleles
    vars_alleles$percent_recovered = vars_alleles$n_recovered/vars_alleles$total*100
    vars_alleles$pool_id = paste0('Pool_', vars_alleles$pool)
    vars_alleles$variantcaller = variantcaller
    return(vars_alleles)
  }
}

# Group by number of non-reference alleles at locus for recall
group_proband_nonref_alleles = function(variants, variantcaller) {
  variants %>% group_by(pool, nonref_alleles_probands) %>% 
    summarise(n_recovered = sum(recovered_all), total = length(recovered_all)) -> vars_alleles
  vars_alleles$percent_recovered = vars_alleles$n_recovered/vars_alleles$total*100
  vars_alleles$pool_id = paste0('Pool_', vars_alleles$pool)
  vars_alleles$variantcaller = variantcaller
  return(vars_alleles)
}

# Group by number of non-reference alleles at locus for recall for read filter
group_proband_nonref_reads = function(variants, variantcaller) {
  variants %>% group_by(pool, nonref_alleles_probands) %>% 
    summarise(n_recovered = sum(recovered_reads), total = length(recovered_reads)) -> vars_alleles
  vars_alleles$percent_recovered = vars_alleles$n_recovered/vars_alleles$total*100
  vars_alleles$pool_id = paste0('Pool_', vars_alleles$pool)
  vars_alleles$variantcaller = variantcaller
  return(vars_alleles)
}

# Calculate recall over all variants in all probands
recall_all_probands = function(variants) {
  variants %>% group_by(pool, simulation) %>% 
    summarise(n_recovered = sum(recovered, na.rm = T), total = length(recovered)) -> vars_alleles
  vars_alleles$percent_recovered = vars_alleles$n_recovered/vars_alleles$total*100
  vars_alleles$pool_id = paste0('Pool_', vars_alleles$pool)
  return(vars_alleles)
}

# Calculate recall with each variant locus once (even if multiple alleles and/or in multiple probands)
# Calculate the recovery rate over loci genotyped non-ref in any proband
# but only counting each locus once (even if there were multiple non-ref probands)
recall_by_locus = function(variants) {
  variants %>% group_by(pool, simulation, variant) %>% 
    summarise(any_recovered = any(recovered, na.rm = T)) -> any_recovered
  
  any_recovered %>% group_by(pool, simulation) %>% 
    summarise(n_recovered = sum(any_recovered, na.rm = T), 
              total = length(any_recovered)) -> vars_alleles
  
  vars_alleles$percent_recovered = vars_alleles$n_recovered/vars_alleles$total*100
  vars_alleles$pool_id = paste0('Pool_', vars_alleles$pool)
  return(vars_alleles)
}

# Calculate recall per proband
recall_by_proband = function(variants) {
  variants %>% group_by(pool, simulation, proband) %>% 
    summarise(n_recovered = sum(recovered, na.rm = T), 
              total = length(recovered)) -> vars_alleles
  
  vars_alleles$percent_recovered = vars_alleles$n_recovered/vars_alleles$total*100
  vars_alleles$pool_id = paste0('Pool_', vars_alleles$pool)
  return(vars_alleles)
}

falsepos_rate = function(variants) {
  # Only count each false positive locus (remove multiple entries for additional probands)
  variants = subset(variants, select=-c(proband, recovered))
  dup_rows = duplicated(variants[,c('variant', 'pool', 'simulation')])
  variants = variants[!dup_rows,]
  variants %>% group_by(pool, simulation) %>% 
    summarise(n_fp = sum(falsepos, na.rm = T), total = length(falsepos)) -> vars_alleles
  vars_alleles$percent_fp = vars_alleles$n_fp/vars_alleles$total*100
  return(vars_alleles)
}

# Extract and combine data from multiple csv files output by filter_multiVCF.py
extract_multi_csvs = function(data_dir, replicates, sim_type) {
  start.time <- Sys.time()
  
  all_data_list = list()
  for (simulation in 1:length(replicates)) {
    sim_files = list.files(path = data_dir, pattern = paste0(replicates[simulation],".*csv$"))
    print(sim_files)
    for (f in sim_files) {
      pool = strsplit(f, '.', fixed = T)[[1]][3]
      data = read.csv(paste0(data_dir, f), stringsAsFactors = F)
      data$simulation = simulation
      data$pool = as.integer(pool)
      data$total_reads_pool = NA
      #all_data_list[[simulation]] = data
      all_data_list[[(length(all_data_list) +1)]] = data
    }
  }
  all_data = bind_rows(all_data_list)
  all_data$simulation = as.factor(all_data$simulation)
  all_data$sim_type = sim_type
  colnames(all_data)[colnames(all_data)=="recovered_in_proband"] <- "recovered"
  colnames(all_data)[colnames(all_data)=="nonref_alleles_probands"] <- "nonref_allele_count_truth"
  
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  print(time.taken)
  
  return(all_data)
}

# Extract and combine data from multiple csv files output by filter_individualVCF.py
extract_ind_csvs = function(data_dir, replicates, sim_type) {
  start.time <- Sys.time()
  
  all_data_list = list()
  for (simulation in 1:length(replicates)) {
    sim_files = list.files(path = data_dir, pattern = paste0(replicates[simulation],".*csv$"))
    print(sim_files)
    for (f in sim_files) {
      data = read.csv(paste0(data_dir, f), stringsAsFactors = F)
      data$simulation = simulation
      all_data_list[[(length(all_data_list) +1)]] = data
    }
  }
  all_data = bind_rows(all_data_list)
  all_data$simulation = as.factor(all_data$simulation)
  all_data$sim_type = sim_type
  colnames(all_data)[colnames(all_data)=="recovered_in_proband"] <- "recovered"
  colnames(all_data)[colnames(all_data)=="nonref_alleles_probands"] <- "nonref_allele_count_truth"
  
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  print(time.taken)
  
  return(all_data)
}
