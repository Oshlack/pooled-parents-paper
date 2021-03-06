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
      data$pool = as.integer(strsplit(strsplit(f, '.', fixed = T)[[1]][2], '_')[[1]][2])
      # Calculate nonref_allele_count_truth over all probands
      group_by(data, position) %>% 
        summarise(nonref_allele_count_truth = sum(nonref_alleles_proband, na.rm = T),
                  total_alleles_probands = sum(total_alleles_proband, na.rm = T),
                  nonref_reads_probands = sum(nonref_reads_proband, na.rm = T)
        ) -> summed_alleles
      data = merge(data, summed_alleles, all.x = T)
      all_data_list[[(length(all_data_list) +1)]] = data
    }
  }
  all_data = bind_rows(all_data_list)
  all_data$simulation = as.factor(all_data$simulation)
  all_data$sim_type = sim_type
  all_data$recovered = all_data$recovered_proband
  
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  print(time.taken)
  
  return(all_data)
}

get_vep_header = function(vep_tsv_file) {
  all_lines = readLines(vep_tsv_file)
  header_row = all_lines[grep("#Uploaded_variation", all_lines)]
  header_vect = strsplit(header_row, split='\t')[[1]]
  header_vect = sub('#', '', header_vect)
  return(header_vect)
}

convert_location = function(l) {
  l_split = strsplit(l, ':')
  chrom = l_split[[1]][1]
  pos = as.integer(l_split[[1]][2])
  pos_string = formatC(pos, width = 9, format = "d", flag = "0")
  return(paste0(chrom,'_',pos_string))
}

variant_to_location = function(variant) {
  v_split = strsplit(variant, '_')
  return(paste0(v_split[[1]][c(1,2,4)], collapse = '_'))
}

# Extract vep tsv annotations from multiple individuals
extract_vep_tsvs = function(data_dir) {
  files = list.files(path = data_dir, pattern = ".*vep.tsv$")
  all_data_list = list()
  for(i in 1:length(files)) {
    sample_id = strsplit(files[i], '.', fixed=T)[[1]][1]
    this_file = paste0(data_dir,files[i])
    data = read.table(this_file, stringsAsFactors = F, comment.char = "#", 
                      col.names=get_vep_header(this_file), na.strings = "-")
    data$proband = sample_id
    all_data_list[[(length(all_data_list) +1)]] = data
  }
  all_data = bind_rows(all_data_list)
  all_data$location_str = sapply(all_data$Location, convert_location)
  all_data$location_str = paste(all_data$location_str, all_data$Allele, sep = '_')
  return(all_data)
}

# simplify consequences by taking first if multiple and removing "_variant"
simplify_consequence = function(consequence) {
  first = strsplit(consequence, ',', fixed=T)[[1]][1]
  simplified = sub('_variant', '', first)
  return(simplified)
}

# Set order of factor by frequency
factor_by_freq = function(vect, decreasing=T) {
  factor(vect, levels=names(sort(table(vect), decreasing=decreasing)))
}