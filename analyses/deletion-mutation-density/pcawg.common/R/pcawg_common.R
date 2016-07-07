library(VariantAnnotation)
library(splitstackshape)


#Read in sample metadata
get_pcawg_metadata <- function(sample_file_location, blacklist_file_location, split_multi_tumors = T){
  sample_meta = read.table(sample_file_location, header=TRUE, sep="\t", stringsAsFactors = F)
  blacklist_meta = read.table(blacklist_file_location, header=TRUE, sep="\t", stringsAsFactors = F)
  
  sample_meta = sample_meta[-na.omit(match(blacklist_meta$donor_unique_id, sample_meta$donor_unique_id)),]
  
  if(split_multi_tumors){
    library(splitstackshape)
    
    #Split samples with multiple tumors into separate rows
    sample_meta = cSplit(sample_meta, "tumor_wgs_aliquot_id", direction="long")
  }
  
  return(sample_meta)
}

keep_first_of_multi_tumors <- function(sample_metadata){
  sample_metadata$tumor_wgs_aliquot_id = unlist(lapply(strsplit(sample_metadata$tumor_wgs_aliquot_id, ","), function(x) x[1]))
  return(sample_metadata)
}

#Read SNV samples from a specified directory
read_snv_samples_from_vcf <- function(sample_location_directory, sample_meta, range){
  
  files = list.files(sample_location_directory, pattern="*.vcf.gz$", full.names=T, recursive=F)
  
  snv_samples = list()
  
  for(cur_file in files){
    #Extract tumor_wgs_aliquot_id from filename
    cur_tumor_uuid = gsub(".*([0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}).*", "\\1", x=cur_file)
    
    if(!missing(range)){
      tab = TabixFile(cur_file)
      snv_samples[cur_tumor_uuid] = readVcf(tab, "hs37d5", param=range)
    }else{
      snv_samples[cur_tumor_uuid] = readVcf(cur_file, genome="hs37d5")
    }
  }
  
  return(snv_samples)
}

#Save SNV samples in R format to a specified location
save_snv_samples <- function(snv_samples, file_location){
  save(snv_samples, file=file_location)
}

#When deletions are imported they don't seem to have a proper end recorded
set_deletion_range_ends <- function(dels){
  eval.parent(substitute(end(ranges(rowRanges(dels))) <- info(dels)$END))
}

#Read in germline deletions
read_germline_deletions_from_vcf <- function(file_location, range){
  if(!missing(range)){
    tab = TabixFile(file_location)
    dels = readVcf(tab, "hs37d5", param=range)
  }else{
    dels = readVcf(file_location, genome="hs37d5")
  }
  
  set_deletion_range_ends(dels)
  return(dels)
}

#Read in germline deletions
read_germline_snvs_from_vcf <- function(file_location, range){
  if(!missing(range)){
    tab = TabixFile(file_location)
    germ_snv = readVcf(tab, "hs37d5", param=range)
  }else{
    germ_snv = readVcf(file_location, genome="hs37d5")
  }
  
  return(germ_snv)
}

#Save deletions in R format to a specified location
save_deletions_to_r_file <- function(deletions, file_location){
  save(deletions, file=file_location)
}

#If genotype of element is in genotypes_list record it as 1, otherwise record it as 0
genotype_mask <- function(element, genotypes_list){
  if(element %in% genotypes_list){
    return(1)
  } else {
    return(0)
  }
}

#If variant genotype is "0/1" or "1/0" mark it 1
#If variant genotype is "0/0" mark it 0
#Otherwise mark NA
mark_carriers_and_non_carriers <- function(element){
  if(element == "0/1" || element == "1/0"){
    return(1)
  } else if (element == "0/0") {
    return(0)
  } else {
    return(NA)
  }
}

get_het_carrier_mask <- function(my_data){
  return(as.data.frame(apply(my_data, MARGIN=c(1,2), FUN=mark_carriers_and_non_carriers)))
}

get_clinical_metadata <- function(file_location, sample_metadata){
  clinical_metadata = read.csv(file_location, sep = "\t", stringsAsFactors = F)
  colnames(clinical_metadata)[1] = "donor_unique_id"
  clinical_metadata = clinical_metadata[match(sample_metadata$donor_unique_id, clinical_metadata$donor_unique_id),]
  clinical_metadata[clinical_metadata$donor_sex == "" & !is.na(clinical_metadata$donor_sex),]$donor_sex = NA
  clinical_metadata$donor_unique_id = sample_metadata$donor_unique_id
  
  to_remap = cbind(c("25", "25.1", "25.0", "25.2", "25.0/25.9", "25.7", "25.1/25.2", "c92.0", "C92.00", "C26.8"), c("C25", "C25.1", "C25.0", "C25.2", "C25.0", "C25.7", "C25.1", "C92.0", "C92.0", "C26.9"))
  apply(to_remap, 1, function(x) clinical_metadata[which(clinical_metadata$donor_diagnosis_icd10 == x[1]),]$donor_diagnosis_icd10 = x[2])
  
  
  return(clinical_metadata)
}

get_histology_metadata <- function(file_location, sample_metadata){
  hist_meta = read.csv(file_location, sep = "\t", stringsAsFactors = F)
  hist_meta = hist_meta[match(sample_metadata$donor_unique_id, hist_meta$donor_unique_id),]
  
  return(hist_meta)
  
}