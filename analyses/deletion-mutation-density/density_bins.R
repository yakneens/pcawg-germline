#!/usr/bin/env Rscript

library(VariantAnnotation)
library(data.table)
library(GenomicFeatures)
library(progress)
library(logging)

'usage: density_bins_by_project.R [options]

options:
 --chrom CHROM  Chromosome
 --non_carriers NON_CARRIERS  Non Carriers [default: FALSE]
 --num_bins NUM_BINS  Number of bins per deletion (Total number of bins will be 3 * NUM_BINS)
 --donor_meta_path DONOR_META_PATH Donor  Meta Path
 --deletion_ranges_path DELETION_RANGES_PATH  Deletion Ranges Path
 --snv_ranges_path SNV_RANGES_PATH SNV  Ranges Path
 --carrier_mask_path CARRIER_MASK_PATH  Carrier Mask Path
 --result_path RESULT_PATH  Result Path' -> doc
library(docopt)

basicConfig()

opts = docopt(doc)
selected_chrom = opts$chrom
num_bins = opts$num_bins
donor_meta_path = opts$donor_meta_path
deletion_ranges_path = opts$deletion_ranges_path
snv_ranges_path = opts$snv_ranges_path
carrier_mask_path = opts$carrier_mask_path
result_path = opts$result_path
non_carriers = as.logical(opts$non_carriers)

loginfo("Parsed command-line options: %s", opts)

if(non_carriers){
  carrier_str = "non_carriers"
}else{
  carrier_str = "carriers"
}


load(donor_meta_path)
load(deletion_ranges_path)
load(snv_ranges_path)
load(carrier_mask_path)

loginfo("Loaded input data")

get_binned_densities <- function(deletion_filter, snv_filter){
  snv_counts = unlist(lapply(snv_ranges, length))
  
  filtered_deletion_carrier_mask = deletion_carrier_mask
  filtered_snv_counts = snv_counts
  
  if(!is.null(deletion_filter)){
    filtered_deletion_ranges = deletion_ranges[-deletion_filter]
    filtered_deletion_carrier_mask = filtered_deletion_carrier_mask[-deletion_filter,]
  }else{
    filtered_deletion_ranges = deletion_ranges
    filtered_deletion_carrier_mask = deletion_carrier_mask
  }
  loginfo("Filtered deletions. Retaining %d of total %d deletions.", length(filtered_deletion_ranges), length(deletion_ranges))
  
  if(!is.null(snv_filter)){
    filtered_snv_ranges = snv_ranges[-snv_filter]
    filtered_snv_counts = snv_counts[-snv_filter]
    filtered_deletion_carrier_mask = filtered_deletion_carrier_mask[,-snv_filter]
  }else{
    filtered_snv_ranges = snv_ranges
    filtered_snv_counts = snv_counts
  }
  loginfo("Filtered SNV samples. Retaining %d of total %d samples", length(filtered_snv_ranges), length(snv_ranges))
  
  pre_deletion_ranges = trim(flank(filtered_deletion_ranges, width(filtered_deletion_ranges), start = T))
  post_deletion_ranges = trim(flank(filtered_deletion_ranges, width(filtered_deletion_ranges), start = F))
  
  deletion_tiles = tile(filtered_deletion_ranges, n = as.integer(num_bins))
  pre_deletion_tiles = tile(pre_deletion_ranges, n = as.integer(num_bins))
  post_deletion_tiles = tile(post_deletion_ranges, n = as.integer(num_bins))
  
  all_tiles = pc(pre_deletion_tiles, deletion_tiles, post_deletion_tiles)
  
  loginfo("Created %s  density bins.", dim(all_tiles))
  
  if(non_carriers){
    hits = which(filtered_deletion_carrier_mask == 0, arr.ind = T)
  }else{
    hits = which(filtered_deletion_carrier_mask > 0, arr.ind = T)
  }
 
  loginfo("Analyzing %s hits", dim(hits))
  
  pb = progress_bar$new(format=":current/:total [:bar] :percent :elapsed :eta",total = dim(hits)[1])
  
  binned_densities = apply(hits, 1, function(x){pb$tick(); countOverlaps(all_tiles[[x[1]]], filtered_snv_ranges[[x[2]]]) / filtered_snv_counts[x[2]];})
  loginfo("Completed %s density bins", dim(binned_densities))
  
  return(binned_densities)
  
}

var_name = paste("density_bins_", selected_chrom, "_", carrier_str, sep="")
assign(var_name, get_binned_densities(which(as.character(seqnames(deletion_ranges)) != selected_chrom), NULL))
full_path = paste(result_path, "/", var_name, ".RData", sep="")
save(list=c(var_name), file=full_path)
loginfo("Saved density bins to %s", full_path)

