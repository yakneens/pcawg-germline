#!/usr/bin/env Rscript

library(VariantAnnotation)
library(data.table)
library(GenomicFeatures)


'usage: density_bins_by_project.R [-p -m -d -v -c -r]

options:
 -p Project Name
 -m Donor Meta Path
 -d Deletion Ranges Path
 -v SNV Ranges Path
 -c Carrier Mask Path
 -r Result Path' -> doc
library(docopt)

opts = docopt(doc)
project_name = opts$p
donor_meta_path = opts$m
deletion_ranges_path = opts$d
snv_ranges_path = opts$v
carrier_mask_path = opts$c
result_path = opts$r

load(donor_meta_path)
load(deletion_ranges_path)
load(snv_ranges_path)
load(carrier_mask_path)


get_binned_densities <- function(deletion_filter, snv_filter){
snv_counts = unlist(lapply(snv_ranges, length))
  
  if(!is.null(deletion_filter)){
    filtered_deletion_ranges = deletion_ranges[-deletion_filter]
    filtered_deletion_carrier_mask = deletion_carrier_mask[-deletion_filter,]
  }else{
    filtered_deletion_ranges = deletion_ranges
    filtered_deletion_carrier_mask = deletion_carrier_mask
  }
  
  if(!is.null(snv_filter)){
    filtered_snv_ranges = snv_ranges[-snv_filter]
    filtered_snv_counts = snv_counts[-snv_filter]
    filtered_deletion_carrier_mask = filtered_deletion_carrier_mask[,-snv_filter]
  }else{
    filtered_snv_ranges = snv_ranges
    filtered_snv_counts = snv_counts
  }
  
  pre_deletion_ranges = trim(flank(filtered_deletion_ranges, width(filtered_deletion_ranges), start = T))
  post_deletion_ranges = trim(flank(filtered_deletion_ranges, width(filtered_deletion_ranges), start = F))
  
  deletion_tiles = tile(filtered_deletion_ranges, n = 10)
  pre_deletion_tiles = tile(pre_deletion_ranges, n = 10)
  post_deletion_tiles = tile(post_deletion_ranges, n = 10)
  
  all_tiles = pc(pre_deletion_tiles, deletion_tiles, post_deletion_tiles)
  
  hits = which(filtered_deletion_carrier_mask[,] > 0, arr.ind = T)
  pb = progress_bar$new(format=":current/:total [:bar] :percent :elapsed :eta",total = dim(hits)[1])
  
  binned_densities = apply(hits, 1, function(x){pb$tick(); countOverlaps(all_tiles[[x[1]]], filtered_snv_ranges[[x[2]]]) / filtered_snv_counts[x[2]];})
  
  return(binned_densities)
  
}

var_name = paste("density_bins_", project_name)
assign(var_name, get_binned_densities(NULL, which(donor_meta$dcc_project_code != project_name)))
save(get(var_name), file=paste(result_path, "/", var_name, ".RData"))

