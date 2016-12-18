#!/usr/bin/env Rscript

library(VariantAnnotation)
library(data.table)
library(GenomicFeatures)
library(progress)
library(logging)

selected_chrom = 20
donor_meta_path = "~/Downloads/pcawg_data/del_density/input_data/donor_meta.RData"
deletion_ranges_path = "~/Downloads/pcawg_data/del_density/input_data/deletion_ranges.RData"
snv_ranges_path = "~/Downloads/pcawg_data/del_density/input_data/snv_ranges.RData"
carrier_mask_path = "~/Downloads/pcawg_data/del_density/input_data/deletion_carrier_mask.RData"
result_path = "~/Downloads/pcawg_data/del_density/"
non_carriers = FALSE



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


#deletion_filter = which(as.character(seqnames(deletion_ranges)) != selected_chrom)
num_bins = 10

deletion_filter = which(width(deletion_ranges) < 1400)
deletion_filter = NULL
snv_filter = NULL

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

dels_by_donor = as.data.table(hits)[, list(lapply(.SD, c)), by=col]

get_overlap_counts<- function(donor_index, del_list, all_tiles){
  tiles_by_donor = unlist(all_tiles[unlist(del_list)])
  variants_by_donor = filtered_snv_ranges[[donor_index]]
  overlaps = as.data.table(findOverlaps(tiles_by_donor, variants_by_donor))
  overlap_counts = overlaps[, .N, queryHits]
  x = vector(mode="integer", length = 30 * length(del_list))
  names(x) = unlist(lapply(names(deletion_ranges)[del_list], function(x){rep(x,30)}))
  x[overlap_counts$queryHits] = overlap_counts$N / length(variants_by_donor) / width(deletion_ranges[1 + overlap_counts$queryHits %/% 30])
  return(x)
}

pb = progress_bar$new(format=":current/:total [:bar] :percent :elapsed :eta",total = dim(dels_by_donor)[1])


binned_densities = do.call(cbind, lapply(apply(dels_by_donor, 1, function(x){pb$tick(); get_overlap_counts(x[[1]], x[[2]], all_tiles)}), function(x){matrix(x, nrow=30);}))

loginfo("Completed %s density bins", dim(binned_densities))

#var_name = paste("density_bins_", selected_chrom, "_", carrier_str, sep="")
#assign(var_name, get_binned_densities(which(as.character(seqnames(deletion_ranges)) != selected_chrom), NULL))
#full_path = paste(result_path, "/", var_name, ".RData", sep="")
#save(list=c(var_name), file=full_path)
#loginfo("Saved density bins to %s", full_path)

