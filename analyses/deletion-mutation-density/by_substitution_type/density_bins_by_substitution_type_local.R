#!/usr/bin/env Rscript

library(VariantAnnotation)
library(data.table)
library(GenomicFeatures)
library(progress)
library(logging)



donor_meta_path = "~/Downloads/pcawg_data/del_density/input_data/donor_meta.RData"
deletion_ranges_path = "~/Downloads/pcawg_data/del_density/input_data/deletion_ranges.RData"
deletion_info_path = "~/Downloads/pcawg_data/del_density/input_data/deletion_info.RData"
snv_ranges_path = "~/Downloads/pcawg_data/del_density/input_data/snv_ranges.RData"
carrier_mask_path = "~/Downloads/pcawg_data/del_density/input_data/deletion_carrier_mask.RData"
non_carriers = FALSE



if(non_carriers){
  carrier_str = "non_carriers"
}else{
  carrier_str = "carriers"
}

load(donor_meta_path)
load(deletion_ranges_path)
load(deletion_info_path)
load(snv_ranges_path)
load(carrier_mask_path)

sub_types = c("c_to_a", "c_to_t", "c_to_g", "t_to_g", "t_to_c", "t_to_a")


sub_type_mapping = list()
sub_type_mapping[["c_to_a"]] = list(c("C", "A"), c("G", "T"))
sub_type_mapping[["c_to_g"]] = list(c("C", "G"), c("G", "C"))
sub_type_mapping[["c_to_t"]] = list(c("C", "T"), c("G", "A"))
sub_type_mapping[["t_to_a"]] = list(c("T", "A"), c("A", "T"))
sub_type_mapping[["t_to_g"]] = list(c("T", "G"), c("A", "C"))
sub_type_mapping[["t_to_c"]] = list(c("T", "C"), c("A", "G"))

#deletion_filter = which(as.character(seqnames(deletion_ranges)) != selected_chrom)
deletion_filter = NULL
snv_filter = NULL

filtered_snv_ranges = snv_ranges

filtered_deletion_carrier_mask = deletion_carrier_mask
filtered_snv_counts = snv_counts
if(!is.null(deletion_filter)){
  filtered_deletion_ranges = deletion_ranges[-deletion_filter]
  filtered_deletion_carrier_mask = deletion_carrier_mask[-deletion_filter,]
}else{
  filtered_deletion_ranges = deletion_ranges
  filtered_deletion_carrier_mask = deletion_carrier_mask
}
loginfo("Filtered deletions. Retaining %d of total %d deletions.", length(filtered_deletion_ranges), length(deletion_ranges))

pre_deletion_ranges = flank(filtered_deletion_ranges, width(filtered_deletion_ranges), start = T)
post_deletion_ranges = flank(filtered_deletion_ranges, width(filtered_deletion_ranges), start = F)

deletion_tiles = tile(filtered_deletion_ranges, n = 10)
pre_deletion_tiles = tile(pre_deletion_ranges, n = 10)
post_deletion_tiles = tile(post_deletion_ranges, n = 10)

all_tiles = pc(pre_deletion_tiles, deletion_tiles, post_deletion_tiles)
loginfo("Created %s density bins",  dim(all_tiles))

if(non_carriers){
  hits = which(filtered_deletion_carrier_mask == 0, arr.ind = T)
}else{
  hits = which(filtered_deletion_carrier_mask > 0, arr.ind = T)
}

empty_samples = which(lapply(filtered_snv_ranges, length) == 0)

if(length(empty_samples) > 0){
  hits = hits[-which(hits[,2] %in% empty_samples),]
  loginfo("Filtering out samples that do not have substitutions of %s type. Filtered out samples: %s", sub_type, empty_samples)
}

dels_by_donor = as.data.table(hits)[, list(lapply(.SD, c)), by=col]

get_overlap_counts<- function(donor_index, del_list, all_tiles){
  tiles_by_donor = unlist(all_tiles[unlist(del_list)])
  density_bins_by_sub_type = list()
  for(sub_type in sub_types){
    variants_by_donor = filtered_snv_ranges[[donor_index]]
    sub_type_variants_by_donor = variants_by_donor[which((as.character(variants_by_donor$REF) == sub_type_mapping[[sub_type]][[1]][1] & 
                                                          as.character(unlist(variants_by_donor$ALT)) == sub_type_mapping[[sub_type]][[1]][2]) | 
                                                         (as.character(variants_by_donor$REF) == sub_type_mapping[[sub_type]][[2]][1] & 
                                                          as.character(unlist(variants_by_donor$ALT)) == sub_type_mapping[[sub_type]][[2]][2]))]
    overlaps = as.data.table(findOverlaps(tiles_by_donor, sub_type_variants_by_donor))
    overlap_counts = overlaps[, .N, queryHits]
    x = vector(mode="integer", length = 30 * length(del_list))
    names(x) = unlist(lapply(names(deletion_ranges)[del_list], function(x){rep(x,30)}))
    x[overlap_counts$queryHits] = overlap_counts$N / length(sub_type_variants_by_donor) / width(deletion_ranges[1 + overlap_counts$queryHits %/% 30])
    density_bins_by_sub_type[[sub_type]] = matrix(x, nrow=30)  
  }
  return(density_bins_by_sub_type)
}

pb = progress_bar$new(format=":current/:total [:bar] :percent :elapsed :eta",total = dim(dels_by_donor)[1])
binned_densities = apply(dels_by_donor, 1, function(x){pb$tick(); get_overlap_counts(x[[1]], x[[2]], all_tiles)}) %>% pmap(cbind)



