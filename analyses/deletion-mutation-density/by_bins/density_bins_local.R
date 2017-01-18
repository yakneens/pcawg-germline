#!/usr/bin/env Rscript

library(VariantAnnotation)
library(data.table)
library(GenomicFeatures)
library(progress)
library(logging)
library(purrr)

donor_meta_path = "~/Downloads/pcawg_data/del_density/input_data/donor_meta.RData"
deletion_ranges_path = "~/Downloads/pcawg_data/del_density/input_data/deletion_ranges.RData"
snv_ranges_path = "~/Downloads/pcawg_data/del_density/input_data/snv_ranges.RData"
carrier_mask_path = "~/Downloads/pcawg_data/del_density/input_data/deletion_carrier_mask.RData"
result_path = "~/Downloads/pcawg_data/del_density/"
non_carriers = FALSE

skip_breakpoints = T
breakpoint_margin = 1000
scaling_factor = 1e09


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

exclusion_list = c("DEL00185014","DEL00185540","DEL00185702","DEL00185656","DEL00185183","DEL00185635",
                   "DEL00184932","DEL00184909","DEL00185063","DEL00185776","DEL00184729","DEL00184485",
                   "DEL00184528","DEL00185533","DEL00184807","DEL00105731")

deletion_filter = which(width(deletion_ranges) < 3000 | names(deletion_ranges) %in% exclusion_list)
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



if(skip_breakpoints){
  filtered_deletion_ranges = narrow(filtered_deletion_ranges, start=breakpoint_margin + 1, end=-(breakpoint_margin + 1))
  pre_deletion_ranges = filtered_deletion_ranges %>% 
    flank(., width(.), start = T) %>% 
    trim(.) %>% 
    GenomicRanges::shift(., -2*breakpoint_margin)
  post_deletion_ranges = filtered_deletion_ranges %>% 
    flank(., width(.), start = F) %>% 
    trim(.) %>% 
    GenomicRanges::shift(., -2*breakpoint_margin)
}else{
  pre_deletion_ranges = filtered_deletion_ranges %>%
    flank(., width(.), start = T) %>%
    trim(.)
  post_deletion_ranges = filtered_deletion_ranges %>%
    flank(., width(.), start = F) %>%
    trim(.)
}


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
  names(x) = unlist(lapply(names(filtered_deletion_ranges)[del_list], function(x){rep(x,30)}))
  x[overlap_counts$queryHits] = overlap_counts$N / length(variants_by_donor) / width(filtered_deletion_ranges[1 + overlap_counts$queryHits %/% 30]) / length(filtered_deletion_ranges) * scaling_factor
  return(x)
}

pb = progress_bar$new(format=":current/:total [:bar] :percent :elapsed :eta",total = dim(dels_by_donor)[1])


binned_densities = do.call(cbind, lapply(apply(dels_by_donor, 1, function(x){pb$tick(); get_overlap_counts(x[[1]], x[[2]], all_tiles)}), function(x){matrix(x, nrow=30, dimnames = list(NULL, names(x)[1 + (seq(length(x)/30) - 1)*30]));}))

loginfo("Completed %s density bins", dim(binned_densities))



