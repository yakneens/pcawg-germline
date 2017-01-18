#!/usr/bin/env Rscript

library(VariantAnnotation)
library(data.table)
library(GenomicFeatures)
library(progress)


bin_width = 100
num_bins = 30
donor_meta_path = "~/Downloads/pcawg_data/del_density/input_data/donor_meta.RData"
deletion_ranges_path = "~/Downloads/pcawg_data/del_density/input_data/deletion_ranges.RData"
deletion_info_path = "~/Downloads/pcawg_data/del_density/input_data/deletion_info.RData"
snv_ranges_path = "~/Downloads/pcawg_data/del_density/input_data/snv_ranges.RData"
carrier_mask_path = "~/Downloads/pcawg_data/del_density/input_data/deletion_carrier_mask.RData"
breakpoint_margin = 200

load(donor_meta_path)
load(deletion_ranges_path)
load(deletion_info_path)
load(snv_ranges_path)
load(carrier_mask_path)

if(non_carriers){
  carrier_str = "non_carriers"
}else{
  carrier_str = "carriers"
}

loginfo("Loaded input data")

del_widths = width(deletion_ranges)

bin_based_del_size_cutoff = (bin_width * num_bins) / 2 + breakpoint_margin
loginfo("Set deletion size cutoff at: %d", bin_based_del_size_cutoff)

deletion_filter = which(del_widths < bin_based_del_size_cutoff)

snv_filter = NULL

filtered_deletion_carrier_mask = deletion_carrier_mask
snv_counts = unlist(lapply(snv_ranges, length))
filtered_snv_counts = snv_counts

if(!is.null(deletion_filter)){
  filtered_deletion_ranges = deletion_ranges[-deletion_filter]
  filtered_deletion_carrier_mask = deletion_carrier_mask[-deletion_filter,]
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

start_breakpoints = GRanges(IRanges(start(filtered_deletion_ranges), width=1), seqnames = seqnames(filtered_deletion_ranges))
start_breakpoint_bins = flank(start_breakpoints, bin_width * num_bins, both=T)
start_tiles = tile(start_breakpoint_bins, n=num_bins)

end_breakpoints = GRanges(IRanges(end(filtered_deletion_ranges), width=1), seqnames = seqnames(filtered_deletion_ranges))
end_breakpoint_bins = flank(end_breakpoints, bin_width * num_bins, both=T)
end_tiles = tile(end_breakpoint_bins, n=num_bins)

all_tiles = list()
all_tiles[["start"]] = start_tiles
all_tiles[["end"]] = end_tiles

loginfo("Created breakpoint tiles. %d start tiles and %d end tiles", length(start_tiles), length(end_tiles))

if(non_carriers){
  hits = which(filtered_deletion_carrier_mask == 0, arr.ind = T)
}else{
  hits = which(filtered_deletion_carrier_mask > 0, arr.ind = T)
}
loginfo("Analyzing %s hits", dim(hits))

dels_by_donor = as.data.table(hits)[, list(lapply(.SD, c)), by=col]

get_overlap_counts<- function(donor_index, del_list, all_tiles){
  
  density_bins_by_breakpoint_type = list()
  breakpoint_types = c("start", "end")
  for(breakpoint_type in breakpoint_types){
    tiles_by_donor = unlist(all_tiles[[breakpoint_type]][unlist(del_list)])
    variants_by_donor = filtered_snv_ranges[[donor_index]]
    overlaps = as.data.table(findOverlaps(tiles_by_donor, variants_by_donor))
    overlap_counts = overlaps[, .N, queryHits]
    x = vector(mode="integer", length = 30 * length(del_list))
    names(x) = unlist(lapply(names(deletion_ranges)[del_list], function(x){rep(x,30)}))
    x[overlap_counts$queryHits] = overlap_counts$N / length(variants_by_donor)
    density_bins_by_breakpoint_type[[breakpoint_type]] = matrix(x, nrow=30, dimnames = list(NULL, names(x)[1 + (seq(length(x)/30) - 1)*30]))  
  }
  return(density_bins_by_breakpoint_type)
}

pb = progress_bar$new(format=":current/:total [:bar] :percent :elapsed :eta",total = dim(dels_by_donor)[1])
binned_densities = apply(dels_by_donor, 1, function(x){pb$tick(); get_overlap_counts(x[[1]], x[[2]], all_tiles)}) %>% pmap(cbind)


pb = progress_bar$new(format=":current/:total [:bar] :percent :elapsed :eta",total = dim(hits)[1])

start_binned_densities = apply(hits, 1, function(x){pb$tick(); countOverlaps(start_tiles[[x[1]]], filtered_snv_ranges[[x[2]]]) / filtered_snv_counts[x[2]];})
#start_binned_densities = apply(hits, 1, function(x){pb$tick(); countOverlaps(start_tiles[[x[1]]], filtered_snv_ranges[[x[2]]]);})
loginfo("Completed %s start densities", dim(start_binned_densities))

#end_binned_densities = apply(hits, 1, function(x){countOverlaps(end_tiles[[x[1]]], filtered_snv_ranges[[x[2]]]) / filtered_snv_counts[x[2]];})
end_binned_densities = apply(hits, 1, function(x){countOverlaps(end_tiles[[x[1]]], filtered_snv_ranges[[x[2]]]);})
loginfo("Completed %s end densities.", dim(start_binned_densities))

base_var_name = paste("breakpoint_density_bins_chrom_", selected_chrom, "_width_", bin_width, "_num_", num_bins, "_margin_", breakpoint_margin, "_", carrier_str, sep="")

start_var_name = paste("start_", base_var_name, sep="")

assign(start_var_name, start_binned_densities)
full_path = paste(result_path, "/", start_var_name, ".RData", sep="")
save(list=c(start_var_name), file=full_path)
loginfo("Saved start densities to %s.", full_path)

end_var_name = paste("end_", base_var_name, sep="")
assign(end_var_name, end_binned_densities)
full_path = paste(result_path, "/", end_var_name, ".RData", sep="")
save(list=c(end_var_name), file=full_path)
loginfo("Saved end densities to %s", full_path)
