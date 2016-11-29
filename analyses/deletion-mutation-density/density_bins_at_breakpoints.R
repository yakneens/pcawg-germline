#!/usr/bin/env Rscript

library(VariantAnnotation)
library(data.table)
library(GenomicFeatures)
library(progress)
library(logging)

'usage: density_bins_at_breakpoints.R [options]

options:
 --chrom CHROM  Chromosome
 --bin_width BIN_WIDTH  Bin Width
 --num_bins NUM_BINS  Number of Bins
 --donor_meta_path DONOR_META_PATH  Donor Meta Path
 --deletion_ranges_path DELETION_RANGES_PATH  Deletion Ranges Path
 --deletion_info_path DELETION_INFO_PATH  Deletion Info Path
 --snv_ranges_path SNV_RANGES_PATH SNV  Ranges Path
 --carrier_mask_path CARRIER_MASK_PATH  Carrier Mask Path
 --result_path RESULT_PATH  Result Path
 --breakpoint_margin BREAKPOINT_MARGIN  Breakpoint Margin
 --non_carriers NON_CARRIERS  Non Carriers [default: FALSE]' -> doc

library(docopt)

basicConfig()


opts = docopt(doc)
selected_chrom = opts$chrom
bin_width = as.integer(opts$bin_width)
num_bins = as.integer(opts$num_bins)
breakpoint_margin = as.integer(opts$breakpoint_margin)
donor_meta_path = opts$donor_meta_path
deletion_ranges_path = opts$deletion_ranges_path
deletion_info_path = opts$deletion_info_path
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
load(deletion_info_path)
load(snv_ranges_path)
load(carrier_mask_path)

loginfo("Loaded input data")


del_widths = width(deletion_ranges)


bin_based_del_size_cutoff = (bin_width * num_bins) / 2 + breakpoint_margin
loginfo("Set deletion size cutoff at: %d", bin_based_del_size_cutoff)


deletion_filter = which(as.character(seqnames(deletion_ranges)) != selected_chrom | del_widths < bin_based_del_size_cutoff)

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

start_breakpoints = GRanges(IRanges(start(filtered_deletion_ranges), width=1), seqnames = as.integer(selected_chrom))
start_breakpoint_bins = flank(start_breakpoints, bin_width * num_bins, both=T)
start_tiles = tile(start_breakpoint_bins, n=num_bins)

end_breakpoints = GRanges(IRanges(end(filtered_deletion_ranges), width=1), seqnames = as.integer(selected_chrom))
end_breakpoint_bins = flank(end_breakpoints, bin_width * num_bins, both=T)
end_tiles = tile(end_breakpoint_bins, n=num_bins)

loginfo("Created breakpoint tiles. %d start tiles and %d end tiles", length(start_tiles), length(end_tiles))

if(non_carriers){
  hits = which(filtered_deletion_carrier_mask == 0, arr.ind = T)
}else{
  hits = which(filtered_deletion_carrier_mask > 0, arr.ind = T)
}
loginfo("Analyzing %s hits", dim(hits))

pb = progress_bar$new(format=":current/:total [:bar] :percent :elapsed :eta",total = dim(hits)[1])

start_binned_densities = apply(hits, 1, function(x){pb$tick(); countOverlaps(start_tiles[[x[1]]], filtered_snv_ranges[[x[2]]]) / filtered_snv_counts[x[2]];})
loginfo("Completed %s start densities", dim(start_binned_densities))

end_binned_densities = apply(hits, 1, function(x){countOverlaps(end_tiles[[x[1]]], filtered_snv_ranges[[x[2]]]) / filtered_snv_counts[x[2]];})
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
