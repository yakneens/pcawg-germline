#!/usr/bin/env Rscript

'usage: density_bins_at_breakpoints.R [options]

options:
--chrom CHROM  Chromosome
--non_carriers NON_CARRIERS  Non Carriers [default: FALSE]
--sub_type SUB_TYPE  Sub Type
--bin_width BIN_WIDTH  Bin Width
--num_bins NUM_BINS  Number of Bins
--donor_meta_path DONOR_META_PATH  Donor Meta Path
--deletion_ranges_path DELETION_RANGES_PATH  Deletion Ranges Path
--deletion_info_path DELETION_INFO_PATH  Deletion Info Path
--snv_ranges_path SNV_RANGES_PATH  SNV Ranges Path
--carrier_mask_path CARRIER_MASK_PATH  Carrier Mask Path
--result_path RESULT_PATH  Result Path' -> doc

library(docopt)
library(logging)

opts = docopt(doc)
selected_chrom = opts$chrom
non_carriers = as.logical(opts$non_carriers)
sub_type = opts$sub_type
bin_width = as.integer(opts$bin_width)
num_bins = as.integer(opts$num_bins)
donor_meta_path = opts$donor_meta_path
deletion_ranges_path = opts$deletion_ranges_path
deletion_info_path = opts$deletion_info_path
snv_ranges_path = opts$snv_ranges_path
carrier_mask_path = opts$carrier_mask_path
result_path = opts$result_path

loginfo("Parsed command-line options: %s", opts)

if(non_carriers){
  carrier_str = "non_carriers"
}else{
  carrier_str = "carriers"
}

library(VariantAnnotation)
library(data.table)
library(GenomicFeatures)
library(progress)

load(donor_meta_path)
load(deletion_ranges_path)
load(deletion_info_path)
load(snv_ranges_path)
load(carrier_mask_path)

loginfo("Loaded input data")

filtered_snv_ranges = list()

if(sub_type == "c_to_a"){
  filtered_snv_ranges = lapply(snv_ranges, function(x) x[which((as.character(x$REF) == "C" & as.character(unlist(x$ALT)) == "A") | (as.character(x$REF) == "G" & as.character(unlist(x$ALT)) == "T"))])
} else if(sub_type == "c_to_g"){
  filtered_snv_ranges = lapply(snv_ranges, function(x) x[which((as.character(x$REF) == "C" & as.character(unlist(x$ALT)) == "G") | (as.character(x$REF) == "G" & as.character(unlist(x$ALT)) == "C"))])
} else if(sub_type == "c_to_t"){
  filtered_snv_ranges = lapply(snv_ranges, function(x) x[which((as.character(x$REF) == "C" & as.character(unlist(x$ALT)) == "T") | (as.character(x$REF) == "G" & as.character(unlist(x$ALT)) == "A"))])  
} else if(sub_type == "t_to_g"){
  filtered_snv_ranges = lapply(snv_ranges, function(x) x[which((as.character(x$REF) == "T" & as.character(unlist(x$ALT)) == "G") | (as.character(x$REF) == "A" & as.character(unlist(x$ALT)) == "C"))])
} else if(sub_type == "t_to_a"){
  filtered_snv_ranges = lapply(snv_ranges, function(x) x[which((as.character(x$REF) == "T" & as.character(unlist(x$ALT)) == "A") | (as.character(x$REF) == "A" & as.character(unlist(x$ALT)) == "T"))])
} else if(sub_type == "t_to_c"){
  filtered_snv_ranges = lapply(snv_ranges, function(x) x[which((as.character(x$REF) == "T" & as.character(unlist(x$ALT)) == "C") | (as.character(x$REF) == "A" & as.character(unlist(x$ALT)) == "G"))])
}
loginfo("Filtered SNV Ranges to selected substitution type %s", sub_type)

binned_densities = list()
snv_counts = unlist(lapply(snv_ranges, length))

deletion_filter = which(as.character(seqnames(deletion_ranges)) != selected_chrom)
snv_filter = NULL

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

hits = which(filtered_deletion_carrier_mask[,] > 0, arr.ind = T)
loginfo("Analyzing %s hits", dim(hits))

#Filter out samples that don't have any variants left after filtering
empty_samples = which(lapply(filtered_snv_ranges, length) == 0)

if(length(empty_samples) > 0){
  hits = hits[-which(hits[,2] %in% empty_samples),]
  loginfo("Filtering out samples that do not have substitutions of %s type. Filtered out samples: %s", sub_type, empty_samples)
}
pb = progress_bar$new(format=":current/:total [:bar] :percent :elapsed :eta",total = dim(hits)[1])

binned_densities[[selected_chrom]] = apply(hits, 1, function(x){pb$tick(); countOverlaps(all_tiles[[x[1]]], filtered_snv_ranges[[x[2]]]) / filtered_snv_counts[x[2]];})

loginfo("Completed %s density bins.", dim(binned_densities[[selected_chrom]]))

var_name = paste("density_bins_chrom_", selected_chrom, "_", sub_type, "_", carrier_str, "_all", sep="")
assign(var_name, binned_densities)
full_path = paste(result_path, "/", var_name, ".RData", sep="")
save(list=c(var_name), file=full_path)
loginfo("Saved density bins to %s.", full_path)
