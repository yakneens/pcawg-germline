#!/usr/bin/env Rscript

'usage: density_bins_at_breakpoints.R [options]

options:
-o CHROM Chromosome
-s SUB_TYPE  Sub Type
-b BIN_WIDTH  Bin Width
-n NUM_BINS  Number of Bins
-m DONOR_META_PATH  Donor Meta Path
-d DELETION_RANGES_PATH  Deletion Ranges Path
-i DELETION_INFO_PATH  Deletion Info Path
-v SNV_RANGES_PATH  SNV Ranges Path
-c CARRIER_MASK_PATH  Carrier Mask Path
-r RESULT_PATH  Result Path' -> doc

library(docopt)

opts = docopt(doc)
selected_chrom = opts$o
sub_type = opts$s
bin_width = as.integer(opts$b)
num_bins = as.integer(opts$n)
donor_meta_path = opts$m
deletion_ranges_path = opts$d
deletion_info_path = opts$i
snv_ranges_path = opts$v
carrier_mask_path = opts$c
result_path = opts$r

library(VariantAnnotation)
library(data.table)
library(GenomicFeatures)
library(progress)


# donor_meta_path = "~/Downloads/pcawg_data/del_density/by_project/donor_meta.RData"
# deletion_ranges_path = "~/Downloads/pcawg_data/del_density/by_project/deletion_ranges.RData"
# deletion_info_path = "~/Downloads/pcawg_data/del_density/by_project/deletion_info.RData"
# snv_ranges_path = "~/Downloads/pcawg_data/del_density/by_project/snv_ranges.RData"
# carrier_mask_path = "~/Downloads/pcawg_data/del_density/by_project/deletion_carrier_mask.RData"

load(donor_meta_path)
load(deletion_ranges_path)
load(deletion_info_path)
load(snv_ranges_path)
load(carrier_mask_path)

filtered_snv_ranges = list()

if(sub_type == "c_to_a"){
  filtered_snv_ranges = lapply(snv_ranges, function(x) x[which((as.character(x$REF) == "C" & as.character(unlist(x$ALT)) == "A") | (as.character(x$REF) == "G" & as.character(unlist(x$ALT)) == "T"))])
} else if(sub_type == "c_to_g"){
  filtered_snv_ranges = lapply(snv_ranges, function(x) x[which((as.character(x$REF) == "C" & as.character(unlist(x$ALT)) == "G") | (as.character(x$REF) == "G" & as.character(unlist(x$ALT)) == "C"))])
} else if(sub_type == "c_to_t"){
  c_to_t_snv_ranges = lapply(snv_ranges, function(x) x[which((as.character(x$REF) == "C" & as.character(unlist(x$ALT)) == "T") | (as.character(x$REF) == "G" & as.character(unlist(x$ALT)) == "A"))])  
} else if(sub_type == "t_to_g"){
  t_to_g_snv_ranges = lapply(snv_ranges, function(x) x[which((as.character(x$REF) == "T" & as.character(unlist(x$ALT)) == "G") | (as.character(x$REF) == "A" & as.character(unlist(x$ALT)) == "C"))])
} else if(sub_type == "t_to_a"){
  t_to_a_snv_ranges = lapply(snv_ranges, function(x) x[which((as.character(x$REF) == "T" & as.character(unlist(x$ALT)) == "A") | (as.character(x$REF) == "A" & as.character(unlist(x$ALT)) == "T"))])
} else if(sub_type == "t_to_c"){
  t_to_c_snv_ranges = lapply(snv_ranges, function(x) x[which((as.character(x$REF) == "T" & as.character(unlist(x$ALT)) == "C") | (as.character(x$REF) == "A" & as.character(unlist(x$ALT)) == "G"))])
}

binned_densities = list()
snv_counts = unlist(lapply(snv_ranges, length))


print(paste("Processing chromosome:", selected_chrom))

deletion_filter = which(deletion_info$IMPRECISE == T |  as.character(seqnames(deletion_ranges)) != selected_chrom)
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

pre_deletion_ranges = flank(filtered_deletion_ranges, width(filtered_deletion_ranges), start = T)
post_deletion_ranges = flank(filtered_deletion_ranges, width(filtered_deletion_ranges), start = F)

deletion_tiles = tile(filtered_deletion_ranges, n = 10)
pre_deletion_tiles = tile(pre_deletion_ranges, n = 10)
post_deletion_tiles = tile(post_deletion_ranges, n = 10)

all_tiles = pc(pre_deletion_tiles, deletion_tiles, post_deletion_tiles)

hits = which(filtered_deletion_carrier_mask[,] > 0, arr.ind = T)

#Filter out samples that don't have any variants left after filtering
empty_samples = which(lapply(filtered_snv_ranges, length) == 0)

if(length(empty_samples) > 0){
  hits = hits[-which(hits[,2] %in% empty_samples),]
}
pb = progress_bar$new(format=":current/:total [:bar] :percent :elapsed :eta",total = dim(hits)[1])

binned_densities[[selected_chrom]] = apply(hits, 1, function(x){pb$tick(); countOverlaps(all_tiles[[x[1]]], filtered_snv_ranges[[x[2]]]) / filtered_snv_counts[x[2]];})

var_name = paste("density_bins_chrom_", selected_chrom, "_", sub_type,sep="")
assign(var_name, binned_densities)
save(list=c(var_name), file=paste(result_path, "/", var_name, ".RData", sep=""))
