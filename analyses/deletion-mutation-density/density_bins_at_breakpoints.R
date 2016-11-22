#!/usr/bin/env Rscript

library(VariantAnnotation)
library(data.table)
library(GenomicFeatures)
library(progress)

'usage: density_bins_at_breakpoints.R -o <x> -b <l> -n <a> -m <y> -d <z> -i <j> -v <w> -c <q> -r <u>

options:
 -o <x> Chromosome
 -b <l> Bin Width
 -n <a> Number of Bins
 -m <y> Donor Meta Path
 -d <z> Deletion Ranges Path
 -i <j> Deletion Info Path
 -v <w> SNV Ranges Path
 -c <q> Carrier Mask Path
 -r <u> Result Path' -> doc

library(docopt)

opts = docopt(doc)
selected_chrom = opts$h
bin_width = as.integer(opts$b)
num_bins = as.integer(opts$n)
donor_meta_path = opts$m
deletion_ranges_path = opts$d
deletion_info_path = opts$i
snv_ranges_path = opts$v
carrier_mask_path = opts$c
result_path = opts$r

load(donor_meta_path)
load(deletion_ranges_path)
load(deletion_info_path)
load(snv_ranges_path)
load(carrier_mask_path)

#Temporary
#selected_chrom = "20"
#bin_width = 100
# Number of bins to each side of breakpoint
#num_bins = 11

deletion_filter = which(deletion_info$IMPRECISE == T | as.character(seqnames(deletion_ranges)) != selected_chrom)
snv_filter = NULL

filtered_deletion_carrier_mask = deletion_carrier_mask

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

breakpoints = c(IRanges(start(filtered_deletion_ranges), width=1), IRanges(end(filtered_deletion_ranges), width=1))

breakpoint_bins = flank(breakpoints, bin_width * num_bins, both=T)

all_tiles = tile(breakpoint_bins, n=num_bins)

hits = which(filtered_deletion_carrier_mask[,] > 0, arr.ind = T)
pb = progress_bar$new(format=":current/:total [:bar] :percent :elapsed :eta",total = dim(hits)[1])

binned_densities = apply(hits, 1, function(x){pb$tick(); countOverlaps(all_tiles[[x[1]]], ranges(filtered_snv_ranges[[x[2]]])) / filtered_snv_counts[x[2]];})

var_name = paste("breakpoint_density_bins_chrom_", selected_chrom, "_width_", bin_width, "_num_", num_bins, sep="")
assign(var_name, binned_densities)
save(list=c(var_name), file=paste(result_path, "/", var_name, ".RData", sep=""))
