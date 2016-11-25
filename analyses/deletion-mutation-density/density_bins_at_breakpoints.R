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
selected_chrom = opts$o
bin_width = as.integer(opts$b)
num_bins = as.integer(opts$n)
donor_meta_path = opts$m
deletion_ranges_path = opts$d
deletion_info_path = opts$i
snv_ranges_path = opts$v
carrier_mask_path = opts$c
result_path = opts$r

#donor_meta_path = "~/Downloads/pcawg_data/del_density/by_project/donor_meta.RData"
#deletion_ranges_path = "~/Downloads/pcawg_data/del_density/by_project/deletion_ranges.RData"
#deletion_info_path = "~/Downloads/pcawg_data/del_density/by_project/deletion_info.RData"
#snv_ranges_path = "~/Downloads/pcawg_data/del_density/by_project/snv_ranges.RData"
#carrier_mask_path = "~/Downloads/pcawg_data/del_density/by_project/deletion_carrier_mask.RData"

load(donor_meta_path)
load(deletion_ranges_path)
load(deletion_info_path)
load(snv_ranges_path)
load(carrier_mask_path)

#Temporary
#selected_chrom = "20"
#bin_width = 100
# Number of bins to each side of breakpoint
#num_bins = 31

del_widths = width(deletion_ranges)

breakpoint_margin = 200
bin_based_del_size_cutoff = (bin_width * num_bins) / 2 + breakpoint_margin
print(paste("Deletion size cutoff:", bin_based_del_size_cutoff))


deletion_filter = which(deletion_info$IMPRECISE == T | as.character(seqnames(deletion_ranges)) != selected_chrom | del_widths < bin_based_del_size_cutoff)

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

if(!is.null(snv_filter)){
  filtered_snv_ranges = snv_ranges[-snv_filter]
  filtered_snv_counts = snv_counts[-snv_filter]
  filtered_deletion_carrier_mask = filtered_deletion_carrier_mask[,-snv_filter]
}else{
  filtered_snv_ranges = snv_ranges
  filtered_snv_counts = snv_counts
}

print(paste(length(filtered_deletion_ranges), "our of", length(deletion_ranges), "deletions left after filtration"))
print(paste(length(filtered_snv_ranges), "our of", length(snv_ranges), "samples left after filtration"))
      
start_breakpoints = GRanges(IRanges(start(filtered_deletion_ranges), width=1), seqnames = as.integer(selected_chrom))
start_breakpoint_bins = flank(start_breakpoints, bin_width * num_bins, both=T)
start_tiles = tile(start_breakpoint_bins, n=num_bins)

end_breakpoints = GRanges(IRanges(end(filtered_deletion_ranges), width=1), seqnames = as.integer(selected_chrom))
end_breakpoint_bins = flank(end_breakpoints, bin_width * num_bins, both=T)
end_tiles = tile(end_breakpoint_bins, n=num_bins)

hits = which(filtered_deletion_carrier_mask[,] > 0, arr.ind = T)
pb = progress_bar$new(format=":current/:total [:bar] :percent :elapsed :eta",total = dim(hits)[1])

start_binned_densities = apply(hits, 1, function(x){pb$tick(); countOverlaps(start_tiles[[x[1]]], filtered_snv_ranges[[x[2]]]) / filtered_snv_counts[x[2]];})
end_binned_densities = apply(hits, 1, function(x){pb$tick(); countOverlaps(end_tiles[[x[1]]], filtered_snv_ranges[[x[2]]]) / filtered_snv_counts[x[2]];})

start_var_name = paste("start_breakpoint_density_bins_chrom_", selected_chrom, "_width_", bin_width, "_num_", num_bins, sep="")
assign(start_var_name, start_binned_densities)
save(list=c(start_var_name), file=paste(result_path, "/", start_var_name, ".RData", sep=""))

end_var_name = paste("end_breakpoint_density_bins_chrom_", selected_chrom, "_width_", bin_width, "_num_", num_bins, sep="")
assign(end_var_name, end_binned_densities)
save(list=c(end_var_name), file=paste(result_path, "/", end_var_name, ".RData", sep=""))
