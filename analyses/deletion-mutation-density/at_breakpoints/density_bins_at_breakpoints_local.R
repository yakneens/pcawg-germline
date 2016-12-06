#!/usr/bin/env Rscript

library(VariantAnnotation)
library(data.table)
library(GenomicFeatures)
library(progress)


bin_width = 100
num_bins = 31
donor_meta_path = "~/Downloads/pcawg_data/del_density/input_data/donor_meta.RData"
deletion_ranges_path = "~/Downloads/pcawg_data/del_density/input_data/deletion_ranges.RData"
deletion_info_path = "~/Downloads/pcawg_data/del_density/input_data/deletion_info.RData"
snv_ranges_path = "~/Downloads/pcawg_data/del_density/input_data/snv_ranges.RData"
carrier_mask_path = "~/Downloads/pcawg_data/del_density/input_data/deletion_carrier_mask.RData"

load(donor_meta_path)
load(deletion_ranges_path)
load(deletion_info_path)
load(snv_ranges_path)
load(carrier_mask_path)

chroms = c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22")

breakpoint_margin = 200
bin_based_del_size_cutoff = (bin_width * num_bins) / 2 + breakpoint_margin
print(paste("Deletion size cutoff:", bin_based_del_size_cutoff))



start_bin_list = list()
end_bin_list = list()
for(selected_chrom in chroms){
  print(paste("Processing Chromosome:", selected_chrom))
  #deletion_filter = which(deletion_info$IMPRECISE == T | as.character(seqnames(deletion_ranges)) != selected_chrom)
  del_widths = width(deletion_ranges)
  
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
  
  start_bin_list[[selected_chrom]] = apply(hits, 1, function(x){pb$tick(); countOverlaps(start_tiles[[x[1]]], filtered_snv_ranges[[x[2]]]) / filtered_snv_counts[x[2]];})
  end_bin_list[[selected_chrom]] = apply(hits, 1, function(x){countOverlaps(end_tiles[[x[1]]], filtered_snv_ranges[[x[2]]]) / filtered_snv_counts[x[2]];})
}

breakpoint_densities_100_31_size_10kb_to_100kb = list()
for(i in seq(22)){
  #filepath = paste("~/Downloads/pcawg_data/del_density/at_breakpoints/breakpoint_density_bins_chrom_", i, "_width_1000_num_31.RData", sep="")
  #load(filepath, verbose=T)
  obj_name = paste("breakpoint_density_bins_chrom_", i, "_width_100_num_31", sep="")
  mat_obj = matrix(get(obj_name), nrow=31)
  print(dim(mat_obj))
  breakpoint_densities_100_31_size_10kb_to_100kb[[i]] = mat_obj
}

breakpoint_densities = matrix(unlist(lapply(breakpoint_densities_100_31_del_size_1k_to_10kb, rowSums)), nrow=22, byrow = T)
#breakpoint_densities = matrix(unlist(lapply(breakpoint_densities_100_31_size_10kb_to_100kb, rowSums)), nrow=22, byrow = T)
del_bins = colSums(breakpoint_densities)
del_bins = data.table(seq(1,31), del_bins)

setnames(del_bins, c("index", "bins"))

max_density = max(del_bins$bins) 
midpoint = ceiling(num_bins / 2)
ggplot(del_bins, aes(x=(index - midpoint) * bin_width, y=bins)) + geom_bar(stat="identity", aes(fill=abs(midpoint - index) * bin_width))  + guides(fill=guide_legend(title="Distance from breakpoint")) + xlab("Distance from breakpoint") + ylab("SNV Density") + annotate("text", x=0, y=1.1*max_density, label="Breakpoint")