#!/usr/bin/env Rscript

library(VariantAnnotation)
library(data.table)
library(GenomicFeatures)
library(progress)



donor_meta_path = "~/Downloads/pcawg_data/del_density/by_project/donor_meta.RData"
deletion_ranges_path = "~/Downloads/pcawg_data/del_density/by_project/deletion_ranges.RData"
deletion_info_path = "~/Downloads/pcawg_data/del_density/by_project/deletion_info.RData"
snv_ranges_path = "~/Downloads/pcawg_data/del_density/by_project/snv_ranges.RData"
carrier_mask_path = "~/Downloads/pcawg_data/del_density/by_project/deletion_carrier_mask.RData"

load(donor_meta_path)
load(deletion_ranges_path)
load(deletion_info_path)
load(snv_ranges_path)
load(carrier_mask_path)


c_to_a_snv_ranges = lapply(snv_ranges, function(x) x[which((as.character(x$REF) == "C" & as.character(unlist(x$ALT)) == "A") | (as.character(x$REF) == "G" & as.character(unlist(x$ALT)) == "T"))])
c_to_a_snv_counts = unlist(lapply(c_to_a_snv_ranges, length))

c_to_g_snv_ranges = lapply(snv_ranges, function(x) x[which((as.character(x$REF) == "C" & as.character(unlist(x$ALT)) == "G") | (as.character(x$REF) == "G" & as.character(unlist(x$ALT)) == "C"))])
c_to_g_snv_counts = unlist(lapply(c_to_g_snv_ranges, length))

c_to_t_snv_ranges = lapply(snv_ranges, function(x) x[which((as.character(x$REF) == "C" & as.character(unlist(x$ALT)) == "T") | (as.character(x$REF) == "G" & as.character(unlist(x$ALT)) == "A"))])
c_to_t_snv_counts = unlist(lapply(c_to_t_snv_ranges, length))

t_to_g_snv_ranges = lapply(snv_ranges, function(x) x[which((as.character(x$REF) == "T" & as.character(unlist(x$ALT)) == "G") | (as.character(x$REF) == "A" & as.character(unlist(x$ALT)) == "C"))])
t_to_g_snv_counts = unlist(lapply(t_to_g_snv_ranges, length))

t_to_a_snv_ranges = lapply(snv_ranges, function(x) x[which((as.character(x$REF) == "T" & as.character(unlist(x$ALT)) == "A") | (as.character(x$REF) == "A" & as.character(unlist(x$ALT)) == "T"))])
t_to_a_snv_counts = unlist(lapply(t_to_a_snv_ranges, length))

t_to_c_snv_ranges = lapply(snv_ranges, function(x) x[which((as.character(x$REF) == "T" & as.character(unlist(x$ALT)) == "C") | (as.character(x$REF) == "A" & as.character(unlist(x$ALT)) == "G"))])
t_to_c_snv_counts = unlist(lapply(t_to_c_snv_ranges, length))

chroms = c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22")


c_to_a_binned_densities = list()
c_to_t_binned_densities = list()
c_to_g_binned_densities = list()
t_to_g_binned_densities = list()
t_to_c_binned_densities = list()
t_to_a_binned_densities = list()
snv_counts = unlist(lapply(snv_ranges, length))

for(selected_chrom in chroms){
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
  hits = hits[-which(hits[,2] == 2519),] 
  pb = progress_bar$new(format=":current/:total [:bar] :percent :elapsed :eta",total = dim(hits)[1])
  
  c_to_a_binned_densities[[selected_chrom]] = apply(hits, 1, function(x){pb$tick(); countOverlaps(all_tiles[[x[1]]], c_to_a_snv_ranges[[x[2]]]) / snv_counts[x[2]];})
  c_to_t_binned_densities[[selected_chrom]] = apply(hits, 1, function(x){countOverlaps(all_tiles[[x[1]]], c_to_t_snv_ranges[[x[2]]]) / snv_counts[x[2]];})
  c_to_g_binned_densities[[selected_chrom]] = apply(hits, 1, function(x){countOverlaps(all_tiles[[x[1]]], c_to_g_snv_ranges[[x[2]]]) / snv_counts[x[2]];})
  t_to_g_binned_densities[[selected_chrom]] = apply(hits, 1, function(x){countOverlaps(all_tiles[[x[1]]], t_to_g_snv_ranges[[x[2]]]) / snv_counts[x[2]];})
  t_to_c_binned_densities[[selected_chrom]] = apply(hits, 1, function(x){countOverlaps(all_tiles[[x[1]]], t_to_c_snv_ranges[[x[2]]]) / snv_counts[x[2]];})
  t_to_a_binned_densities[[selected_chrom]] = apply(hits, 1, function(x){countOverlaps(all_tiles[[x[1]]], t_to_a_snv_ranges[[x[2]]]) / snv_counts[x[2]];})
}


c_to_a_aggregated_densities = data.table(colSums(matrix(unlist(lapply(c_to_a_binned_densities, rowSums)), nrow=22, byrow = T)), seq(30), "C_to_A")
c_to_g_aggregated_densities = data.table(colSums(matrix(unlist(lapply(c_to_g_binned_densities, rowSums)), nrow=22, byrow = T)), seq(30), "C_to_G")
c_to_t_aggregated_densities = data.table(colSums(matrix(unlist(lapply(c_to_t_binned_densities, rowSums)), nrow=22, byrow = T)), seq(30), "C_to_T")
t_to_g_aggregated_densities = data.table(colSums(matrix(unlist(lapply(t_to_g_binned_densities, rowSums)), nrow=22, byrow = T)), seq(30), "T_to_G")
t_to_c_aggregated_densities = data.table(colSums(matrix(unlist(lapply(t_to_c_binned_densities, rowSums)), nrow=22, byrow = T)), seq(30), "T_to_C")
t_to_a_aggregated_densities = data.table(colSums(matrix(unlist(lapply(t_to_a_binned_densities, rowSums)), nrow=22, byrow = T)), seq(30), "T_to_A")

aggregated_densities = rbind(c_to_a_aggregated_densities,
                             c_to_g_aggregated_densities,
                             c_to_t_aggregated_densities,
                             t_to_g_aggregated_densities,
                             t_to_c_aggregated_densities,
                             t_to_a_aggregated_densities)
setnames(aggregated_densities, c("bins", "index", "sub_type"))

max_density = max(aggregated_densities$bins)
ggplot(aggregated_densities, aes(x=index, y=bins)) + geom_bar(stat="identity", fill="grey70", aes(colour=sub_type)) + geom_vline(xintercept=10.5, linetype="dashed") + geom_vline(xintercept=20.5, linetype="dashed") + xlab("Bin") + ylab("SNV Density") + annotate("text", x=5, y=1.1*max_density, label="Left Flank") + annotate("text", x=15, y=1.1*max_density, label="Deletion") + annotate("text", x=25, y=1.1*max_density, label="Right Flank")
