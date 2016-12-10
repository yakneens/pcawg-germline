library(ggplot2)
library(data.table)
num_bins = 30
bin_width = 100

#Raw data
load("~/Downloads/pcawg_data/del_density/overall/density_bins_carriers_norm.RData", verbose=T)
binned_densities_carriers_norm_dt = as.data.table(do.call(cbind, binned_densities_carriers_norm))
del_list = copy(colnames(binned_densities_carriers_norm_dt))
binned_densities_carriers_norm_dt[, bin_index := seq(num_bins)]
melted_binned_densities_carriers_norm_dt = melt(binned_densities_carriers_norm_dt, id.vars="bin_index",variable.name="del", value.name="snv_density")
del_list = unlist(lapply(del_list, function(x) rep(x, num_bins)))
melted_binned_densities_carriers_norm_dt$del = del_list
save(melted_binned_densities_carriers_norm_dt, file="~/Downloads/temp.RData")
load("~/Downloads/temp.RData")
reduced_melted_binned_densities_carriers_norm_dt = melted_binned_densities_carriers_norm_dt[snv_density > 0][order(bin_index)]
reduced_melted_binned_densities_carriers_norm_dt[,bin_index := as.factor(bin_index)]

del_hit_counts = reduced_melted_binned_densities_carriers_norm_dt[,.N, by=del][order(-N)]
ggplot(del_hit_counts, aes(del, N)) + geom_density(stat="identity")

bin_hits_counts = reduced_melted_binned_densities_carriers_norm_dt[, .N, by=bin_index]
bin_hits_counts[bin_index %in% (10 + seq(10)), mean(N)]
bin_hits_counts[!(bin_index %in% (10 + seq(10))), mean(N)]

ggplot(bin_hits_counts, aes(x=bin_index, y=N)) + 
  geom_bar(stat="identity")  + 
  xlab("Bin") + 
  ylab("# Bins that overlap >=1 SNV") + 
  theme(text=element_text(size=20))

reduced_melted_binned_densities_carriers_norm_dt[order(-snv_density)][1:20]

look_up_deletion_info <- function(my_del){
  print(deletion_ranges[my_del])
  print(width(deletion_ranges[my_del]))
  
  print(donor_meta[which(deletion_carrier_mask[my_del,] > 0)]$dcc_project_code)
  print(reduced_melted_binned_densities_carriers_norm_dt[del == my_del])
  
  donor_hits = as.vector(colnames(deletion_carrier_mask[which(deletion_carrier_mask[my_del,] > 0)])[ceiling(which(melted_binned_densities_carriers_norm_dt[del == my_del]$snv_density > 0) / 30)])
  print(donor_hits)
  print(donor_meta[donor_unique_id %in% donor_hits])
  
  pre_deletion_ranges = trim(flank(deletion_ranges[my_del], width(deletion_ranges[my_del]), start = T))
  post_deletion_ranges = trim(flank(deletion_ranges[my_del], width(deletion_ranges[my_del]), start = F))
  
  deletion_tiles = tile(deletion_ranges[my_del], n = as.integer(10))
  pre_deletion_tiles = tile(pre_deletion_ranges, n = as.integer(10))
  post_deletion_tiles = tile(post_deletion_ranges, n = as.integer(10))
  
  all_tiles = pc(pre_deletion_tiles, deletion_tiles, post_deletion_tiles)
  
  for(my_donor in donor_hits){
    overlaps = countOverlaps(all_tiles[[1]], snv_ranges[my_donor][[1]])
    print(overlaps)
    print(snv_ranges[my_donor][[1]][which(as.character(seqnames(snv_ranges[my_donor][[1]])) == seqnames(deletion_ranges[my_del]))])
  }
}

#Bin 10
look_up_deletion_info("DEL00047078")
look_up_deletion_info("DEL00022110")

#Bin 20
look_up_deletion_info("DEL00221829")
look_up_deletion_info("DEL00132251")
look_up_deletion_info("DEL00222020")

#Bin 26
look_up_deletion_info("DEL00068818")



ggplot(reduced_melted_binned_densities_carriers_norm_dt[, mean(snv_density), by=bin_index], aes(x=bin_index, y=V1)) + 
  geom_bar(stat="identity") + 
  xlab("Bin") + 
  ylab("Mean SNV density")
ggplot(reduced_melted_binned_densities_carriers_norm_dt[, median(snv_density), by=bin_index], aes(x=bin_index, y=V1)) + 
  geom_bar(stat="identity") + 
  xlab("Bin") + 
  ylab("Median SNV density")
ggplot(reduced_melted_binned_densities_carriers_norm_dt, aes(bin_index, log10(snv_density))) + 
  geom_boxplot() + 
  xlab("Bin") + 
  ylab("Log10 SNV density")

reduced_melted_binned_densities_carriers_norm_dt[as.integer(bin_index) >= 1 & as.integer(bin_index) <= 10, bin_type := "Left Flank"]
reduced_melted_binned_densities_carriers_norm_dt[as.integer(bin_index) >= 11 & as.integer(bin_index) <= 20, bin_type := "Deletion"]
reduced_melted_binned_densities_carriers_norm_dt[as.integer(bin_index) >= 21 & as.integer(bin_index) <= 30, bin_type := "Right Flank"]
reduced_melted_binned_densities_carriers_norm_dt[, bin_type := factor(bin_type, levels=c("Left Flank", "Deletion", "Right Flank"), ordered=T)]
ggplot(reduced_melted_binned_densities_carriers_norm_dt) + 
  geom_density(aes(log10(snv_density), colour=bin_index)) + 
  facet_grid(.~bin_type) +
  xlab("Log10 SNV density") + 
  ylab("Density")





load("~/Downloads/pcawg_data/del_density/overall/summed_density_bins_carriers_norm.RData", verbose=T)
summed_binned_densities_carriers_for_plot = as.data.table(summed_binned_densities_carriers_norm)
setnames(summed_binned_densities_carriers_for_plot, as.character(seq(num_bins)))

#Normalize to 1
#summed_binned_densities_carriers_for_plot = summed_binned_densities_carriers_for_plot / summed_binned_densities_carriers_for_plot[,do.call(pmax, .SD)]

summed_binned_densities_carriers_for_plot[, chrom := seq(22)]

melted_summed_binned_densities_carriers_for_plot = melt(summed_binned_densities_carriers_for_plot, id.vars = c("chrom"))
setnames(melted_summed_binned_densities_carriers_for_plot, c("chrom", "bin_index", "snv_density"))
melted_summed_binned_densities_carriers_for_plot[, chrom := as.factor(chrom)][,bin_index := as.integer(bin_index)]

max_density = max(melted_summed_binned_densities_carriers_for_plot$snv_density) 
ggplot(melted_summed_binned_densities_carriers_for_plot, aes(x=bin_index, y=snv_density)) + 
  geom_bar(stat="identity") + 
  geom_vline(xintercept=10.5, linetype="dashed") + 
  geom_vline(xintercept=20.5, linetype="dashed") + 
  facet_wrap(~chrom, ncol=3) + 
  xlab("Bin") + 
  ylab("SNV Density") + 
  annotate("text", x=5, y=1.1*max_density, label="Left Flank") + 
  annotate("text", x=15, y=1.1*max_density, label="Deletion") + 
  annotate("text", x=25, y=1.1*max_density, label="Right Flank")

#Raw counts

load("~/Downloads/pcawg_data/del_density/overall/density_bins_raw_counts_carriers_norm.RData", verbose=T)

binned_densities_raw_counts_carriers_norm_dt = as.data.table(do.call(cbind, binned_densities_raw_counts_carriers_norm))
del_list = copy(colnames(binned_densities_raw_counts_carriers_norm_dt))
binned_densities_raw_counts_carriers_norm_dt[, bin_index := seq(num_bins)]
melted_binned_densities_raw_counts_carriers_norm_dt = melt(binned_densities_raw_counts_carriers_norm_dt, id.vars="bin_index",variable.name="del", value.name="snv_count")
del_list = unlist(lapply(del_list, function(x) rep(x, num_bins)))
melted_binned_densities_raw_counts_carriers_norm_dt$del = del_list
reduced_melted_binned_densities_raw_counts_carriers_norm_dt = melted_binned_densities_raw_counts_carriers_norm_dt[snv_count > 0][order(bin_index)]
reduced_melted_binned_densities_raw_counts_carriers_norm_dt[,bin_index := as.factor(bin_index)]

del_hit_counts = reduced_melted_binned_densities_raw_counts_carriers_norm_dt[,.N, by=del][order(-N)]
ggplot(del_hit_counts, aes(reorder(del, -N), N)) + geom_density(stat="identity")

bin_hits_counts = reduced_melted_binned_densities_raw_counts_carriers_norm_dt[, .N, by=bin_index]
bin_hits_counts[bin_index %in% (10 + seq(10)), mean(N)]
bin_hits_counts[!(bin_index %in% (10 + seq(10))), mean(N)]

ggplot(bin_hits_counts, aes(x=bin_index, y=N)) + 
  geom_bar(stat="identity")  + 
  xlab("Bin") + 
  ylab("# Bins that overlap >=1 SNV") + 
  theme(text=element_text(size=20))

ggplot(reduced_melted_binned_densities_raw_counts_carriers_norm_dt[, mean(snv_count), by=bin_index], aes(x=bin_index, y=V1)) + 
  geom_bar(stat="identity") + 
  xlab("Bin") + 
  ylab("Mean SNV density")
ggplot(reduced_melted_binned_densities_raw_counts_carriers_norm_dt[, median(snv_count), by=bin_index], aes(x=bin_index, y=V1)) + 
  geom_bar(stat="identity") + 
  xlab("Bin") + 
  ylab("Median SNV density")
ggplot(reduced_melted_binned_densities_raw_counts_carriers_norm_dt, aes(factor(bin_index), snv_count)) + 
  geom_boxplot() + 
  xlab("Bin") + 
  ylab("Log10 SNV density")

load("~/Downloads/pcawg_data/del_density/overall/summed_density_bins_raw_counts_carriers.RData", verbose=T)

#Donors with more than 1K SNV
load("~/Downloads/pcawg_data/del_density/overall/density_bins_carriers_1k_snv_norm.RData", verbose=T)
binned_densities_carriers_1k_snv_norm_dt = as.data.table(do.call(cbind, binned_densities_carriers_1k_snv_norm))
del_list = copy(colnames(binned_densities_carriers_1k_snv_norm_dt))
binned_densities_carriers_1k_snv_norm_dt[, bin_index := seq(num_bins)]
melted_binned_densities_carriers_1k_snv_norm_dt = melt(binned_densities_carriers_1k_snv_norm_dt, id.vars="bin_index",variable.name="del", value.name="snv_density")
del_list = unlist(lapply(del_list, function(x) rep(x, num_bins)))
melted_binned_densities_carriers_1k_snv_norm_dt$del = del_list
save(melted_binned_densities_carriers_1k_snv_norm_dt, file="~/Downloads/temp.RData")
load("~/Downloads/temp.RData")
reduced_melted_binned_densities_carriers_1k_snv_norm_dt = melted_binned_densities_carriers_1k_snv_norm_dt[snv_density > 0][order(bin_index)]
reduced_melted_binned_densities_carriers_1k_snv_norm_dt[,bin_index := as.factor(bin_index)]

del_hit_counts = reduced_melted_binned_densities_carriers_1k_snv_norm_dt[,.N, by=del][order(-N)]
ggplot(del_hit_counts, aes(del, N)) + geom_density(stat="identity")

bin_hits_counts = reduced_melted_binned_densities_carriers_1k_snv_norm_dt[, .N, by=bin_index]
bin_hits_counts[bin_index %in% (10 + seq(10)), mean(N)]
bin_hits_counts[!(bin_index %in% (10 + seq(10))), mean(N)]

ggplot(bin_hits_counts, aes(x=bin_index, y=N)) + 
  geom_bar(stat="identity")  + 
  xlab("Bin") + 
  ylab("# Bins that overlap >=1 SNV") + 
  theme(text=element_text(size=20))

reduced_melted_binned_densities_carriers_1k_snv_norm_dt[order(-snv_density)][1:20]

look_up_deletion_info <- function(my_del){
  print(deletion_ranges[my_del])
  print(width(deletion_ranges[my_del]))
  
  print(donor_meta[which(deletion_carrier_mask[my_del,] > 0)]$dcc_project_code)
  print(reduced_melted_binned_densities_carriers_1k_snv_norm_dt[del == my_del])
  
  donor_hits = as.vector(colnames(deletion_carrier_mask[which(deletion_carrier_mask[my_del,] > 0)])[ceiling(which(melted_binned_densities_carriers_1k_snv_norm_dt[del == my_del]$snv_density > 0) / 30)])
  print(donor_hits)
  print(donor_meta[donor_unique_id %in% donor_hits])
  
  pre_deletion_ranges = trim(flank(deletion_ranges[my_del], width(deletion_ranges[my_del]), start = T))
  post_deletion_ranges = trim(flank(deletion_ranges[my_del], width(deletion_ranges[my_del]), start = F))
  
  deletion_tiles = tile(deletion_ranges[my_del], n = as.integer(10))
  pre_deletion_tiles = tile(pre_deletion_ranges, n = as.integer(10))
  post_deletion_tiles = tile(post_deletion_ranges, n = as.integer(10))
  
  all_tiles = pc(pre_deletion_tiles, deletion_tiles, post_deletion_tiles)
  
  for(my_donor in donor_hits){
    overlaps = countOverlaps(all_tiles[[1]], snv_ranges[my_donor][[1]])
    print(overlaps)
    print(snv_ranges[my_donor][[1]][which(as.character(seqnames(snv_ranges[my_donor][[1]])) == seqnames(deletion_ranges[my_del]))])
  }
}

#Bin 10
look_up_deletion_info("DEL00047078")
look_up_deletion_info("DEL00022110")

#Bin 20
look_up_deletion_info("DEL00221829")
look_up_deletion_info("DEL00132251")
look_up_deletion_info("DEL00222020")

#Bin 26
look_up_deletion_info("DEL00068818")



ggplot(reduced_melted_binned_densities_carriers_1k_snv_norm_dt[, mean(snv_density), by=bin_index], aes(x=bin_index, y=V1)) + 
  geom_bar(stat="identity") + 
  xlab("Bin") + 
  ylab("Mean SNV density")
ggplot(reduced_melted_binned_densities_carriers_1k_snv_norm_dt[, median(snv_density), by=bin_index], aes(x=bin_index, y=V1)) + 
  geom_bar(stat="identity") + 
  xlab("Bin") + 
  ylab("Median SNV density")
ggplot(reduced_melted_binned_densities_carriers_1k_snv_norm_dt, aes(bin_index, log10(snv_density))) + 
  geom_boxplot() + 
  xlab("Bin") + 
  ylab("Log10 SNV density")

reduced_melted_binned_densities_carriers_1k_snv_norm_dt[as.integer(bin_index) >= 1 & as.integer(bin_index) <= 10, bin_type := "Left Flank"]
reduced_melted_binned_densities_carriers_1k_snv_norm_dt[as.integer(bin_index) >= 11 & as.integer(bin_index) <= 20, bin_type := "Deletion"]
reduced_melted_binned_densities_carriers_1k_snv_norm_dt[as.integer(bin_index) >= 21 & as.integer(bin_index) <= 30, bin_type := "Right Flank"]
reduced_melted_binned_densities_carriers_1k_snv_norm_dt[, bin_type := factor(bin_type, levels=c("Left Flank", "Deletion", "Right Flank"), ordered=T)]
ggplot(reduced_melted_binned_densities_carriers_1k_snv_norm_dt) + 
  geom_density(aes(log10(snv_density), colour=bin_index)) + 
  facet_grid(.~bin_type) +
  xlab("Log10 SNV density") + 
  ylab("Density")





load("~/Downloads/pcawg_data/del_density/overall/summed_density_bins_carriers_1k_snv_norm.RData", verbose=T)
summed_binned_densities_carriers_1k_snv_for_plot = as.data.table(summed_binned_densities_carriers_1k_snv_norm)
setnames(summed_binned_densities_carriers_1k_snv_for_plot, as.character(seq(num_bins)))

#Normalize to 1
#summed_binned_densities_carriers_1k_snv_for_plot = summed_binned_densities_carriers_1k_snv_for_plot / summed_binned_densities_carriers_1k_snv_for_plot[,do.call(pmax, .SD)]

summed_binned_densities_carriers_1k_snv_for_plot[, chrom := seq(22)]

melted_summed_binned_densities_carriers_1k_snv_for_plot = melt(summed_binned_densities_carriers_1k_snv_for_plot, id.vars = c("chrom"))
setnames(melted_summed_binned_densities_carriers_1k_snv_for_plot, c("chrom", "bin_index", "snv_density"))
melted_summed_binned_densities_carriers_1k_snv_for_plot[, chrom := as.factor(chrom)][,bin_index := as.integer(bin_index)]

max_density = max(melted_summed_binned_densities_carriers_1k_snv_for_plot$snv_density) 
ggplot(melted_summed_binned_densities_carriers_1k_snv_for_plot, aes(x=bin_index, y=snv_density)) + 
  geom_bar(stat="identity") + 
  geom_vline(xintercept=10.5, linetype="dashed") + 
  geom_vline(xintercept=20.5, linetype="dashed") + 
  facet_wrap(~chrom, ncol=3) + 
  xlab("Bin") + 
  ylab("SNV Density") + 
  annotate("text", x=5, y=1.1*max_density, label="Left Flank") + 
  annotate("text", x=15, y=1.1*max_density, label="Deletion") + 
  annotate("text", x=25, y=1.1*max_density, label="Right Flank")

