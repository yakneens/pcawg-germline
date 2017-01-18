library(ggplot2)
library(scales)
library(data.table)

num_bins = 30
midpoint = ceiling(num_bins / 2)
bin_width = 100

bar_width = 90

load("~/Downloads/pcawg_data/del_density/at_breakpoints/breakpoint_densities_all_100_30_200_carriers_start.RData", verbose = T)


breakpoint_densities_100_30_200_carriers_start_dt = as.data.table(do.call(cbind, breakpoint_densities_100_30_200_carriers_start))
del_list = copy(colnames(breakpoint_densities_100_30_200_carriers_start_dt))
breakpoint_densities_100_30_200_carriers_start_dt[, bin_index := seq(num_bins)]
melted_breakpoint_densities_100_30_200_carriers_start_dt = melt(breakpoint_densities_100_30_200_carriers_start_dt, id.vars="bin_index",variable.name="del", value.name="snv_density")
del_list = unlist(lapply(del_list, function(x) rep(x, num_bins)))
melted_breakpoint_densities_100_30_200_carriers_start_dt$del = del_list
melted_breakpoint_densities_100_30_200_carriers_start_dt[, breakpoint_type := "start"]
#save(melted_breakpoint_densities_100_30_200_carriers_start_dt, file="temp_start.RData")
load("temp_start.RData", verbose=T)

load("~/Downloads/pcawg_data/del_density/at_breakpoints/breakpoint_densities_all_100_30_200_carriers_end.RData", verbose = T)

breakpoint_densities_100_30_200_carriers_end_dt = as.data.table(do.call(cbind, breakpoint_densities_100_30_200_carriers_end))
del_list = copy(colnames(breakpoint_densities_100_30_200_carriers_end_dt))
breakpoint_densities_100_30_200_carriers_end_dt[, bin_index := seq(num_bins)]
melted_breakpoint_densities_100_30_200_carriers_end_dt = melt(breakpoint_densities_100_30_200_carriers_end_dt, id.vars="bin_index",variable.name="del", value.name="snv_density")
del_list = unlist(lapply(del_list, function(x) rep(x, num_bins)))
melted_breakpoint_densities_100_30_200_carriers_end_dt$del = del_list
melted_breakpoint_densities_100_30_200_carriers_end_dt[, breakpoint_type := "end"]
#save(melted_breakpoint_densities_100_30_200_carriers_end_dt, file="temp_end.RData")
load("temp_end.RData", verbose=T)

melted_breakpoint_densities = rbind(melted_breakpoint_densities_100_30_200_carriers_start_dt, melted_breakpoint_densities_100_30_200_carriers_end_dt)
melted_breakpoint_densities$breakpoint_type = factor(melted_breakpoint_densities$breakpoint_type, levels=c("start", "end"))

reduced_melted_breakpoint_densities = melted_breakpoint_densities[snv_density > 0][order(bin_index)]
reduced_melted_breakpoint_densities[,bin_index := as.factor(bin_index)]

del_hit_counts = reduced_melted_breakpoint_densities[,.N, by=c("del","breakpoint_type")][order(-N)]
ggplot(del_hit_counts, aes(reorder(del, -N), N, fill=breakpoint_type)) + 
  geom_bar(stat = "identity") + 
  facet_grid(breakpoint_type~.) +
  xlab("Deletion") +
  ylab("# of bins that overlap SNVs")

bin_hits_counts = reduced_melted_breakpoint_densities[, .N, by=c("bin_index","breakpoint_type")]
bin_hits_counts[bin_index %in% (seq(15)), mean(N), by=breakpoint_type]
bin_hits_counts[!(bin_index %in% (seq(15))), mean(N), by=breakpoint_type]

max_bin_hit_count = max(bin_hits_counts$N)

ggplot(bin_hits_counts, aes(x=(as.integer(bin_index) - midpoint) * bin_width - (bar_width/2), y=N, fill=breakpoint_type, width=bar_width)) + 
  geom_bar(stat="identity")  +
  geom_segment(x = 0, xend=0, y=0, yend=1.05*max_bin_hit_count, linetype="dotted") +
  annotate("text", x=0, y= 1.08*max_bin_hit_count, label="Breakpoint") +9
  scale_x_continuous("Distance from breakpoint", sec.axis = sec_axis(~ (. + (bar_width/2))/bin_width + midpoint, name = "Bin")) +
  facet_grid(.~breakpoint_type) +
  ylab("# Deletions that overlap >=1 SNV") + 
  theme(text=element_text(size=10), legend.position = "bottom") 

bin_16 = reduced_melted_breakpoint_densities[order(-snv_density)][bin_index == 16 & breakpoint_type =="start"]
mean(bin_16$snv_density)
var(bin_16$snv_density)
bin_15 = reduced_melted_breakpoint_densities[order(-snv_density)][bin_index == 15 & breakpoint_type =="start"]

data_for_plot = reduced_melted_breakpoint_densities[, .(mean_density = mean(snv_density)), by=c("bin_index","breakpoint_type")]
max_y_val = max(data_for_plot$mean_density)

data_for_plot[, .(mean(mean_density), var(mean_density)), by=breakpoint_type]

ggplot(data_for_plot, aes(x=(as.integer(bin_index) - midpoint) * bin_width - (bar_width/2), y=mean_density, fill=breakpoint_type, width=bar_width)) + 
  geom_bar(stat="identity") +
  geom_segment(x = 0, xend=0, y=0, yend=1.05*max_y_val, linetype="dotted") +
  annotate("text", x=0, y=1.08*max_y_val, label="Breakpoint") +
  scale_x_continuous("Distance from breakpoint", sec.axis = sec_axis(~ (. + (bar_width/2))/bin_width + midpoint, name = "Bin")) +
  facet_grid(.~breakpoint_type) +
  ylab("Mean SNV density")

#Distribution of mean density
ggplot(data_for_plot, aes(mean_density, colour=breakpoint_type)) + geom_density() + facet_grid(.~breakpoint_type)

ggplot(data_for_plot[breakpoint_type == "start"]) + 
  geom_density(aes(mean_density)) +
  geom_vline(xintercept = 7.293306e-05)

data_for_plot = reduced_melted_breakpoint_densities[, .(median_density = median(snv_density)), by=c("bin_index","breakpoint_type")]
max_y_val = max(data_for_plot$median_density)
ggplot(data_for_plot, aes(x=(as.integer(bin_index) - midpoint) * bin_width - (bar_width/2), y=median_density, fill=breakpoint_type, width=bar_width)) + 
  geom_bar(stat="identity") + 
  geom_segment(x = 0, xend=0, y=0, yend=1.05*max_y_val, linetype="dotted") +
  annotate("text", x=0, y=1.08*max_y_val, label="Breakpoint") +
  scale_x_continuous("Distance from breakpoint", sec.axis = sec_axis(~ (. + (bar_width/2))/bin_width + midpoint, name = "Bin")) +
  facet_grid(.~breakpoint_type) +
  ylab("Median SNV density")


ggplot(reduced_melted_breakpoint_densities, aes(bin_index, log(snv_density))) + 
  geom_boxplot() + 
  facet_grid(.~breakpoint_type) +
  xlab("Bin") + 
  ylab("Log10 SNV density")

reduced_melted_breakpoint_densities[as.integer(bin_index) >= 1 & as.integer(bin_index) <= 10, bin_type := "Left Flank"]
reduced_melted_breakpoint_densities[as.integer(bin_index) >= 11 & as.integer(bin_index) <= 20, bin_type := "Breakpoint"]
reduced_melted_breakpoint_densities[as.integer(bin_index) >= 21 & as.integer(bin_index) <= 30, bin_type := "Right Flank"]
reduced_melted_breakpoint_densities[, bin_type := factor(bin_type, levels=c("Left Flank", "Breakpoint", "Right Flank"), ordered=T)]

ggplot(reduced_melted_breakpoint_densities[as.integer(bin_index) %in% c(1,13,14,15,16,17,18, 30)]) + 
  geom_density(aes(log10(snv_density), colour=bin_index)) + 
  facet_wrap(breakpoint_type~bin_index, ncol=8) +
  xlab("Log10 SNV density") + 
  ylab("Density")

data_for_plot = reduced_melted_breakpoint_densities[, .(snv_density = sum(snv_density)), by=c("del","bin_index", "breakpoint_type")][, total_snv_density := sum(snv_density), by=c("bin_index", "breakpoint_type")][, scaled_snv_density := snv_density / total_snv_density][order(bin_index, breakpoint_type, -scaled_snv_density)]
ggplot(data_for_plot[breakpoint_type=="end",.SD[1:5],by=c("bin_index", "breakpoint_type")], aes(x=bin_index, y=scaled_snv_density, fill=reorder(del, -scaled_snv_density), label=scales::percent(round(scaled_snv_density, digits = 2)))) + 
  geom_bar(stat="identity", position="stack") + 
  geom_text(size = 3, position = position_stack(vjust = 0.5)) + 
  guides(fill=FALSE) +
  facet_grid(.~breakpoint_type) +
  scale_y_continuous(labels = scales::percent) +
  xlab("Bin") + 
  ylab("% Total Density")

#Raw counts
load("~/Downloads/pcawg_data/del_density/at_breakpoints/breakpoint_densities_all_100_30_200_raw_counts_carriers_start.RData", verbose = T)


breakpoint_densities_100_30_200_raw_counts_carriers_start_dt = as.data.table(do.call(cbind, breakpoint_densities_100_30_200_raw_counts_carriers_start))
del_list = copy(colnames(breakpoint_densities_100_30_200_raw_counts_carriers_start_dt))
breakpoint_densities_100_30_200_raw_counts_carriers_start_dt[, bin_index := seq(num_bins)]
melted_breakpoint_densities_100_30_200_raw_counts_carriers_start_dt = melt(breakpoint_densities_100_30_200_raw_counts_carriers_start_dt, id.vars="bin_index",variable.name="del", value.name="snv_count")
del_list = unlist(lapply(del_list, function(x) rep(x, num_bins)))
melted_breakpoint_densities_100_30_200_raw_counts_carriers_start_dt$del = del_list
melted_breakpoint_densities_100_30_200_raw_counts_carriers_start_dt[, breakpoint_type := "start"]
#save(melted_breakpoint_densities_100_30_200_carriers_start_dt, file="temp_start.RData")
load("temp_start.RData", verbose=T)

load("~/Downloads/pcawg_data/del_density/at_breakpoints/breakpoint_densities_all_100_30_200_raw_counts_carriers_end.RData", verbose = T)

breakpoint_densities_100_30_200_raw_counts_carriers_end_dt = as.data.table(do.call(cbind, breakpoint_densities_100_30_200_raw_counts_carriers_end))
del_list = copy(colnames(breakpoint_densities_100_30_200_raw_counts_carriers_end_dt))
breakpoint_densities_100_30_200_raw_counts_carriers_end_dt[, bin_index := seq(num_bins)]
melted_breakpoint_densities_100_30_200_raw_counts_carriers_end_dt = melt(breakpoint_densities_100_30_200_raw_counts_carriers_end_dt, id.vars="bin_index",variable.name="del", value.name="snv_count")
del_list = unlist(lapply(del_list, function(x) rep(x, num_bins)))
melted_breakpoint_densities_100_30_200_raw_counts_carriers_end_dt$del = del_list
melted_breakpoint_densities_100_30_200_raw_counts_carriers_end_dt[, breakpoint_type := "end"]
#save(melted_breakpoint_densities_100_30_200_carriers_end_dt, file="temp_end.RData")
load("temp_end.RData", verbose=T)

melted_breakpoint_densities_raw_counts = rbind(melted_breakpoint_densities_100_30_200_raw_counts_carriers_start_dt, melted_breakpoint_densities_100_30_200_raw_counts_carriers_end_dt)
melted_breakpoint_densities_raw_counts$breakpoint_type = factor(melted_breakpoint_densities_raw_counts$breakpoint_type, levels=c("start", "end"))

reduced_melted_breakpoint_densities_raw_counts = melted_breakpoint_densities_raw_counts[snv_count > 0][order(bin_index)]
reduced_melted_breakpoint_densities_raw_counts[,bin_index := as.factor(bin_index)]

ggplot(reduced_melted_breakpoint_densities_raw_counts, aes(x=snv_count)) + geom_histogram() + scale_y_log10()
