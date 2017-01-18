library(ggplot2)
library(data.table)

#By substitution type
load("~/Downloads/pcawg_data/del_density/by_substitution/tabular_summed_binned_densities_by_sub_type_carriers_norm.RData", verbose=T)
binned_densities_by_sub_type_carriers_norm_for_plot = tabular_summed_binned_densities_by_sub_type_carriers_norm[,.(snv_density=rowSums(.SD)),.(sub_type=x, bin_index = V2),.SDcols=as.character(seq(22))]

binned_densities_by_sub_type_carriers_norm_for_plot = melt(as.data.table(summed_merged_binned_densities_by_sub_type)[,bin_index:=seq(30)], id.vars = "bin_index", variable.name = "sub_type", value.name="snv_density")


summed_binned_densities_by_sub_type_carriers_dt = as.data.table(do.call(rbind, lapply(binned_densities, rowSums)), c("c_to_a", "c_to_t", "c_to_g", "t_to_g", "t_to_c", "t_to_a"))
setnames(summed_binned_densities_by_sub_type_carriers_dt, c("sub_type", as.character(seq(30))))
melted_summed_binned_densities_by_sub_type_carriers_dt = melt(summed_binned_densities_by_sub_type_carriers_dt, id.vars = c("sub_type"), variable.name = "bin_index", value.name = "snv_density")

reduced_melted_binned_densities_by_sub_type_carriers_dt = melted_binned_densities_by_sub_type_carriers_dt[snv_density > 0]


max_density = max(melted_summed_binned_densities_by_sub_type_carriers_dt$snv_density) 
ggplot(melted_summed_binned_densities_by_sub_type_carriers_dt, aes(x=bin_index, y=snv_density,fill=sub_type)) + 
  geom_bar(stat="identity") + 
  geom_vline(xintercept=10.5, linetype="dashed") + 
  geom_vline(xintercept=20.5, linetype="dashed") + 
  facet_wrap(~sub_type, nrow=2, ncol=3) + 
  guides(fill=guide_legend(title="Substitution Type")) + 
  xlab("Bin") + 
  ylab("SNV Density") + 
  annotate("text", x=5, y=1.1*max_density, label="Left Flank") + 
  annotate("text", x=15, y=1.1*max_density, label="Deletion") + 
  annotate("text", x=25, y=1.1*max_density, label="Right Flank")

ggplot(melted_summed_binned_densities_by_sub_type_carriers_dt, aes(x=bin_index, y=snv_density,fill=sub_type)) + 
  geom_bar(stat="identity", position="stack") + 
  geom_vline(xintercept=10.5, linetype="dashed") + 
  geom_vline(xintercept=20.5, linetype="dashed") + 
  guides(fill=guide_legend(title="Substitution Type")) + 
  xlab("Bin") + 
  ylab("SNV Density") + 
  annotate("text", x=5, y=1.1*max_density, label="Left Flank") + 
  annotate("text", x=15, y=1.1*max_density, label="Deletion") + 
  annotate("text", x=25, y=1.1*max_density, label="Right Flank")

#Corrected densities
corrected_binned_densities_by_sub_type_carriers_norm_for_plot = copy(binned_densities_by_sub_type_carriers_norm_for_plot)
corrected_binned_densities_by_sub_type_carriers_norm_for_plot = corrected_binned_densities_by_sub_type_carriers_norm_for_plot[bin_index > 10 & bin_index <= 20, snv_density := 2 * snv_density]

max_density = max(corrected_binned_densities_by_sub_type_carriers_norm_for_plot$snv_density) 
ggplot(corrected_binned_densities_by_sub_type_carriers_norm_for_plot, aes(x=bin_index, y=snv_density,fill=sub_type)) + 
  geom_bar(stat="identity") + 
  geom_vline(xintercept=10.5, linetype="dashed") + 
  geom_vline(xintercept=20.5, linetype="dashed") + 
  facet_wrap(~sub_type, nrow=2, ncol=3) + 
  guides(fill=guide_legend(title="Substitution Type")) + 
  xlab("Bin") + 
  ylab("SNV Density") + 
  annotate("text", x=5, y=1.1*max_density, label="Left Flank") + 
  annotate("text", x=15, y=1.1*max_density, label="Deletion") + 
  annotate("text", x=25, y=1.1*max_density, label="Right Flank")


ggplot(corrected_binned_densities_by_sub_type_carriers_norm_for_plot, aes(x=bin_index, y=snv_density,fill=sub_type)) + 
  geom_bar(stat="identity", position="stack") + 
  geom_vline(xintercept=10.5, linetype="dashed") + 
  geom_vline(xintercept=20.5, linetype="dashed") + 
  guides(fill=guide_legend(title="Substitution Type")) + 
  xlab("Bin") + 
  ylab("SNV Density") + 
  annotate("text", x=5, y=1.1*max_density, label="Left Flank") + 
  annotate("text", x=15, y=1.1*max_density, label="Deletion") + 
  annotate("text", x=25, y=1.1*max_density, label="Right Flank")

#Non-carriers
load("~/Downloads/pcawg_data/del_density/by_substitution/tabular_summed_binned_densities_by_sub_type_non_carriers_norm.RData", verbose=T)
binned_densities_by_sub_type_non_carriers_norm_for_plot = tabular_summed_binned_densities_by_sub_type_non_carriers_norm[,.(snv_density=rowSums(.SD)),.(sub_type=x, bin_index = V2),.SDcols=as.character(seq(22))]

max_density = max(binned_densities_by_sub_type_non_carriers_norm_for_plot$snv_density) 
ggplot(binned_densities_by_sub_type_non_carriers_norm_for_plot, aes(x=bin_index, y=snv_density,fill=sub_type)) + 
  geom_bar(stat="identity") + 
  geom_vline(xintercept=10.5, linetype="dashed") + 
  geom_vline(xintercept=20.5, linetype="dashed") + 
  facet_wrap(~sub_type, nrow=2, ncol=3) + 
  guides(fill=guide_legend(title="Substitution Type")) + 
  xlab("Bin") + 
  ylab("SNV Density") + 
  annotate("text", x=5, y=1.1*max_density, label="Left Flank") + 
  annotate("text", x=15, y=1.1*max_density, label="Deletion") + 
  annotate("text", x=25, y=1.1*max_density, label="Right Flank")


binned_densities_by_sub_type_dt = data.table(cbind(do.call(rbind, binned_densities_by_sub_type), names(binned_densities_by_sub_type) %>% map(~rep(.,30)) %>% unlist, rep(seq(30), 6)))
melted_binned_densities_by_sub_type_dt = melt(variable.factor=F, binned_densities_by_sub_type_dt, id.vars = c("V1180412", "V1180413"), value.name = "snv_density", variable.name = "del")
setnames(melted_binned_densities_by_sub_type_dt, "V1180412", "sub_type")
setnames(melted_binned_densities_by_sub_type_dt, "V1180413", "bin_index")
melted_binned_densities_by_sub_type_dt$snv_density = as.double(melted_binned_densities_by_sub_type_dt$snv_density)


data_for_plot = melted_binned_densities_by_sub_type_dt[, .(snv_density = sum(snv_density)), by=c("del","bin_index", "sub_type")][, total_snv_density := sum(snv_density), by=c("bin_index", "sub_type")][, scaled_snv_density := snv_density / total_snv_density][order(bin_index, -scaled_snv_density),]
ggplot(data_for_plot[,.SD[1:5],by=c("bin_index", "sub_type")], aes(x=as.integer(bin_index), y=scaled_snv_density, fill=reorder(-scaled_snv_density,del), label=scales::percent(round(scaled_snv_density, digits = 2)))) + 
  geom_bar(stat="identity", position="stack") + 
  geom_text(size = 3, position = position_stack(vjust = 0.5)) + 
  guides(fill=FALSE) +
  scale_y_continuous(labels = scales::percent) +
  xlab("Bin") + 
  ylab("% Total Density") +
  facet_wrap(~sub_type, nrow=2)


sort(table(data_for_plot[,.SD[1:5],by="bin_index"]$del), decreasing = T)
data_for_plot[,.SD[1:5],by="bin_index"][, .(mean_snv_density = format(mean(scaled_snv_density), scientific=F),max_snv_density = max(scaled_snv_density),num_hits=.N, del_width=width(filtered_deletion_ranges[del]), num_carriers =  length(which(deletion_carrier_mask[del,] > 0))), by=del][order(-num_hits)]
top_deletions = data_for_plot[sub_type=="c_to_g",.SD[1:5],by="bin_index"][, .(mean_snv_density = format(mean(scaled_snv_density), scientific=F),max_snv_density = max(scaled_snv_density),num_hits=.N, del_width=width(deletion_ranges[del]), num_carriers =  length(which(deletion_carrier_mask[del,] > 0)), chrom = as.character(seqnames(deletion_ranges[del])), start = start(deletion_ranges[del]), end = end(deletion_ranges[del])), by=del][order(-max_snv_density)]
plot.new()
grid.table(top_deletions)