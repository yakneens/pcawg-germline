library(ggplot2)
library(data.table)

#By substitution type
load("~/Downloads/pcawg_data/del_density/by_substitution/tabular_summed_binned_densities_by_sub_type_carriers_norm.RData", verbose=T)
binned_densities_by_sub_type_carriers_norm_for_plot = tabular_summed_binned_densities_by_sub_type_carriers_norm[,.(snv_density=rowSums(.SD)),.(sub_type=x, bin_index = V2),.SDcols=as.character(seq(22))]

max_density = max(binned_densities_by_sub_type_carriers_norm_for_plot$snv_density) 
ggplot(binned_densities_by_sub_type_carriers_norm_for_plot, aes(x=bin_index, y=snv_density,fill=sub_type)) + 
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

ggplot(binned_densities_by_sub_type_carriers_norm_for_plot, aes(x=bin_index, y=snv_density,fill=sub_type)) + 
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

