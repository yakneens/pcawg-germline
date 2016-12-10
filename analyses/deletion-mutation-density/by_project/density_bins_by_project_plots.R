load("~/Downloads/pcawg_data/del_density/by_project/tabular_summed_density_bins_by_project_carriers_norm.RData")
melted_tabular_summed_density_bins_by_project_carriers_norm = melt(tabular_summed_density_bins_by_project_carriers_norm, id.vars = "bin_index", variable.name = "project", value.name = "snv_density")

max_density = max(binned_densities_by_sub_type_non_carriers_norm_for_plot$snv_density) 
ggplot(melted_tabular_summed_density_bins_by_project_carriers_norm, aes(x=bin_index, y=snv_density,fill=project)) + 
  geom_bar(stat="identity") + 
  geom_vline(xintercept=10.5, linetype="dashed") + 
  geom_vline(xintercept=20.5, linetype="dashed") + 
  facet_wrap(~project, ncol=3) + 
  guides(fill=guide_legend(title="Substitution Type")) + 
  xlab("Bin") + 
  ylab("SNV Density") + 
  annotate("text", x=5, y=1.1*max_density, label="Left Flank") + 
  annotate("text", x=15, y=1.1*max_density, label="Deletion") + 
  annotate("text", x=25, y=1.1*max_density, label="Right Flank")

load("~/Downloads/pcawg_data/del_density/by_project/tabular_summed_density_bins_by_project_non_carriers_norm.RData")
a = as.data.table(rowSums(tabular_summed_density_bins_by_project_non_carriers_norm[,1:46]))
a[, bin_index:=seq(30)]
ggplot(a, aes(bin_index, V1)) + 
  geom_bar(stat="identity")
