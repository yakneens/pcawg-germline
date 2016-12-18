#Overall density
load("~/Downloads/pcawg_data/del_density/overall/summed_density_bins_carriers.RData", verbose=T)
num_bins = 30
#binned_densities_carriers_for_plot = data.table(colSums(summed_binned_densities_carriers), seq(num_bins))
binned_densities_carriers_for_plot = data.table(rowSums(binned_densities), seq(num_bins))
setnames(binned_densities_carriers_for_plot, c("snv_density", "bin_index"))

pre_del = binned_densities_carriers_for_plot[bin_index %in% seq(10), snv_density]
pre_del_t = t.test(pre_del)
in_del = binned_densities_carriers_for_plot[bin_index %in% (10 + seq(10)), snv_density]
in_del_t = t.test(in_del)
post_del = binned_densities_carriers_for_plot[bin_index %in% (20 + seq(10)), snv_density]
post_del_t = t.test(post_del)

max_density = max(binned_densities_carriers_for_plot$snv_density) 

ggplot(binned_densities_carriers_for_plot, aes(x=bin_index, y=snv_density)) + geom_bar(stat="identity") + geom_vline(xintercept=10.5, linetype="dashed") + geom_vline(xintercept=20.5, linetype="dashed") + xlab("Bin") + ylab("SNV Density") + annotate("text", x=5, y=1.1*max_density, label="Left Flank") + annotate("text", x=15, y=1.1*max_density, label="Deletion") + annotate("text", x=25, y=1.1*max_density, label="Right Flank")



#Corrected densities
corrected_binned_densities_carriers_for_plot = copy(binned_densities_carriers_for_plot)
corrected_binned_densities_carriers_for_plot = corrected_binned_densities_carriers_for_plot[bin_index > 10 & bin_index <= 20, snv_density := 2*snv_density ]

cor_in_del = corrected_binned_densities_carriers_for_plot[bin_index %in% (10 + seq(10)), snv_density]
cor_in_del_t = t.test(cor_in_del)



max_density = max(corrected_binned_densities_carriers_for_plot$snv_density) 

ggplot(corrected_binned_densities_carriers_for_plot, aes(x=bin_index, y=snv_density)) + geom_bar(stat="identity") + geom_vline(xintercept=10.5, linetype="dashed") + geom_vline(xintercept=20.5, linetype="dashed") + xlab("Bin") + ylab("SNV Density") + annotate("text", x=5, y=1.1*max_density, label="Left Flank") + annotate("text", x=15, y=1.1*max_density, label="Deletion") + annotate("text", x=25, y=1.1*max_density, label="Right Flank")
