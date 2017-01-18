#Overall density
load("~/Downloads/pcawg_data/del_density/overall/summed_melted_binned_densities_carriers_dt.RData", verbose=T)
bar_width=0.8

pre_del = summed_melted_binned_densities_carriers_dt[bin_index %in% seq(10), snv_density]
pre_del_t = t.test(pre_del)
pre_del_t$estimate
pre_del_t$estimate - pre_del_t$conf.int[1]

in_del = summed_melted_binned_densities_carriers_dt[bin_index %in% (10 + seq(10)), snv_density]
in_del_t = t.test(in_del)

in_del_t$estimate
in_del_t$estimate - in_del_t$conf.int[1]


post_del = summed_melted_binned_densities_carriers_dt[bin_index %in% (20 + seq(10)), snv_density]
post_del_t = t.test(post_del)
post_del_t$estimate
post_del_t$estimate - post_del_t$conf.int[1]


max_density = max(summed_melted_binned_densities_carriers_dt$snv_density) 

ggplot(summed_melted_binned_densities_carriers_dt, aes(x=bin_index - (bar_width/2), y=snv_density, width=bar_width)) + 
  geom_bar(stat="identity") + 
  geom_vline(xintercept=10.5 - (bar_width/2), linetype="dashed") + 
  geom_vline(xintercept=20.5 - (bar_width/2), linetype="dashed") + 
  xlab("Bin") + 
  ylab("SNV Density") + 
  annotate("text", x=5, y=1.1*max_density, label="Left Flank") + 
  annotate("text", x=15, y=1.1*max_density, label="Deletion") + 
  annotate("text", x=25, y=1.1*max_density, label="Right Flank")



#Corrected densities
corrected_summed_melted_binned_densities_carriers_dt = copy(summed_melted_binned_densities_carriers_dt)
corrected_summed_melted_binned_densities_carriers_dt = corrected_summed_melted_binned_densities_carriers_dt[bin_index > 10 & bin_index <= 20, snv_density := 2*snv_density ]

cor_in_del = corrected_summed_melted_binned_densities_carriers_dt[bin_index %in% (10 + seq(10)), snv_density]
cor_in_del_t = t.test(cor_in_del)
cor_in_del_t$estimate
cor_in_del_t$estimate - cor_in_del_t$conf.int[1]


max_density = max(corrected_summed_melted_binned_densities_carriers_dt$snv_density) 

ggplot(corrected_summed_melted_binned_densities_carriers_dt, aes(x=bin_index - (bar_width/2), y=snv_density)) + 
  geom_bar(stat="identity") + 
  geom_vline(xintercept=10.5 - (bar_width/2), linetype="dashed") + 
  geom_vline(xintercept=20.5 - (bar_width/2), linetype="dashed") + 
  xlab("Bin") + 
  ylab("SNV Density") + 
  annotate("text", x=5, y=1.1*max_density, label="Left Flank") + 
  annotate("text", x=15, y=1.1*max_density, label="Deletion") + 
  annotate("text", x=25, y=1.1*max_density, label="Right Flank")

