library(ggplot2)
library(data.table)
num_bins = 30
bin_width = 100


#All breakpoints
load("~/Downloads/pcawg_data/del_density/at_breakpoints/summed_breakpoint_densities_100_30_200_carriers_start.RData", verbose=T)
load("~/Downloads/pcawg_data/del_density/at_breakpoints/summed_breakpoint_densities_100_30_200_carriers_end.RData", verbose=T)

breakpoint_densities_100_30_200_carriers_for_plot = rbind(data.table(colSums(summed_breakpoint_densities_100_30_200_carriers_start), seq(num_bins), "start"),
                                                     data.table(colSums(summed_breakpoint_densities_100_30_200_carriers_end), seq(num_bins), "end"))

setnames(breakpoint_densities_100_30_200_carriers_for_plot, c("snv_density", "bin_index", "breakpoint_type"))
breakpoint_densities_100_30_200_carriers_for_plot$breakpoint_type = factor(breakpoint_densities_100_30_200_carriers_for_plot$breakpoint_type, levels=c("start", "end"))

breakpoint_densities_100_30_200_carriers_for_plot = melt(as.data.table(lapply(binned_densities, rowSums))[,bin_index:=seq(30)], id.vars = "bin_index", variable.name = "breakpoint_type", value.name = "snv_density")



max_density = max(breakpoint_densities_100_30_200_carriers_for_plot$snv_density) 
midpoint = ceiling(num_bins / 2)
ggplot(breakpoint_densities_100_30_200_carriers_for_plot, aes(x=(bin_index - midpoint) * bin_width, y=snv_density,fill=abs(midpoint - bin_index) * bin_width)) + 
  geom_bar(stat="identity") + 
  facet_grid(.~breakpoint_type) + 
  guides(fill=guide_legend(title="Distance from breakpoint")) + 
  xlab("Distance from breakpoint") + 
  ylab("SNV Density") + 
  annotate("text", x=0, y=1.1*max_density, label="Breakpoint")

#Corrected densities
corrected_breakpoint_densities_100_30_200_carriers_for_plot = copy(breakpoint_densities_100_30_200_carriers_for_plot)
corrected_breakpoint_densities_100_30_200_carriers_for_plot[(breakpoint_type == "start" & bin_index > num_bins/2) | (breakpoint_type == "end" & bin_index <= num_bins/2), snv_density := 2*snv_density]

max_density = max(corrected_breakpoint_densities_100_30_200_carriers_for_plot$snv_density) 
midpoint = ceiling(num_bins / 2)
ggplot(corrected_breakpoint_densities_100_30_200_carriers_for_plot, aes(x=(bin_index - midpoint) * bin_width, y=snv_density,fill=abs(midpoint - bin_index) * bin_width)) + 
  geom_bar(stat="identity") + 
  facet_grid(.~breakpoint_type) + 
  guides(fill=guide_legend(title="Distance from breakpoint")) + 
  xlab("Distance from breakpoint") + 
  ylab("SNV Density") + 
  annotate("text", x=0, y=1.1*max_density, label="Breakpoint")



#Non-carriers
load("~/Downloads/pcawg_data/del_density/at_breakpoints/summed_breakpoint_densities_100_30_200_non_carriers_start.RData", verbose=T)
load("~/Downloads/pcawg_data/del_density/at_breakpoints/summed_breakpoint_densities_100_30_200_non_carriers_end.RData", verbose=T)

breakpoint_densities_100_30_200_non_carriers_for_plot = rbind(data.table(colSums(summed_breakpoint_densities_100_30_200_non_carriers_start), seq(num_bins), "start"),
                                                          data.table(colSums(summed_breakpoint_densities_100_30_200_non_carriers_end), seq(num_bins), "end"))
setnames(breakpoint_densities_100_30_200_non_carriers_for_plot, c("snv_density", "bin_index", "breakpoint_type"))
breakpoint_densities_100_30_200_non_carriers_for_plot$breakpoint_type = factor(breakpoint_densities_100_30_200_non_carriers_for_plot$breakpoint_type, levels=c("start", "end"))


max_density = max(breakpoint_densities_100_30_200_non_carriers_for_plot$snv_density) 
midpoint = ceiling(num_bins / 2)
ggplot(breakpoint_densities_100_30_200_non_carriers_for_plot, aes(x=(bin_index - midpoint) * bin_width, y=snv_density,fill=abs(midpoint - bin_index) * bin_width)) + 
  geom_bar(stat="identity") + facet_grid(.~breakpoint_type) + 
  guides(fill=guide_legend(title="Distance from breakpoint")) + 
  xlab("Distance from breakpoint") + 
  ylab("SNV Density") + 
  annotate("text", x=0, y=1.1*max_density, label="Breakpoint")

#By project
load("~/Downloads/pcawg_data/del_density/at_breakpoints/old/proj/summed_breakpoint_densities_precise_by_project_start.RData", verbose=T)
load("~/Downloads/pcawg_data/del_density/at_breakpoints/old/proj/summed_breakpoint_densities_precise_by_project_end.RData", verbose=T)

breakpoint_densities_precise_by_project_for_plot = rbind(do.call(rbind, lapply(names(summed_breakpoint_densities_precise_by_project_start), function(x){data.table(rowSums(t(summed_breakpoint_densities_precise_by_project_start[[x]])),seq(31),x, "start")})),
                                                                do.call(rbind, lapply(names(summed_breakpoint_densities_precise_by_project_end), function(x){data.table(rowSums(t(summed_breakpoint_densities_precise_by_project_end[[x]])),seq(31),x, "end")})))
setnames(breakpoint_densities_precise_by_project_for_plot, c("snv_density", "bin_index", "project_name", "breakpoint_type"))
breakpoint_densities_precise_by_project_for_plot$breakpoint_type = factor(breakpoint_densities_precise_by_project_for_plot$breakpoint_type, levels=c("start", "end"))

max_density = max(breakpoint_densities_precise_by_project_for_plot$snv_density) 
midpoint = ceiling(num_bins / 2)
ggplot(breakpoint_densities_precise_by_project_for_plot, aes(x=(bin_index - midpoint) * bin_width, y=snv_density)) + 
  geom_bar(stat="identity", aes(fill=abs(midpoint - bin_index) * bin_width)) + 
  facet_grid(project_name~breakpoint_type) + 
  guides(fill=guide_legend(title="Distance from breakpoint")) + 
  xlab("Distance from breakpoint") + 
  ylab("SNV Density") + 
  annotate("text", x=0, y=1.1*max_density, label="Breakpoint")

#Corrected densities
corrected_breakpoint_densities_precise_by_project_for_plot = copy(breakpoint_densities_precise_by_project_for_plot)
corrected_breakpoint_densities_precise_by_project_for_plot[(breakpoint_type == "start" & bin_index > 15) | (breakpoint_type == "end" & bin_index < 16), snv_density := 2*snv_density]

max_density = max(corrected_breakpoint_densities_precise_by_project_for_plot[project_name=="PBCA_DE"]$snv_density) 
midpoint = ceiling(num_bins / 2)
ggplot(corrected_breakpoint_densities_precise_by_project_for_plot[project_name=="PBCA_DE"], aes(x=(bin_index - midpoint) * bin_width, y=snv_density)) + 
  geom_bar(stat="identity", aes(fill=abs(midpoint - bin_index) * bin_width)) + 
  facet_grid(.~breakpoint_type) + 
  guides(fill=guide_legend(title="Distance from breakpoint")) + 
  xlab("Distance from breakpoint") + 
  ylab("SNV Density") + annotate("text", x=0, y=1.1*max_density, label="Breakpoint")




#All

load("~/Downloads/pcawg_data/del_density/at_breakpoints/old/proj/summed_breakpoint_densities_all_by_project_start.RData", verbose=T)
load("~/Downloads/pcawg_data/del_density/at_breakpoints/old/proj/summed_breakpoint_densities_all_by_project_end.RData", verbose=T)

breakpoint_densities_all_by_project_for_plot = rbind(do.call(rbind, lapply(names(summed_breakpoint_densities_all_by_project_start), function(x){data.table(rowSums(t(summed_breakpoint_densities_all_by_project_start[[x]])),seq(31),x, "start")})),
                                                         do.call(rbind, lapply(names(summed_breakpoint_densities_all_by_project_end), function(x){data.table(rowSums(t(summed_breakpoint_densities_all_by_project_end[[x]])),seq(31),x, "end")})))
setnames(breakpoint_densities_all_by_project_for_plot, c("snv_density", "bin_index", "project_name", "breakpoint_type"))
breakpoint_densities_all_by_project_for_plot$breakpoint_type = factor(breakpoint_densities_all_by_project_for_plot$breakpoint_type, levels=c("start", "end"))

max_density = max(breakpoint_densities_all_by_project_for_plot$snv_density) 
midpoint = ceiling(num_bins / 2)
ggplot(breakpoint_densities_all_by_project_for_plot, aes(x=(bin_index - midpoint) * bin_width, y=snv_density)) + 
  geom_bar(stat="identity", aes(fill=abs(midpoint - bin_index) * bin_width)) + 
  facet_grid(project_name~breakpoint_type) + 
  guides(fill=guide_legend(title="Distance from breakpoint")) + 
  xlab("Distance from breakpoint") + 
  ylab("SNV Density") + 
  annotate("text", x=0, y=1.1*max_density, label="Breakpoint")

#Corrected densities
corrected_breakpoint_densities_all_by_project_for_plot = copy(breakpoint_densities_all_by_project_for_plot)
corrected_breakpoint_densities_all_by_project_for_plot[(breakpoint_type == "start" & bin_index > 15) | (breakpoint_type == "end" & bin_index < 16), snv_density := 2*snv_density]

max_density = max(corrected_breakpoint_densities_all_by_project_for_plot$snv_density) 
midpoint = ceiling(num_bins / 2)
ggplot(corrected_breakpoint_densities_all_by_project_for_plot, aes(x=(bin_index - midpoint) * bin_width, y=snv_density)) + 
  geom_bar(stat="identity", aes(fill=abs(midpoint - bin_index) * bin_width)) + 
  facet_grid(project_name~breakpoint_type) + 
  guides(fill=guide_legend(title="Distance from breakpoint")) + 
  xlab("Distance from breakpoint") + 
  ylab("SNV Density") + 
  annotate("text", x=0, y=1.1*max_density, label="Breakpoint")








