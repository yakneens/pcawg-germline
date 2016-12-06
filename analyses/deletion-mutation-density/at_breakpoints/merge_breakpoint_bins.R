##Breakpoints carriers
breakpoint_densities_100_30_200_carriers_start = list()
breakpoint_densities_100_30_200_carriers_end = list()
for(i in seq(22)){
  start_filepath = paste("~/Downloads/pcawg_data/del_density/at_breakpoints/start_breakpoint_density_bins_chrom_", i, "_width_100_num_30_margin_200_carriers.RData", sep="")
  load(start_filepath, verbose=T)
  obj_name = paste("start_breakpoint_density_bins_chrom_", i, "_width_100_num_30_margin_200_carriers", sep="")
  breakpoint_densities_100_30_200_carriers_start[[i]] = get(obj_name)
 
  
  end_filepath = paste("~/Downloads/pcawg_data/del_density/at_breakpoints/end_breakpoint_density_bins_chrom_", i, "_width_100_num_30_margin_200_carriers.RData", sep="")
  load(end_filepath, verbose=T)
  obj_name = paste("end_breakpoint_density_bins_chrom_", i, "_width_100_num_30_margin_200_carriers", sep="")
  breakpoint_densities_100_30_200_carriers_end[[i]] = get(obj_name)
  rm(obj_name)
}

rm(list=ls(pattern = ".*start_breakpoint_density_bins.*"))
rm(list=ls(pattern = ".*end_breakpoint_density_bins.*"))

save(breakpoint_densities_100_30_200_carriers_start, file="~/Downloads/pcawg_data/del_density/at_breakpoints/breakpoint_densities_all_100_30_200_carriers_start.RData")
save(breakpoint_densities_100_30_200_carriers_end, file="~/Downloads/pcawg_data/del_density/at_breakpoints/breakpoint_densities_all_100_30_200_carriers_end.RData")

summed_breakpoint_densities_100_30_200_carriers_start = matrix(unlist(lapply(breakpoint_densities_100_30_200_carriers_start, rowSums)), nrow=22, byrow = T)
summed_breakpoint_densities_100_30_200_carriers_end = matrix(unlist(lapply(breakpoint_densities_100_30_200_carriers_end, rowSums)), nrow=22, byrow = T)

save(summed_breakpoint_densities_100_30_200_carriers_start, file="~/Downloads/pcawg_data/del_density/at_breakpoints/summed_breakpoint_densities_100_30_200_carriers_start.RData")
save(summed_breakpoint_densities_100_30_200_carriers_end, file="~/Downloads/pcawg_data/del_density/at_breakpoints/summed_breakpoint_densities_100_30_200_carriers_end.RData")

tabular_summed_breakpoint_densities_100_30_200_carriers_start = data.table(summed_breakpoint_densities_100_30_200_carriers_start)
tabular_summed_breakpoint_densities_100_30_200_carriers_start[,chrom:=seq(22)]

tabular_summed_breakpoint_densities_100_30_200_carriers_end = data.table(summed_breakpoint_densities_100_30_200_carriers_end)
tabular_summed_breakpoint_densities_100_30_200_carriers_end[,chrom:=seq(22)]

save(tabular_summed_breakpoint_densities_100_30_200_carriers_start, file="~/Downloads/pcawg_data/del_density/at_breakpoints/tabular_summed_breakpoint_densities_100_30_200_carriers_start.RData")
save(tabular_summed_breakpoint_densities_100_30_200_carriers_end, file="~/Downloads/pcawg_data/del_density/at_breakpoints/tabular_summed_breakpoint_densities_100_30_200_carriers_end.RData")

##Breakpoints non_carriers
breakpoint_densities_100_30_200_non_carriers_start = list()
breakpoint_densities_100_30_200_non_carriers_end = list()
for(i in seq(22)){
  start_filepath = paste("~/Downloads/pcawg_data/del_density/at_breakpoints/start_breakpoint_density_bins_chrom_", i, "_width_100_num_30_margin_200_non_carriers.RData", sep="")
  load(start_filepath, verbose=T)
  obj_name = paste("start_breakpoint_density_bins_chrom_", i, "_width_100_num_30_margin_200_non_carriers", sep="")
  breakpoint_densities_100_30_200_non_carriers_start[[i]] = get(obj_name)
  
  
  end_filepath = paste("~/Downloads/pcawg_data/del_density/at_breakpoints/end_breakpoint_density_bins_chrom_", i, "_width_100_num_30_margin_200_non_carriers.RData", sep="")
  load(end_filepath, verbose=T)
  obj_name = paste("end_breakpoint_density_bins_chrom_", i, "_width_100_num_30_margin_200_non_carriers", sep="")
  breakpoint_densities_100_30_200_non_carriers_end[[i]] = get(obj_name)
  rm(obj_name)
}

rm(list=ls(pattern = ".*start_breakpoint_density_bins.*"))
rm(list=ls(pattern = ".*end_breakpoint_density_bins.*"))

save(breakpoint_densities_100_30_200_non_carriers_start, file="~/Downloads/pcawg_data/del_density/at_breakpoints/breakpoint_densities_all_100_30_200_non_carriers_start.RData")
save(breakpoint_densities_100_30_200_non_carriers_end, file="~/Downloads/pcawg_data/del_density/at_breakpoints/breakpoint_densities_all_100_30_200_non_carriers_end.RData")

summed_breakpoint_densities_100_30_200_non_carriers_start = matrix(unlist(lapply(breakpoint_densities_100_30_200_non_carriers_start, rowSums)), nrow=22, byrow = T)
summed_breakpoint_densities_100_30_200_non_carriers_end = matrix(unlist(lapply(breakpoint_densities_100_30_200_non_carriers_end, rowSums)), nrow=22, byrow = T)

save(summed_breakpoint_densities_100_30_200_non_carriers_start, file="~/Downloads/pcawg_data/del_density/at_breakpoints/summed_breakpoint_densities_100_30_200_non_carriers_start.RData")
save(summed_breakpoint_densities_100_30_200_non_carriers_end, file="~/Downloads/pcawg_data/del_density/at_breakpoints/summed_breakpoint_densities_100_30_200_non_carriers_end.RData")

tabular_summed_breakpoint_densities_100_30_200_non_carriers_start = data.table(summed_breakpoint_densities_100_30_200_non_carriers_start)
tabular_summed_breakpoint_densities_100_30_200_non_carriers_start[,chrom:=seq(22)]

tabular_summed_breakpoint_densities_100_30_200_non_carriers_end = data.table(summed_breakpoint_densities_100_30_200_non_carriers_end)
tabular_summed_breakpoint_densities_100_30_200_non_carriers_end[,chrom:=seq(22)]

save(tabular_summed_breakpoint_densities_100_30_200_non_carriers_start, file="~/Downloads/pcawg_data/del_density/at_breakpoints/tabular_summed_breakpoint_densities_100_30_200_non_carriers_start.RData")
save(tabular_summed_breakpoint_densities_100_30_200_non_carriers_end, file="~/Downloads/pcawg_data/del_density/at_breakpoints/tabular_summed_breakpoint_densities_100_30_200_non_carriers_end.RData")


#Breakpoints by project all
projects = c("BRCA_US", "PRAD_CA", "PBCA_DE", "LIRI_JP")
summed_breakpoint_densities_all_by_project_start = list()
summed_breakpoint_densities_all_by_project_end = list()
for(i in seq_along(projects)){
  start_filepath = paste("~/Downloads/pcawg_data/del_density/at_breakpoints/proj/breakpoint_densities_100_31_",projects[i], "_all_start.RData", sep="")
  load(start_filepath, verbose=T)
  obj_name = paste("breakpoint_densities_100_31_", projects[i], "_all_start", sep="")
  print(obj_name)
  summed_breakpoint_densities_all_by_project_start[[i]] = matrix(unlist(lapply(get(obj_name), rowSums)), nrow=22, byrow = T)
  
  
  end_filepath = paste("~/Downloads/pcawg_data/del_density/at_breakpoints/proj/breakpoint_densities_100_31_",projects[i], "_all_end.RData", sep="")
  load(end_filepath, verbose=T)
  obj_name = paste("breakpoint_densities_100_31_", projects[i], "_all_end", sep="")
  summed_breakpoint_densities_all_by_project_end[[i]] = matrix(unlist(lapply(get(obj_name), rowSums)), nrow=22, byrow = T)
}

rm(list=ls(pattern="breakpoint_densities_100_31.*"))

names(summed_breakpoint_densities_all_by_project_start) = projects
names(summed_breakpoint_densities_all_by_project_end) = projects

save(summed_breakpoint_densities_all_by_project_start, file="~/Downloads/pcawg_data/del_density/at_breakpoints/proj/summed_breakpoint_densities_all_by_project_start.RData")
save(summed_breakpoint_densities_all_by_project_end, file="~/Downloads/pcawg_data/del_density/at_breakpoints/proj/summed_breakpoint_densities_all_by_project_end.RData")