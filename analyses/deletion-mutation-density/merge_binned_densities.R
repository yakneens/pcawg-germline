##Precise breakpoints
breakpoint_densities_precise_100_31_start = list()
breakpoint_densities_precise_100_31_end = list()
for(i in seq(22)){
  filepath = paste("~/Downloads/pcawg_data/del_density/at_breakpoints/precise/end_breakpoint_density_bins_chrom_", i, "_width_100_num_31.RData", sep="")
  load(filepath, verbose=T)
  obj_name = paste("end_breakpoint_density_bins_chrom_", i, "_width_100_num_31", sep="")
  mat_obj = matrix(get(obj_name), nrow=31)
  breakpoint_densities_precise_100_31_end[[i]] = mat_obj
  rm(obj_name)
}
rm(list=ls(pattern = "*end_breakpoint_density_bins*"))

save(breakpoint_densities_precise_100_31_start, file="~/Downloads/pcawg_data/del_density/at_breakpoints/precise/breakpoint_densities_precise_100_31_start.RData")
save(breakpoint_densities_precise_100_31_end, file="~/Downloads/pcawg_data/del_density/at_breakpoints/precise/breakpoint_densities_precise_100_31_end.RData")

summed_breakpoint_densities_precise_100_31_start = matrix(unlist(lapply(breakpoint_densities_precise_100_31_start, rowSums)), nrow=22, byrow = T)
summed_breakpoint_densities_precise_100_31_end = matrix(unlist(lapply(breakpoint_densities_precise_100_31_end, rowSums)), nrow=22, byrow = T)
rm(list=ls(pattern="breakpoint_densities_precise_100_31_*"))

save(summed_breakpoint_densities_precise_100_31_start, file="~/Downloads/pcawg_data/del_density/at_breakpoints/precise/summed_breakpoint_densities_precise_100_31_start.RData")
save(summed_breakpoint_densities_precise_100_31_end, file="~/Downloads/pcawg_data/del_density/at_breakpoints/precise/summed_breakpoint_densities_precise_100_31_end.RData")

##All breakpoints
breakpoint_densities_all_100_30_start = list()
breakpoint_densities_all_100_30_end = list()
for(i in seq(22)){
  filepath = paste("~/Downloads/pcawg_data/del_density/at_breakpoints/precise_and_imprecise/end_breakpoint_density_bins_chrom_", i, "_width_100_num_30.RData", sep="")
  load(filepath, verbose=T)
  obj_name = paste("end_breakpoint_density_bins_chrom_", i, "_width_100_num_30", sep="")
  mat_obj = matrix(get(obj_name), nrow=30)
  breakpoint_densities_all_100_30_end[[i]] = mat_obj
  rm(obj_name)
}
rm(list=ls(pattern = "*end_breakpoint_density_bins*"))

save(breakpoint_densities_all_100_30_start, file="~/Downloads/pcawg_data/del_density/at_breakpoints/precise_and_imprecise/breakpoint_densities_all_100_30_start.RData")
save(breakpoint_densities_all_100_30_end, file="~/Downloads/pcawg_data/del_density/at_breakpoints/precise_and_imprecise/breakpoint_densities_all_100_30_end.RData")

summed_breakpoint_densities_all_100_30_start = matrix(unlist(lapply(breakpoint_densities_all_100_30_start, rowSums)), nrow=22, byrow = T)
summed_breakpoint_densities_all_100_30_end = matrix(unlist(lapply(breakpoint_densities_all_100_30_end, rowSums)), nrow=22, byrow = T)

save(summed_breakpoint_densities_all_100_30_start, file="~/Downloads/pcawg_data/del_density/at_breakpoints/precise_and_imprecise/summed_breakpoint_densities_all_100_30_start.RData")
save(summed_breakpoint_densities_all_100_30_end, file="~/Downloads/pcawg_data/del_density/at_breakpoints/precise_and_imprecise/summed_breakpoint_densities_all_100_30_end.RData")

#Breakpoints by project precise
projects = c("BRCA_US", "PRAD_CA", "PBCA_DE", "LIRI_JP")
summed_breakpoint_densities_precise_by_project_start = list()
summed_breakpoint_densities_precise_by_project_end = list()
for(i in seq_along(projects)){
  start_filepath = paste("~/Downloads/pcawg_data/del_density/at_breakpoints/proj/breakpoint_densities_100_31_",projects[i], "_start.RData", sep="")
  load(start_filepath, verbose=T)
  obj_name = paste("breakpoint_densities_100_31_", projects[i], "_start", sep="")
  print(obj_name)
  summed_breakpoint_densities_precise_by_project_start[[i]] = matrix(unlist(lapply(get(obj_name), rowSums)), nrow=22, byrow = T)
  

  end_filepath = paste("~/Downloads/pcawg_data/del_density/at_breakpoints/proj/breakpoint_densities_100_31_",projects[i], "_end.RData", sep="")
  load(end_filepath, verbose=T)
  obj_name = paste("breakpoint_densities_100_31_", projects[i], "_end", sep="")
  summed_breakpoint_densities_precise_by_project_end[[i]] = matrix(unlist(lapply(get(obj_name), rowSums)), nrow=22, byrow = T)
}

rm(list=ls(pattern="breakpoint_densities_100_31*"))

names(summed_breakpoint_densities_precise_by_project_start) = projects
names(summed_breakpoint_densities_precise_by_project_end) = projects

save(summed_breakpoint_densities_precise_by_project_start, file="~/Downloads/pcawg_data/del_density/at_breakpoints/proj/summed_breakpoint_densities_precise_by_project_start.RData")
save(summed_breakpoint_densities_precise_by_project_end, file="~/Downloads/pcawg_data/del_density/at_breakpoints/proj/summed_breakpoint_densities_precise_by_project_end.RData")

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

rm(list=ls(pattern="breakpoint_densities_100_31*"))

names(summed_breakpoint_densities_all_by_project_start) = projects
names(summed_breakpoint_densities_all_by_project_end) = projects

save(summed_breakpoint_densities_all_by_project_start, file="~/Downloads/pcawg_data/del_density/at_breakpoints/proj/summed_breakpoint_densities_all_by_project_start.RData")
save(summed_breakpoint_densities_all_by_project_end, file="~/Downloads/pcawg_data/del_density/at_breakpoints/proj/summed_breakpoint_densities_all_by_project_end.RData")

#Density by substitution type precise
sub_types = c("c_to_a", "c_to_t", "c_to_g", "t_to_a", "t_to_c", "t_to_g")
binned_densities_by_sub_type_precise = list()
for(j in seq_along(sub_types)){
  binned_densities_by_sub_type_precise[[j]] = list()
  for(i in seq(22)){
    filepath = paste("~/Downloads/pcawg_data/del_density/by_substitution/precise/density_bins_chrom_", i, "_", sub_types[j], ".RData", sep="")
    load(filepath, verbose=T)
    obj_name = paste("density_bins_chrom_", i, "_", sub_types[j], sep="")
    binned_densities_by_sub_type_precise[[j]][[i]] = get(obj_name)[[1]]
  }
  names(binned_densities_by_sub_type_precise[[j]]) = seq(22)
}
names(binned_densities_by_sub_type_precise) = sub_types
rm(list=ls(pattern="density_bins_chrom*"))

save(binned_densities_by_sub_type_precise, file="~/Downloads/pcawg_data/del_density/by_substitution/precise/bined_densitites_by_sub_type_precise.RData")

summed_binned_densities_by_sub_type_precise = lapply(binned_densities_by_sub_type_precise, function(x){do.call(rbind, lapply(x, rowSums));})
save(summed_binned_densities_by_sub_type_precise, file="~/Downloads/pcawg_data/del_density/by_substitution/precise/summed_binned_densities_by_sub_type_precise.RData")

tabular_summed_binned_densities_by_sub_type_precise = do.call(rbind, lapply(sub_types, function(x){data.table(t(summed_binned_densities_by_sub_type_precise[[x]]),seq(30),x)}))

#Density by substitution type all
sub_types = c("c_to_a", "c_to_t", "c_to_g", "t_to_a", "t_to_c", "t_to_g")
binned_densities_by_sub_type_all = list()
for(j in seq_along(sub_types)){
  binned_densities_by_sub_type_all[[j]] = list()
  for(i in seq(22)){
    filepath = paste("~/Downloads/pcawg_data/del_density/by_substitution/precise_and_imprecise/density_bins_chrom_", i, "_", sub_types[j], "_all.RData", sep="")
    load(filepath, verbose=T)
    obj_name = paste("density_bins_chrom_", i, "_", sub_types[j], "_all", sep="")
    binned_densities_by_sub_type_all[[j]][[i]] = get(obj_name)[[1]]
  }
  names(binned_densities_by_sub_type_all[[j]]) = seq(22)
}
names(binned_densities_by_sub_type_all) = sub_types
rm(list=ls(pattern="density_bins_chrom*"))

save(binned_densities_by_sub_type_all, file="~/Downloads/pcawg_data/del_density/by_substitution/precise_and_imprecise/bined_densitites_by_sub_type_all.RData")

summed_binned_densities_by_sub_type_all = lapply(binned_densities_by_sub_type_all, function(x){do.call(rbind, lapply(x, rowSums));})
save(summed_binned_densities_by_sub_type_all, file="~/Downloads/pcawg_data/del_density/by_substitution/precise_and_imprecise/summed_binned_densities_by_sub_type_all.RData")

tabular_summed_binned_densities_by_sub_type_all = do.call(rbind, lapply(sub_types, function(x){data.table(t(summed_binned_densities_by_sub_type_all[[x]]),seq(30),x)}))
save(tabular_summed_binned_densities_by_sub_type_all, file="~/Downloads/pcawg_data/del_density/by_substitution/precise_and_imprecise/tabular_summed_binned_densities_by_sub_type_all.RData")

