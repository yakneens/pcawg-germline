#Overall density
binned_densities_carriers = list()
for(i in seq(22)){
  filepath = paste("~/Downloads/pcawg_data/del_density/overall/density_bins_", i, "_carriers.RData", sep="")
  load(filepath, verbose=T)
  obj_name = paste("density_bins_", i, "_carriers", sep="")
  binned_densities_carriers[[i]] = get(obj_name)
  rm(obj_name)
}
rm(list=ls(pattern = "density_bins.*carriers"))

save(binned_densities_carriers, file="~/Downloads/pcawg_data/del_density/overall/density_bins_carriers.RData")

summed_binned_densities_carriers = matrix(unlist(lapply(binned_densities_carriers, rowSums)), nrow=22, byrow = T)
save(summed_binned_densities_carriers, file="~/Downloads/pcawg_data/del_density/overall/summed_density_bins_carriers.RData")

#Normalized overall density
binned_densities_carriers_norm = list()
for(i in seq(22)){
  filepath = paste("~/Downloads/pcawg_data/del_density/overall/density_bins_", i, "_carriers.RData", sep="")
  load(filepath, verbose=T)
  obj_name = paste("density_bins_", i, "_carriers", sep="")
  binned_densities_carriers_norm[[i]] = t(t(get(obj_name)) / width(deletion_ranges[colnames(get(obj_name))]))
  
}
rm(list=ls(pattern = "density_bins.*carriers"))

save(binned_densities_carriers_norm, file="~/Downloads/pcawg_data/del_density/overall/density_bins_carriers_norm.RData")

summed_binned_densities_carriers_norm = matrix(unlist(lapply(binned_densities_carriers_norm, rowSums)), nrow=22, byrow = T)
save(summed_binned_densities_carriers_norm, file="~/Downloads/pcawg_data/del_density/overall/summed_density_bins_carriers_norm.RData")

