#Density by substitution type all
sub_types = c("c_to_a", "c_to_t", "c_to_g", "t_to_a", "t_to_c", "t_to_g")
binned_densities_by_sub_type_carriers_all = list()
for(j in seq_along(sub_types)){
  binned_densities_by_sub_type_carriers_all[[j]] = list()
  for(i in seq(22)){
    filepath = paste("~/Downloads/pcawg_data/del_density/by_substitution/density_bins_chrom_", i, "_", sub_types[j], "_carriers_all.RData", sep="")
    load(filepath, verbose=T)
    obj_name = paste("density_bins_chrom_", i, "_", sub_types[j], "_carriers_all", sep="")
    binned_densities_by_sub_type_carriers_all[[j]][[i]] = get(obj_name)[[1]]
  }
  names(binned_densities_by_sub_type_carriers_all[[j]]) = seq(22)
}
names(binned_densities_by_sub_type_carriers_all) = sub_types
rm(list=ls(pattern="density_bins_chrom.*"))

save(binned_densities_by_sub_type_carriers_all, file="~/Downloads/pcawg_data/del_density/by_substitution/binned_densitites_by_sub_type_carriers_all.RData")

summed_binned_densities_by_sub_type_carriers_all = lapply(binned_densities_by_sub_type_carriers_all, function(x){do.call(rbind, lapply(x, rowSums));})
save(summed_binned_densities_by_sub_type_carriers_all, file="~/Downloads/pcawg_data/del_density/by_substitution/summed_binned_densities_by_sub_type_carriers_all.RData")

tabular_summed_binned_densities_by_sub_type_carriers_all = do.call(rbind, lapply(sub_types, function(x){data.table(t(summed_binned_densities_by_sub_type_carriers_all[[x]]),seq(30),x)}))
save(tabular_summed_binned_densities_by_sub_type_carriers_all, file="~/Downloads/pcawg_data/del_density/by_substitution/tabular_summed_binned_densities_by_sub_type_carriers_all.RData")

#Normalized density by substitution type
sub_types = c("c_to_a", "c_to_t", "c_to_g", "t_to_a", "t_to_c", "t_to_g")
binned_densities_by_sub_type_carriers_norm = list()
for(j in seq_along(sub_types)){
  binned_densities_by_sub_type_carriers_norm[[j]] = list()
  for(i in seq(22)){
    filepath = paste("~/Downloads/pcawg_data/del_density/by_substitution/density_bins_chrom_", i, "_", sub_types[j], "_carriers_all.RData", sep="")
    load(filepath, verbose=T)
    obj_name = paste("density_bins_chrom_", i, "_", sub_types[j], "_carriers_all", sep="")
    binned_densities_by_sub_type_carriers_norm[[j]][[i]] = t(t(get(obj_name)[[1]]) / width(deletion_ranges[colnames(get(obj_name)[[1]])]))
  }
  names(binned_densities_by_sub_type_carriers_norm[[j]]) = seq(22)
}
names(binned_densities_by_sub_type_carriers_norm) = sub_types
rm(list=ls(pattern="density_bins_chrom.*"))

save(binned_densities_by_sub_type_carriers_norm, file="~/Downloads/pcawg_data/del_density/by_substitution/binned_densitites_by_sub_type_carriers_norm.RData")

summed_binned_densities_by_sub_type_carriers_norm = lapply(binned_densities_by_sub_type_carriers_norm, function(x){do.call(rbind, lapply(x, rowSums));})
save(summed_binned_densities_by_sub_type_carriers_norm, file="~/Downloads/pcawg_data/del_density/by_substitution/summed_binned_densities_by_sub_type_carriers_norm.RData")

tabular_summed_binned_densities_by_sub_type_carriers_norm = do.call(rbind, lapply(sub_types, function(x){data.table(t(summed_binned_densities_by_sub_type_carriers_norm[[x]]),seq(30),x)}))
save(tabular_summed_binned_densities_by_sub_type_carriers_norm, file="~/Downloads/pcawg_data/del_density/by_substitution/tabular_summed_binned_densities_by_sub_type_carriers_norm.RData")

c("BLCA-US","BOCA-UK","BRCA-EU","BRCA-UK","BRCA-US","BTCA-SG","CESC-US","CLLE-ES","CMDI-UK","COAD-US","DLBC-US","EOPC-DE","ESAD-UK","GACA-CN","GBM-US","HNSC-US","KICH-US","KIRC-US","KIRP-US","LAML-KR","LGG-US","LICA-FR","LIHC-US","LINC-JP","LIRI-JP","LUAD-US","LUSC-US","MALY-DE","MELA-AU","ORCA-IN","OV-AU","OV-US","PACA-AU","PACA-CA","PAEN-AU","PBCA-DE","PRAD-CA","PRAD-UK","PRAD-US","READ-US","RECA-EU","SARC-US","SKCM-US","STAD-US","THCA-US","UCEC-US")


#Density by substitution type all
sub_types = c("c_to_a", "c_to_t", "c_to_g", "t_to_a", "t_to_c", "t_to_g")
binned_densities_by_sub_type_non_carriers_all = list()
for(j in seq_along(sub_types)){
  binned_densities_by_sub_type_non_carriers_all[[j]] = list()
  for(i in seq(22)){
    filepath = paste("~/Downloads/pcawg_data/del_density/by_substitution/density_bins_chrom_", i, "_", sub_types[j], "_non_carriers_all.RData", sep="")
    load(filepath, verbose=T)
    obj_name = paste("density_bins_chrom_", i, "_", sub_types[j], "_non_carriers_all", sep="")
    binned_densities_by_sub_type_non_carriers_all[[j]][[i]] = get(obj_name)[[1]]
    rm(list=ls(pattern="density_bins_chrom.*"))
  }
  names(binned_densities_by_sub_type_non_carriers_all[[j]]) = seq(22)
  
}
names(binned_densities_by_sub_type_non_carriers_all) = sub_types


save(binned_densities_by_sub_type_non_carriers_all, file="~/Downloads/pcawg_data/del_density/by_substitution/binned_densitites_by_sub_type_non_carriers_all.RData")

summed_binned_densities_by_sub_type_non_carriers_all = lapply(binned_densities_by_sub_type_non_carriers_all, function(x){do.call(rbind, lapply(x, rowSums));})
save(summed_binned_densities_by_sub_type_non_carriers_all, file="~/Downloads/pcawg_data/del_density/by_substitution/summed_binned_densities_by_sub_type_non_carriers_all.RData")

tabular_summed_binned_densities_by_sub_type_non_carriers_all = do.call(rbind, lapply(sub_types, function(x){data.table(t(summed_binned_densities_by_sub_type_non_carriers_all[[x]]),seq(30),x)}))
save(tabular_summed_binned_densities_by_sub_type_non_carriers_all, file="~/Downloads/pcawg_data/del_density/by_substitution/tabular_summed_binned_densities_by_sub_type_non_carriers_all.RData")

#Normalized density by substitution type
sub_types = c("c_to_a", "c_to_t", "c_to_g", "t_to_a", "t_to_c", "t_to_g")
binned_densities_by_sub_type_carriers_norm = list()
for(j in seq_along(sub_types)){
  binned_densities_by_sub_type_carriers_norm[[j]] = list()
  for(i in seq(22)){
    filepath = paste("~/Downloads/pcawg_data/del_density/by_substitution/density_bins_chrom_", i, "_", sub_types[j], "_carriers_all.RData", sep="")
    load(filepath, verbose=T)
    obj_name = paste("density_bins_chrom_", i, "_", sub_types[j], "_carriers_all", sep="")
    binned_densities_by_sub_type_carriers_norm[[j]][[i]] = t(t(get(obj_name)[[1]]) / width(deletion_ranges[colnames(get(obj_name)[[1]])]))
  }
  names(binned_densities_by_sub_type_carriers_norm[[j]]) = seq(22)
  rm(list=ls(pattern="density_bins_chrom.*"))
}
names(binned_densities_by_sub_type_carriers_norm) = sub_types


save(binned_densities_by_sub_type_carriers_norm, file="~/Downloads/pcawg_data/del_density/by_substitution/binned_densitites_by_sub_type_carriers_norm.RData")

summed_binned_densities_by_sub_type_carriers_norm = lapply(binned_densities_by_sub_type_carriers_norm, function(x){do.call(rbind, lapply(x, rowSums));})
save(summed_binned_densities_by_sub_type_carriers_norm, file="~/Downloads/pcawg_data/del_density/by_substitution/summed_binned_densities_by_sub_type_carriers_norm.RData")

tabular_summed_binned_densities_by_sub_type_carriers_norm = do.call(rbind, lapply(sub_types, function(x){data.table(t(summed_binned_densities_by_sub_type_carriers_norm[[x]]),seq(30),x)}))
save(tabular_summed_binned_densities_by_sub_type_carriers_norm, file="~/Downloads/pcawg_data/del_density/by_substitution/tabular_summed_binned_densities_by_sub_type_carriers_norm.RData")
