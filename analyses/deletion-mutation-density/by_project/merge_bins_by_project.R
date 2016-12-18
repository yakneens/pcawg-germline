library(data.table)
num_bins = 30
load("~/Downloads/pcawg_data/del_density/input_data/deletion_ranges.RData")

#Density bins by project carriers
projects = c("BLCA-US", "BOCA-UK", "BRCA-EU", "BRCA-UK", "BRCA-US", "BTCA-SG", "CESC-US", "CLLE-ES", "CMDI-UK", "COAD-US", "DLBC-US", "EOPC-DE", "ESAD-UK", "GACA-CN", "GBM-US", "HNSC-US", "KICH-US", "KIRC-US", "KIRP-US", "LAML-KR", "LGG-US", "LICA-FR", "LIHC-US", "LINC-JP", "LIRI-JP", "LUAD-US", "LUSC-US", "MALY-DE", "MELA-AU", "ORCA-IN", "OV-AU", "OV-US", "PACA-AU", "PACA-CA", "PAEN-AU", "PBCA-DE", "PRAD-CA", "PRAD-UK", "PRAD-US", "READ-US", "RECA-EU")
density_bins_by_project_carriers = list()
for(i in seq_along(projects)){
  filepath = paste("~/Downloads/pcawg_data/del_density/by_project/exp/density_bins_",projects[i], "_carriers.RData", sep="")
  load(filepath, verbose=T)
  obj_name = paste("density_bins_", projects[i], "_carriers", sep="")
  print(obj_name)
  density_bins_by_project_carriers[[i]] = get(obj_name)
}

rm(list=ls(pattern="density_bins_[A-Z]{2,4}-[A-Z]{2}.*_carriers"))

names(density_bins_by_project_carriers) = projects
save(density_bins_by_project_carriers, file="~/Downloads/pcawg_data/del_density/by_project/density_bins_by_project_carriers.RData")

summed_density_bins_by_project_carriers = matrix(unlist(lapply(density_bins_by_project_carriers, rowSums)), nrow=num_bins, byrow = T)
colnames(summed_density_bins_by_project_carriers) = projects
save(summed_density_bins_by_project_carriers, file="~/Downloads/pcawg_data/del_density/by_project/summed_density_bins_by_project_carriers.RData")

tabular_summed_density_bins_by_project_carriers = data.table(summed_density_bins_by_project_carriers)
tabular_summed_density_bins_by_project_carriers[,bin_index := factor(seq(30))]
save(tabular_summed_density_bins_by_project_carriers, file="~/Downloads/pcawg_data/del_density/by_project/tabular_summed_density_bins_by_project_carriers.RData")

#Normalized Density bins by project carriers
projects = c("BLCA-US", "BOCA-UK", "BRCA-EU", "BRCA-UK", "BRCA-US", "BTCA-SG", "CESC-US", "CLLE-ES", "CMDI-UK", "COAD-US", "DLBC-US", "EOPC-DE", "ESAD-UK", "GACA-CN", "GBM-US", "HNSC-US", "KICH-US", "KIRC-US", "KIRP-US", "LAML-KR", "LGG-US", "LICA-FR", "LIHC-US", "LINC-JP", "LIRI-JP", "LUAD-US", "LUSC-US", "MALY-DE", "MELA-AU", "ORCA-IN", "OV-AU", "OV-US", "PACA-AU", "PACA-CA", "PAEN-AU", "PBCA-DE", "PRAD-CA", "PRAD-UK", "PRAD-US", "READ-US", "RECA-EU")
density_bins_by_project_carriers_norm = list()
for(i in seq_along(projects)){
  filepath = paste("~/Downloads/pcawg_data/del_density/by_project/exp/density_bins_",projects[i], "_carriers.RData", sep="")
  load(filepath, verbose=T)
  obj_name = paste("density_bins_", projects[i], "_carriers", sep="")
  print(obj_name)
  density_bins_by_project_carriers_norm[[i]] = t(t(get(obj_name)) / width(deletion_ranges[colnames(get(obj_name))]))
}

rm(list=ls(pattern="density_bins_[A-Z]{2,4}-[A-Z]{2}.*_carriers"))

names(density_bins_by_project_carriers_norm) = projects
save(density_bins_by_project_carriers_norm, file="~/Downloads/pcawg_data/del_density/by_project/density_bins_by_project_carriers_norm.RData")

summed_density_bins_by_project_carriers_norm = matrix(unlist(lapply(density_bins_by_project_carriers_norm, rowSums)), nrow=num_bins, byrow = T)
colnames(summed_density_bins_by_project_carriers_norm) = projects
save(summed_density_bins_by_project_carriers_norm, file="~/Downloads/pcawg_data/del_density/by_project/summed_density_bins_by_project_carriers_norm.RData")

tabular_summed_density_bins_by_project_carriers_norm = data.table(summed_density_bins_by_project_carriers_norm)
tabular_summed_density_bins_by_project_carriers_norm[,bin_index := factor(seq(30))]
save(tabular_summed_density_bins_by_project_carriers_norm, file="~/Downloads/pcawg_data/del_density/by_project/tabular_summed_density_bins_by_project_carriers_norm.RData")

#Density bins by project non_carriers
projects = c("BLCA-US", "BOCA-UK", "BRCA-EU", "BRCA-UK", "BRCA-US", "BTCA-SG", "CESC-US", "CLLE-ES", "CMDI-UK", "COAD-US", "DLBC-US", "EOPC-DE", "ESAD-UK", "GACA-CN", "GBM-US", "HNSC-US", "KICH-US", "KIRC-US", "KIRP-US", "LAML-KR", "LGG-US", "LICA-FR", "LIHC-US", "LINC-JP", "LIRI-JP", "LUAD-US", "LUSC-US", "MALY-DE", "MELA-AU", "ORCA-IN", "OV-AU", "OV-US", "PACA-AU", "PACA-CA", "PAEN-AU", "PBCA-DE", "PRAD-CA", "PRAD-UK", "PRAD-US", "READ-US", "RECA-EU")
density_bins_by_project_non_carriers = list()
for(i in seq_along(projects)){
  filepath = paste("~/Downloads/pcawg_data/del_density/by_project/density_bins_",projects[i], "_non_carriers.RData", sep="")
  load(filepath, verbose=T)
  obj_name = paste("density_bins_", projects[i], "_non_carriers", sep="")
  print(obj_name)
  density_bins_by_project_non_carriers[[i]] = get(obj_name)
}

rm(list=ls(pattern="density_bins_[A-Z]{2,4}-[A-Z]{2}.*_non_carriers"))

names(density_bins_by_project_non_carriers) = projects
save(density_bins_by_project_non_carriers, file="~/Downloads/pcawg_data/del_density/by_project/density_bins_by_project_non_carriers.RData")

summed_density_bins_by_project_non_carriers = matrix(unlist(lapply(density_bins_by_project_non_carriers, rowSums)), nrow=num_bins, byrow = T)
colnames(summed_density_bins_by_project_non_carriers) = projects
save(summed_density_bins_by_project_non_carriers, file="~/Downloads/pcawg_data/del_density/by_project/summed_density_bins_by_project_non_carriers.RData")

tabular_summed_density_bins_by_project_non_carriers = data.table(summed_density_bins_by_project_non_carriers)
tabular_summed_density_bins_by_project_non_carriers[,bin_index := factor(seq(30))]
save(tabular_summed_density_bins_by_project_non_carriers, file="~/Downloads/pcawg_data/del_density/by_project/tabular_summed_density_bins_by_project_non_carriers.RData")

#Normalized Density bins by project non_carriers
projects = c("BLCA-US", "BOCA-UK", "BRCA-EU", "BRCA-UK", "BRCA-US", "BTCA-SG", "CESC-US", "CLLE-ES", "CMDI-UK", "COAD-US", "DLBC-US", "EOPC-DE", "ESAD-UK", "GACA-CN", "GBM-US", "HNSC-US", "KICH-US", "KIRC-US", "KIRP-US", "LAML-KR", "LGG-US", "LICA-FR", "LIHC-US", "LINC-JP", "LIRI-JP", "LUAD-US", "LUSC-US", "MALY-DE", "MELA-AU", "ORCA-IN", "OV-AU", "OV-US", "PACA-AU", "PACA-CA", "PAEN-AU", "PBCA-DE", "PRAD-CA", "PRAD-UK", "PRAD-US", "READ-US", "RECA-EU")
density_bins_by_project_non_carriers_norm = list()
for(i in seq_along(projects)){
  filepath = paste("~/Downloads/pcawg_data/del_density/by_project/density_bins_",projects[i], "_non_carriers.RData", sep="")
  load(filepath, verbose=T)
  obj_name = paste("density_bins_", projects[i], "_non_carriers", sep="")
  print(obj_name)
  density_bins_by_project_non_carriers_norm[[i]] = t(t(get(obj_name)) / width(deletion_ranges[colnames(get(obj_name))]))
}

rm(list=ls(pattern="density_bins_[A-Z]{2,4}-[A-Z]{2}.*_non_carriers"))

names(density_bins_by_project_non_carriers_norm) = projects
save(density_bins_by_project_non_carriers_norm, file="~/Downloads/pcawg_data/del_density/by_project/density_bins_by_project_non_carriers_norm.RData")

summed_density_bins_by_project_non_carriers_norm = matrix(unlist(lapply(density_bins_by_project_non_carriers_norm, rowSums)), nrow=num_bins, byrow = T)
colnames(summed_density_bins_by_project_non_carriers_norm) = projects
save(summed_density_bins_by_project_non_carriers_norm, file="~/Downloads/pcawg_data/del_density/by_project/summed_density_bins_by_project_non_carriers_norm.RData")

tabular_summed_density_bins_by_project_non_carriers_norm = data.table(summed_density_bins_by_project_non_carriers_norm)
tabular_summed_density_bins_by_project_non_carriers_norm[,bin_index := factor(seq(30))]
save(tabular_summed_density_bins_by_project_non_carriers_norm, file="~/Downloads/pcawg_data/del_density/by_project/tabular_summed_density_bins_by_project_non_carriers_norm.RData")
