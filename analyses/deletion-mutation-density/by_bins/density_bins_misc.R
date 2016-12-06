sample_meta = get_pcawg_metadata("/home/centos/deletion_analysis_data/pcawg_summary.tsv", 
                                 "/home/centos/deletion_analysis_data/excluded_donors.tsv", split_multi_tumors = F)
sample_meta = keep_first_of_multi_tumors(sample_meta)
clinical_meta = get_clinical_metadata("/home/centos/deletion_analysis_data/pcawg_donor_clinical_August2016_v6.tsv", sample_meta)

hist_meta = get_histology_metadata("/home/centos/deletion_analysis_data/pcawg_specimen_histology_August2016_v6.tsv", sample_meta)

load("~/Downloads/pcawg_data/del_density/input_data/deletion_ranges.RData")
load("~/Downloads/pcawg_data/del_density/input_data/donor_meta.RData")
load("~/Downloads/pcawg_data/del_density/input_data/deletion_carrier_mask.RData")
load("~/Downloads/pcawg_data/del_density/input_data/snv_ranges.RData")
load("/home/centos/deletion_analysis_data/dels.Rdata")

#load("~/Downloads/pcawg_data/germline_deletions/dels_chr22.Rdata")

set_deletion_range_ends(deletions)

deletion_genotypes = geno(deletions)$GT[,match(donor_meta$normal_wgs_aliquot_id, colnames(geno(deletions)$GT))] 


#Load somatic SNV data
load("/home/centos/deletion_analysis_data/snv_samples.Rdata")

breakpoint_densities_100_31 = list()
for(i in seq(22)){
  filepath = paste("~/Downloads/pcawg_data/del_density/at_breakpoints/breakpoint_density_bins_chrom_", i, "_width_100_num_31.RData", sep="")
  load(filepath, verbose=T)
  obj_name = paste("breakpoint_density_bins_chrom_", i, "_width_100_num_31", sep="")
  mat_obj = matrix(get(obj_name), nrow=31)
  print(dim(mat_obj))
  breakpoint_densities_100_31[[i]] = mat_obj
}

#breakpoint_densities = matrix(unlist(lapply(breakpoint_densities_100_31, rowSums)), nrow=22, byrow = T)
breakpoint_densities = matrix(unlist(lapply(breakpoint_densities_BRCA_US, rowSums)), nrow=22, byrow = T)
breakpoint_densities = matrix(unlist(lapply(breakpoint_densities_PRAD_CA, rowSums)), nrow=22, byrow = T)
breakpoint_densities = matrix(unlist(lapply(breakpoint_densities_PBCA_DE, rowSums)), nrow=22, byrow = T)
breakpoint_densities = matrix(unlist(lapply(breakpoint_densities_LIRI_JP, rowSums)), nrow=22, byrow = T)
#breakpoint_densities = matrix(unlist(lapply(breakpoint_densities_100_31_size_10kb_to_100kb, rowSums)), nrow=22, byrow = T)
del_bins = colSums(breakpoint_densities)
del_bins = data.table(seq(1,31), del_bins)

setnames(del_bins, c("index", "bins"))

max_density = max(del_bins$bins) 
midpoint = ceiling(num_bins / 2)
ggplot(del_bins, aes(x=(index - midpoint) * bin_width, y=bins)) + geom_bar(stat="identity", aes(fill=abs(midpoint - index) * bin_width))  + guides(fill=guide_legend(title="Distance from breakpoint")) + xlab("Distance from breakpoint") + ylab("SNV Density") + annotate("text", x=0, y=1.1*max_density, label="Breakpoint")


del_bins2 = copy(del_bins)
del_bins2[16:31, bins := bins*2]
max_density = max(del_bins2$bins) 
midpoint = ceiling(num_bins / 2)
ggplot(del_bins2, aes(x=(index - midpoint) * bin_width, y=bins)) + geom_bar(stat="identity", aes(fill=abs(midpoint - index) * bin_width))  + guides(fill=guide_legend(title="Distance from breakpoint")) + xlab("Distance from breakpoint") + ylab("SNV Density") + annotate("text", x=0, y=1.1*max_density, label="Breakpoint")

projects = unique(donor_meta$dcc_project_code)
project_densities = list()
for(project in projects){
  filepath = paste("~/Downloads/pcawg_data/del_density/by_project/density_bins_", project, ".RData", sep="")
  if(file.exists(filepath)){
    load(filepath, verbose=T)
    obj_name = paste("density_bins_", project, sep="")
    mat_obj = matrix(get(obj_name), nrow=30)
    project_densities[[project]] = mat_obj
  }
}

sample_counts_by_project = donor_meta[dcc_project_code %in% names(project_densities),.N, by=dcc_project_code]

aggregated_project_densities = data.table(do.call(rbind, lapply(project_densities, rowSums)))
aggregated_project_densities[,project := names(project_densities)]
aggregated_project_densities[, sample_count := sample_counts_by_project$N]

melted_densities = melt(aggregated_project_densities, id.vars = c("project", "sample_count"))
melted_densities = melted_densities[order(project),]
melted_densities[, index:=seq(30)]
melted_densities[, normalized_value := value / sample_count]


max_density = max(melted_densities[1:90]$normalized_value) 
ggplot(melted_densities, aes(x=index, y=value/sample_count)) + geom_bar(stat="identity", fill="grey70", aes(colour=project)) + geom_vline(xintercept=10.5, linetype="dashed") + geom_vline(xintercept=20.5, linetype="dashed") + xlab("Bin") + ylab("SNV Density") + annotate("text", x=5, y=1.1*max_density, label="Left Flank") + annotate("text", x=15, y=1.1*max_density, label="Deletion") + annotate("text", x=25, y=1.1*max_density, label="Right Flank")

ggplot(melted_densities, aes(x=index, y=normalized_value)) + geom_line(aes(colour=project))


ggplot(melted_densities) + geom_bar(stat="identity", fill="grey70", aes(x=index, y=value, colour=project))

pdf(file="~/Downloads/densities_by_project.pdf")
for(my_project in names(project_densities)){
  max_density = max(melted_densities[project == my_project]$normalized_value) 
  print(ggplot(melted_densities[project == my_project], aes(x=index, y=normalized_value)) + geom_bar(stat="identity", fill="grey70") + geom_vline(xintercept=10.5, linetype="dashed") + geom_vline(xintercept=20.5, linetype="dashed") + xlab("Bin") + ylab("SNV Density") + annotate("text", x=5, y=1.1*max_density, label="Left Flank") + annotate("text", x=15, y=1.1*max_density, label="Deletion") + annotate("text", x=25, y=1.1*max_density, label="Right Flank") + annotate("text", x=15, y=0.1*max_density, label=paste("Project: ", my_project, sep="")))
}
dev.off()

corrected_melted_densities = melted_densities
corrected_melted_densities[index > 10 & index <= 20, value := 2 * value]
corrected_melted_densities[, normalized_value := value / sample_count]
pdf(file="~/Downloads/corrected_densities_by_project.pdf")
for(my_project in names(project_densities)){
  max_density = max(corrected_melted_densities[project == my_project]$normalized_value) 
  print(ggplot(corrected_melted_densities[project == my_project], aes(x=index, y=normalized_value)) + geom_bar(stat="identity", fill="grey70") + geom_vline(xintercept=10.5, linetype="dashed") + geom_vline(xintercept=20.5, linetype="dashed") + xlab("Bin") + ylab("SNV Density") + annotate("text", x=5, y=1.1*max_density, label="Left Flank") + annotate("text", x=15, y=1.1*max_density, label="Deletion") + annotate("text", x=25, y=1.1*max_density, label="Right Flank") + annotate("text", x=15, y=0.1*max_density, label=paste("Project: ", my_project, sep="")))
}
dev.off()