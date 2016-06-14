library(VariantAnnotation)
library(data.table)
library(pcawg.common)

sample_meta = get_pcawg_metadata("~/Downloads/pcawg_data/sample_metadata/pcawg_summary.tsv")

load("~/Downloads/pcawg_data/germline_deletions/dels.Rdata")
load("~/Downloads/pcawg_data/snv_samples.Rdata")



# Normal samples that deletions have been called in
normal_wgs_aliquots = colnames(dels)

# Select donor_unique_id, normal_wgs_aliquot_id, tumor_wgs_aliquot_id from sample metadata
sub_meta = sample_meta[which(sample_meta$normal_wgs_aliquot_id %in% normal_wgs_aliquots),c("normal_wgs_aliquot_id", "donor_unique_id", "tumor_wgs_aliquot_id")]

# For samples with mutliple tumors take the first one in the list
sub_meta$tumor_wgs_aliquot_id = sapply(sub_meta$tumor_wgs_aliquot_id, function(x) as.factor(strsplit(as.character(x), ',')[[1]][1]))

# Rename columns to donors in deletions data frame
colnames(dels) = sapply(colnames(dels), function(x) sub_meta[sub_meta$normal_wgs_aliquot_id == x, "donor_unique_id"])

#If genotype of element is in genotypes_list record it as 1, otherwise record it as 0
genotypeMask <- function(element, genotypes_list){
  if(element %in% genotypes_list){
    return(1)
  } else {
    return(0)
  }
}

#Data frame of heterozygous deletion carriers 
het_carriers = as.data.frame(t(apply(geno(dels)$GT, MARGIN=c(1,2), FUN=genotypeMask, genotypes_list=c("0/1", "1/0"))))
het_carriers$donor_unique_id = sub_meta$donor_unique_id

#Data frame of deletion non-carriers (includes ref, hom alt, and no-calls) 
non_carriers = as.data.frame(t(apply(geno(dels)$GT, MARGIN=c(1,2), FUN=genotypeMask, genotypes_list=c("0/0"))))
non_carriers$donor_unique_id = sub_meta$donor_unique_id


snv_samples = list()
snv_hits = data.frame(matrix(ncol=length(rownames(dels)), nrow=0))
hit_donor_ids = list()

# Go through SNV samples patient-by-patient and calculate overlap with deletions
# Store overlaps in snp_hits (rows ar donors, columns are deletions)
for(aliquot in sub_meta$tumor_wgs_aliquot_id){
  tabix_filename = paste("~/Downloads/pcawg_data/somatic_snv_mnv/",aliquot,".annotated.snv_mnv.vcf.gz", sep="")
  if(file.exists(tabix_filename)){
    tab = TabixFile(tabix_filename)
    snv_samples[aliquot] = readVcf(tab, "hs37d5", param=rng)
    
    num_snv = length(snv_samples[[aliquot]])
    
    #List of donor IDs that have SNVs called in the tumor
    hit_donor_ids = c(hit_donor_ids, as.character(sub_meta[sub_meta$tumor_wgs_aliquot_id == aliquot,]$donor_unique_id))
    
    per_del_hits = as.table(findOverlaps(rowRanges(dels), rowRanges(snv_samples[[aliquot]])))
    # Normalize number of hits by total SNV count for sample and by deletion width
    norm_per_del_hits = per_del_hits / num_snv / width(ranges(rowRanges(dels)))
    
    snv_hits = rbind(snv_hits, norm_per_del_hits)
  }
}

colnames(snv_hits) = as.character(rownames(dels))

#Add column of donor IDs to the snp_hits data frame.
snv_hits$donor_unique_id = unlist(hit_donor_ids)

calculate_hits <- function(my_data){
  hit_list = list()
  
  #Last column has donor IDs, so don't consider it
  for(del_name in head(colnames(my_data), n=-1)){
    #Select donors that have the feature of interest
    selection = my_data[which(my_data[, del_name] == 1), "donor_unique_id"]
    
    #Select number of SNP hits falling on this deletion for each donor in selection
    hits = snv_hits[which(snv_hits$donor_unique_id %in% selection), del_name]
    hit_list = c(hit_list, list(hits))
  }
  
  return(hit_list)
}

carrier_hits =  calculate_hits(het_carriers)

#Dosage correction
#carrier_hits = lapply(carrier_hits, function(x) 2*x)

carrier_means = unlist(lapply(carrier_hits, mean))

non_carrier_hits = calculate_hits(non_carriers)
non_carrier_means = unlist(lapply(non_carrier_hits, mean))

t.test(carrier_means, non_carrier_means)
wilcox.test(carrier_means, non_carrier_means)


# Get rid of observations where there were 0 SNVs observed in deletion carriers
dels_with_carrier_hit_indices = which(lapply(carrier_hits, sum) > 0)
non_zero_carrier_hits = carrier_hits[dels_with_carrier_hit_indices]
non_zero_non_carrier_hits = non_carrier_hits[dels_with_carrier_hit_indices]

non_zero_carrier_hits_means = unlist(lapply(non_zero_carrier_hits, mean))
non_zero_non_carrier_hits_means = unlist(lapply(non_zero_non_carrier_hits, mean))

t.test(non_zero_carrier_hits_means, non_zero_non_carrier_hits_means)
wilcox.test(non_zero_carrier_hits_means, non_zero_non_carrier_hits_means)

boxplot(log(non_zero_carrier_hits_means), log(non_zero_non_carrier_hits_means))

par(mfrow=c(1,2))
hist(log(carrier_means), col=rgb(0.8,0.8,0.8,0.5), ylim=c(0,250))
hist(log(non_carrier_means), col=rgb(0.1,0.1,0.1,0.5),  add=T) 

hist(log(non_zero_carrier_hits_means), col=rgb(0.8,0.8,0.8,0.5), ylim=c(0, 50))
hist(log(non_zero_non_carrier_hits_means), col=rgb(0.1,0.1,0.1,0.5),  add=T) 

max_carrier_density_index = which(lapply(carrier_hits, max) == max(unlist(lapply(carrier_hits, max))))
max(unlist(carrier_hits[max_carrier_density_index]))
max(unlist(non_carrier_hits[max_carrier_density_index]))

max_non_carrier_density_index = which(lapply(non_carrier_hits, max) == max(unlist(lapply(non_carrier_hits, max))))
max(unlist(carrier_hits[max_non_carrier_density_index]))
max(unlist(non_carrier_hits[max_non_carrier_density_index]))

max_mean_carrier_density_index = which(carrier_means == max(carrier_means, na.rm=TRUE))
carrier_means[max_mean_carrier_density_index]
non_carrier_means[max_mean_carrier_density_index]

max_mean_non_carrier_density_index = which(non_carrier_means == max(non_carrier_means, na.rm=TRUE))
carrier_means[max_mean_non_carrier_density_index]
non_carrier_means[max_mean_non_carrier_density_index]
