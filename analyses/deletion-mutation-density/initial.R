library(VariantAnnotation)
library(data.table)
library(pcawg.common)

#Read in metadata
sample_meta = get_pcawg_metadata("~/Downloads/pcawg_data/sample_metadata/pcawg_summary.tsv", split_multi_tumors = F)

sample_meta = keep_first_of_multi_tumors(sample_meta)

#Load germline deletions data
#load("~/Downloads/pcawg_data/germline_deletions/dels.Rdata")
load("~/Downloads/pcawg_data/germline_deletions/dels_chr22.Rdata")

set_deletion_range_ends(deletions)

deletion_genotypes = geno(deletions)$GT[,match(sample_meta$normal_wgs_aliquot_id, colnames(geno(deletions)$GT))] 


#Load somatic SNV data
#load("~/Downloads/pcawg_data/snv_samples.Rdata")
load("~/Downloads/pcawg_data/snv_samples_chr22.Rdata")

snv_samples = snv_samples[match(sample_meta$tumor_wgs_aliquot_id, names(snv_samples))]

#Detect and remove samples that are missing deletions or SNV calls
missing_sample_indices = c(which(is.na(colnames(deletion_genotypes))), which(is.na(names(snv_samples))))

sample_meta = sample_meta[-missing_sample_indices,]
deletion_genotypes = deletion_genotypes[,-missing_sample_indices]
snv_samples = snv_samples[-missing_sample_indices]

colnames(deletion_genotypes) = sample_meta$donor_unique_id
names(snv_samples) = sample_meta$donor_unique_id

#Compute matrix with 1 for heterozygous carriers of a deletion, 0 for non-carriers (hom ref), and NA for others
deletion_carrier_mask = get_het_carrier_mask(deletion_genotypes)

#Number of het carriers per deletion
carrier_counts = rowSums(deletion_carrier_mask, na.rm = T)

min_carrier_count = 5

deletion_filter = which(carrier_counts < min_carrier_count)

deletion_ranges = rowRanges(deletions)[-deletion_filter]
deletion_genotypes = deletion_genotypes[-deletion_filter,]
deletion_carrier_mask = t(deletion_carrier_mask[-deletion_filter,])

#Compute a matrix with count of SNVs overlying each deletion
snv_hits = do.call(rbind, lapply(snv_samples, function(x) as.table(findOverlaps(deletion_ranges, rowRanges(x)))))

colnames(snv_hits) = rownames(deletion_genotypes)

# Normalize each sample by number of SNVs in that sample
# Normalize each deletion by deletion width
del_widths =  width(ranges(deletion_ranges))
snv_counts = unlist(lapply(snv_samples, length))
  
normalized_snv_hits = snv_hits / snv_counts / del_widths[col(snv_hits)]

dim(normalized_snv_hits)

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
