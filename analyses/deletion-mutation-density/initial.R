library(VariantAnnotation)
library(data.table)


#Read in sample metadata
sample_meta = read.table("~/Downloads/pcawg_data/sample_metadata/pcawg_summary.tsv", header=TRUE, sep="\t")

#Read in deletions from Chromosome 22 only 
rng = GRanges(seqnames="22", ranges=IRanges(start=0, end=51304566))
tab = TabixFile("~/Downloads/pcawg_data/germline_deletions/DEL.pcawg.rf.vcf.gz")
dels = readVcf(tab, "hs37d5", param=rng)

#When deletions are imported they don't seem to have a proper end recorded
end(ranges(rowRanges(dels))) = info(dels)$END

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

carrier_hits = calculate_hits(het_carriers)
carrier_means = lapply(carrier_hits, mean) 
#normalized_carrier_means = unlist(carrier_means) 

non_carrier_hits = calculate_hits(non_carriers)
non_carrier_means = lapply(non_carrier_hits, mean) 
#normalized_non_carrier_means = unlist(non_carrier_means)


