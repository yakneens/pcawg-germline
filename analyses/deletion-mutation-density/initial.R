library(VariantAnnotation)
library(data.table)
library(pcawg.common)
library(GenomicFeatures)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(rtracklayer)
library(org.Hs.eg.db)
library(COSMIC.67)
library(calibrate)
library(biomaRt)
library(ggplot2)
library(ggrepel)
library(ggbio)
library(BSgenome.Hsapiens.1000genomes.hs37d5)
library(data.table)
library(icd)

#Read in metadata
sample_meta = get_pcawg_metadata("~/Downloads/pcawg_data/sample_metadata/pcawg_summary.tsv", 
                                 "~/Downloads/pcawg_data/sample_metadata/PCAWG Blacklisted Donors - Blacklist_donors_2016_06_28.tsv", split_multi_tumors = F)

sample_meta = keep_first_of_multi_tumors(sample_meta)

clinical_meta = get_clinical_metadata("~/Downloads/pcawg_data/clinical_metadata/pcawg_donor_clinical_May2016_v2.tsv", sample_meta)

hist_meta = get_histology_metadata("~/Downloads/pcawg_data/clinical_metadata/pcawg_specimen_histology_May2016_v2.tsv", sample_meta)

sample_column_list = c("donor_unique_id", "normal_wgs_aliquot_id", 
                       "tumor_wgs_aliquot_id", "dcc_project_code")
clinical_column_list = c("donor_sex", "donor_vital_status", "donor_diagnosis_icd10",
                         "donor_age_at_diagnosis", "donor_survival_time",
                         "first_therapy_type", "first_therapy_response",
                         "donor_interval_of_last_followup", "tobacco_smoking_history_indicator",
                         "tobacco_smoking_intensity", "alcohol_history")
hist_column_list = c("organ_system", "histology_abbreviation", 
                     "histology_tier1", "histology_tier2",
                     "histology_tier3", "histology_tier4",
                     "tumour_histological_code", "tumour_histological_type",
                     "tumour_stage", "tumour_grade",
                     "percentage_cellularity", "level_of_cellularity")
donor_meta = cbind(sample_meta[,sample_column_list], 
                   clinical_meta[,clinical_column_list], 
                   hist_meta[,hist_column_list])

rm(sample_meta, sample_column_list, clinical_meta, clinical_column_list, hist_meta, hist_column_list)

#Load germline deletions data
load("~/Downloads/pcawg_data/germline_deletions/dels.Rdata")
#load("~/Downloads/pcawg_data/germline_deletions/dels_chr22.Rdata")

set_deletion_range_ends(deletions)

deletion_genotypes = geno(deletions)$GT[,match(donor_meta$normal_wgs_aliquot_id, colnames(geno(deletions)$GT))] 


#Load somatic SNV data
load("~/Downloads/pcawg_data/snv_samples.Rdata")
#load("~/Downloads/pcawg_data/snv_samples_chr22.Rdata")

snv_samples = snv_samples[match(donor_meta$tumor_wgs_aliquot_id, names(snv_samples))]

#Detect and remove samples that are missing deletions or SNV calls
missing_sample_indices = c(which(is.na(colnames(deletion_genotypes))), 
                           which(is.na(names(snv_samples))))

donor_meta = donor_meta[-missing_sample_indices,]
deletion_genotypes = deletion_genotypes[,-missing_sample_indices]
snv_samples = snv_samples[-missing_sample_indices]

colnames(deletion_genotypes) = donor_meta$donor_unique_id
names(snv_samples) = donor_meta$donor_unique_id

#Compute matrix with 1 for heterozygous carriers of a deletion, 0 for non-carriers (hom ref), and NA for others
deletion_carrier_mask = get_het_carrier_mask(deletion_genotypes)

#Number of het carriers per deletion
carrier_counts = rowSums(deletion_carrier_mask, na.rm = T)
donor_deletion_counts = colSums(deletion_carrier_mask, na.rm = T)

carrier_summary = as.data.table(table(carrier_counts))
carrier_summary[,carrier_counts := as.numeric(carrier_counts)]

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

pre_flanks = flank(deletion_ranges, width(deletion_ranges) / 2, start = T)
pre_hits = do.call(rbind, lapply(snv_samples, function(x) as.table(findOverlaps(pre_flanks, rowRanges(x)))))

post_flanks = flank(deletion_ranges, width(deletion_ranges) / 2, start = F)
post_hits = do.call(rbind, lapply(snv_samples, function(x) as.table(findOverlaps(post_flanks, rowRanges(x)))))

flank_hits = pre_hits + post_hits
normalized_flank_hits = flank_hits / snv_counts / del_widths[col(flank_hits)]

rm(pre_flanks, pre_hits, post_flanks, post_hits)


wilcox_pvals = c()
mean_counts = list()
flank_wilcox_pvals = c()
flank_ttest_pvals = c()
flank_cors = c()

num_hypotheses_dels = dim(normalized_snv_hits)[2]
num_hypotheses_donors = dim(normalized_snv_hits)[1]
significance_threshold_dels = 0.05 / num_hypotheses_dels
significance_threshold_donors = 0.05 / num_hypotheses_donors

# Test each deletion individually, among carriers and non-carriers, and among deletions and flanks 
for(i in 1 : dim(normalized_snv_hits)[2]) {
  mean_counts[[i]] = c(mean(normalized_snv_hits[deletion_carrier_mask[,i] == 1,i], na.rm = T), 
                       mean(normalized_snv_hits[deletion_carrier_mask[,i] == 0,i], na.rm = T),
                       mean(normalized_flank_hits[deletion_carrier_mask[,i] == 1,i], na.rm = T))
  wilcox_pvals[i] = wilcox.test(normalized_snv_hits[deletion_carrier_mask[,i] == 1,i], 
                                normalized_snv_hits[deletion_carrier_mask[,i] == 0,i])$p.value
  flank_wilcox_pvals[i] = wilcox.test(na.omit(normalized_snv_hits[deletion_carrier_mask[,i] == 1,i]), 
                                      na.omit(normalized_flank_hits[deletion_carrier_mask[,i] == 1,i]), paired = T)$p.value
  flank_ttest_pvals[i] = t.test(na.omit(normalized_snv_hits[deletion_carrier_mask[,i] == 1,i]), 
                                na.omit(normalized_flank_hits[deletion_carrier_mask[,i] == 1,i]), paired = T)$p.value
  flank_cors[i] = cor(na.omit(normalized_snv_hits[deletion_carrier_mask[,i] == 1,i]), 
                      na.omit(normalized_flank_hits[deletion_carrier_mask[,i] == 1,i]))
}

results_by_del = as.data.table(do.call(rbind, mean_counts))
setnames(results_by_del, c("carriers", "non_carriers", "flanks"))

results_by_del[, `:=` (c_vs_nc_wilcox_pvals = wilcox_pvals,
                       c_vs_nc_mean_difs = carriers - non_carriers,
                       c_vs_nc_mean_ratios = carriers / non_carriers,
                       del_widths = width(deletion_ranges),
                       carrier_counts = carrier_counts[-deletion_filter],
                       chr = as.character(seqnames(deletion_ranges)),
                       d_vs_f_wilcox_pvals = unlist(flank_wilcox_pvals),
                       d_vs_f_ttest_pvals = unlist(flank_ttest_pvals),
                       d_vs_f_cors = unlist(flank_cors),
                       d_vs_f_mean_difs = carriers - flanks,
                       d_vs_f_mean_ratios = carriers / flanks
                       )]


donor_counts = list()
donor_wilcox_pvals = list()
donor_ttest_pvals = list()
donor_cors = list()
for(i in 1 : dim(normalized_snv_hits)[1]){
  donor_counts[[i]] = c(mean(normalized_snv_hits[i, deletion_carrier_mask[i,] == 1], na.rm = T), 
                          mean(normalized_flank_hits[i, deletion_carrier_mask[i,] == 1], na.rm = T))
  donor_wilcox_pvals[[i]] = wilcox.test(na.omit(normalized_snv_hits[i, deletion_carrier_mask[i,] == 1]), 
                                      na.omit(normalized_flank_hits[i, deletion_carrier_mask[i,] == 1]), paired = T)$p.value
  donor_ttest_pvals[[i]] = t.test(na.omit(normalized_snv_hits[i, deletion_carrier_mask[i,] == 1]), 
                                na.omit(normalized_flank_hits[i, deletion_carrier_mask[i,] == 1]), paired = T)$p.value
  donor_cors[[i]] = cor(na.omit(normalized_snv_hits[i, deletion_carrier_mask[i,] == 1]), 
                      na.omit(normalized_flank_hits[i, deletion_carrier_mask[i,] == 1]))
}

results_by_donor = as.data.table(do.call(rbind, donor_counts))
colnames(results_by_donor) = c("dels", "flanks")
results_by_donor = cbind(results_by_donor, donor_meta)
results_by_donor[, `:=` (wilcox_pvals = unlist(donor_wilcox_pvals),
                         ttest_pvals = unlist(donor_ttest_pvals),
                         cors = unlist(donor_cors),
                         mean_ratios = dels / flanks,
                         snv_counts = snv_counts,
                         donor_deletion_counts = donor_deletion_counts,
                         age_bin = findInterval(donor_age_at_diagnosis, seq(0,100,10)))]

rm(donor_counts, donor_wilcox_pvals, donor_ttest_pvals, donor_cors)

filtered_results_by_donor = results_by_donor[which(is.finite(mean_ratios) & donor_diagnosis_icd10 != "")]
donor_count_model = glm(data=filtered_results_by_donor, formula=dels ~ donor_diagnosis_icd10 + tumour_stage + tumour_histological_code + flanks + level_of_cellularity + tumour_grade + dcc_project_code + snv_counts)
anova(donor_count_model, test="Chisq")
