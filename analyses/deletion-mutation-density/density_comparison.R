library(VariantAnnotation)
library(data.table)
library(devtools)
library(GenomicFeatures)
devtools::load_all(pkg="~/Documents/workspace/pcawg-germline/analyses/deletion-mutation-density/pcawg.common")
library(pcawg.common)

library(ggplot2)

source("~/Documents/workspace/pcawg-germline/analyses/deletion-mutation-density/get_binned_densities.R")

sample_meta = get_pcawg_metadata("~/Downloads/pcawg_data/sample_metadata/pcawg_summary.tsv", 
                                 "~/Downloads/pcawg_data/sample_metadata/PCAWG Excluded Donors%2FSamples - Excluded_donors_2016_08_30.tsv", split_multi_tumors = F)

sample_meta = keep_first_of_multi_tumors(sample_meta)

clinical_meta = get_clinical_metadata("~/Downloads/pcawg_data/clinical_metadata/pcawg_donor_clinical_August2016_v6.tsv", sample_meta)

hist_meta = get_histology_metadata("~/Downloads/pcawg_data/clinical_metadata/pcawg_specimen_histology_August2016_v6.tsv", sample_meta)

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
donor_meta = data.table(cbind(sample_meta[,sample_column_list], 
                              clinical_meta[,clinical_column_list], 
                              hist_meta[,hist_column_list]))

rm(sample_meta, sample_column_list, clinical_meta, clinical_column_list, hist_meta, hist_column_list)

#Load germline deletions data
load("~/Downloads/pcawg_data/dels_28_11_2016.Rdata")
#load("~/Downloads/pcawg_data/germline_deletions/dels_chr22.Rdata")
load("~/Downloads/pcawg_data/del_density/input_data/deletion_ranges.RData")
set_deletion_range_ends(deletions)

deletion_genotypes = geno(deletions)$GT[,match(donor_meta$normal_wgs_aliquot_id, colnames(geno(deletions)$GT))] 
deletion_ranges = rowRanges(deletions)
del_widths = width(deletion_ranges)


#Load somatic SNV data
load("~/Downloads/pcawg_data/snv_samples.Rdata")
#load("~/Downloads/pcawg_data/snv_sample_ranges.Rdata")
#snv_samples = snv_ranges
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

snv_ranges = lapply(snv_samples, rowRanges)
snv_counts = unlist(lapply(snv_ranges, length))
rm(snv_samples)

#Filter samples by somatic mutation rate
snv_sample_filter = which(unlist(lapply(snv_ranges, length)) < 1000)
snv_ranges = snv_ranges[-snv_sample_filter]
deletion_genotypes = deletion_genotypes[,-snv_sample_filter]
donor_meta = donor_meta[-snv_sample_filter,]

#Compute matrix with 1 for heterozygous carriers of a deletion, 0 for non-carriers (hom ref), and NA for others
deletion_carrier_mask = get_het_carrier_mask(deletion_genotypes)





deletion_filter = NULL
snv_filter = NULL
#deletion_filter = which(del_widths < 10000 | del_widths > 100000)
#deletion_filter = which(del_widths < 1000 | del_widths > 10000)
#deletion_filter = which(del_widths > 1000)
deletion_filter = which(seqnames(deletion_ranges) != "8")
#deletion_filter = which(del_widths < 100000)


raw_densities = get_binned_densities(deletion_filter, NULL)
del_bins = data.table(seq(1,30), rowSums(raw_densities))
setnames(del_bins, c("index", "bins"))

max_density = max(del_bins$bins) 
ggplot(del_bins, aes(x=index, y=bins)) + geom_bar(stat="identity", fill="grey70") + geom_vline(xintercept=10.5, linetype="dashed") + geom_vline(xintercept=20.5, linetype="dashed") + xlab("Bin") + ylab("SNV Density") + annotate("text", x=5, y=1.1*max_density, label="Left Flank") + annotate("text", x=15, y=1.1*max_density, label="Deletion") + annotate("text", x=25, y=1.1*max_density, label="Right Flank")

del_bins2 = copy(del_bins)
del_bins2[index > 10 & index <= 20, bins := 2 * bins]
max_density = max(del_bins2$bins) 

ggplot(del_bins2, aes(x=index, y=bins)) + geom_bar(stat="identity", fill="grey70") + geom_vline(xintercept=10.5, linetype="dashed") + geom_vline(xintercept=20.5, linetype="dashed") + xlab("Bin") + ylab("SNV Density") + annotate("text", x=5, y=1.1*max_density, label="Left Flank") + annotate("text", x=15, y=1.1*max_density, label="Deletion") + annotate("text", x=25, y=1.1*max_density, label="Right Flank")


bins_by_cancer_project <- function(){
  dir.create(path="~/Downloads/pcawg_data/del_density/by_project")
  project_codes = unique(donor_meta[,dcc_project_code])
  deletion_filter = NULL
  
  sapply(project_codes, 
         function(x){print(x);get_binned_densities(NULL, which(donor_meta$dcc_project_code != x))},
         simplify = F,
         USE.NAMES = T)
  
}

project_based_bins = bins_by_cancer_project()

bins_by_snv_type <- function(){
  
}