library(VariantAnnotation)

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

snv_samples = read_snv_samples_from_vcf("~/Downloads/pcawg_data/test_snv/", sample_meta)
snv_genotypes = geno(unlist(snv_samples))$GT
snv_carrier_mask = get_het_carrier_mask(snv_samples)