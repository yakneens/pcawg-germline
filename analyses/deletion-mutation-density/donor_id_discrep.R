sample_m = read.table("~/Downloads/pcawg_data/sample_metadata/pcawg_summary.tsv", header=TRUE, sep="\t", stringsAsFactors = F)

hist_m = read.csv("~/Downloads/pcawg_data/clinical_metadata/pcawg_specimen_histology_May2016_v2.tsv", sep = "\t", stringsAsFactors = F)

blacklist_m = read.csv("~/Downloads/pcawg_data/sample_metadata/PCAWG Blacklisted Donors - Blacklist_donors_2016_06_28.tsv", sep="\t", stringsAsFactors = F)

sample_m_donors = sample_m[,c("donor_unique_id", "submitter_donor_id")]
hist_m_donors = hist_m[,c("donor_unique_id", "submitted_donor_id")]
colnames(hist_m_donors)[2] = "submitter_donor_id"

sample_m_donors = sample_m_donors[-na.omit(match(blacklist_m$donor_unique_id, sample_m_donors$donor_unique_id)),]
hist_m_donors = hist_m_donors[-na.omit(match(blacklist_m$donor_unique_id, hist_m_donors$donor_unique_id)),]

s_not_in_h = sample_m_donors[which(!(sample_m_donors$donor_unique_id %in% hist_m_donors$donor_unique_id)),]
h_not_in_s = sample_m_donors[which(!(hist_m_donors$donor_unique_id %in% sample_m_donors$donor_unique_id)),]

