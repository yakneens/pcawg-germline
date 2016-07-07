filtered_dtt = results_by_donor[, `:=`(count = .N), by=donor_diagnosis_icd10][count > 2 & donor_diagnosis_icd10 != "" & donor_diagnosis_icd10 != "Unknown" & donor_diagnosis_icd10 != "Unknown (Periampullary?)"]

ggplot(filtered_dtt, aes(x=reorder(donor_diagnosis_icd10, donor_deletion_counts, FUN=median), y=donor_deletion_counts)) + geom_boxplot() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  xlab("Tumor Histology ICD10 Code") + ylab("# Deletions")

icd_lookup = cbind(unique(filtered_dtt[,.(tumor_histology, count)]), icd_explain(filtered_dtt$donor_diagnosis_icd10))
colnames(icd_lookup) = c("ICD10 Code", "# Patients", "ICD10 Description")
setcolorder(icd_lookup, c(1,3,2))
setorder(icd_lookup, -"# Patients")

#Deletions by Histology Tier 4
filtered_dtt = results_by_donor[, `:=`(count = .N), by=histology_tier4][count > 2 & histology_tier4 != "" & histology_tier4 != "Unknown" & histology_tier4 != "Unknown (Periampullary?)"]

ggplot(filtered_dtt, aes(x=reorder(histology_tier4, donor_deletion_counts, FUN=median), y=donor_deletion_counts)) + geom_boxplot() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  xlab("Tumor Histology") + ylab("# Deletions")

#Deletions by Histology Tier 3
filtered_dtt = results_by_donor[, `:=`(count = .N), by=histology_tier3][count > 2 & histology_tier3 != "" & histology_tier3 != "Unknown" & histology_tier3 != "Unknown (Periampullary?)"]

ggplot(filtered_dtt, aes(x=reorder(histology_tier3, donor_deletion_counts, FUN=median), y=donor_deletion_counts)) + geom_boxplot() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  xlab("Tumor Histology") + ylab("# Deletions")

#Deletions by Histology Tier 2
filtered_dtt = results_by_donor[, `:=`(count = .N), by=histology_tier2][count > 2 & histology_tier2 != "" & histology_tier2 != "Unknown" & histology_tier2 != "Unknown (Periampullary?)"]

ggplot(filtered_dtt, aes(x=reorder(histology_tier2, donor_deletion_counts, FUN=median), y=donor_deletion_counts)) + geom_boxplot() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  xlab("Tumor Histology") + ylab("# Deletions")

#Deletions by Histology Tier 1
filtered_dtt = results_by_donor[, `:=`(count = .N), by=histology_tier1][count > 2 & histology_tier1 != "" & histology_tier1 != "Unknown" & histology_tier1 != "Unknown (Periampullary?)"]

ggplot(filtered_dtt, aes(x=reorder(histology_tier1, donor_deletion_counts, FUN=median), y=donor_deletion_counts)) + geom_boxplot() +
  xlab("Tumor Histology") + ylab("# Deletions")

#Deletions by Histology Abbreviation
filtered_dtt = results_by_donor[, `:=`(count = .N), by=histology_abbreviation][count > 2 & histology_abbreviation != "" & histology_abbreviation != "Unknown" & histology_abbreviation != "Unknown (Periampullary?)"]

ggplot(filtered_dtt, aes(x=reorder(histology_abbreviation, donor_deletion_counts, FUN=median), y=donor_deletion_counts)) + geom_boxplot() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  xlab("Tumor Histology") + ylab("# Deletions")

#Deletions by Organ System
filtered_dtt = results_by_donor[, `:=`(count = .N), by=organ_system][count > 2 & organ_system != "" & organ_system != "Unknown" & organ_system != "Unknown (Periampullary?)"]

ggplot(filtered_dtt, aes(x=reorder(organ_system, donor_deletion_counts, FUN=median), y=donor_deletion_counts)) + geom_boxplot() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  xlab("Organ System") + ylab("# Deletions")

#Deletions by Age
#deletions_by_age = deletions_by_age[!is.na(donor_age)]
ggplot(results_by_donor, aes(x=as.factor(age_bin), y=donor_deletion_counts)) + geom_boxplot() +
  xlab("Donor Age at Diagnosis") + ylab("# Deletions") + scale_x_discrete(labels=paste(seq(0,90,10), "-", seq(10,100,10)))

#Deletions by project code
ggplot(results_by_donor, aes(x=reorder(dcc_project_code, donor_deletion_counts, FUN=median), y=donor_deletion_counts)) + 
  geom_boxplot() + xlab("Project Code") + ylab("# Deletions") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) 