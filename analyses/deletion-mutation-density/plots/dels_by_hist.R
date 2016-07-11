exclusions = c("", "Unknown", "Unknown (Periampullary?)")

plot_ordered_boxes(results_by_donor, "donor_diagnosis_icd10", "donor_deletion_counts", median, 2, exclusions, "Tumor Histology ICD10 Code", "# Deletions")

plot_ordered_boxes(results_by_donor, "histology_tier4", "donor_deletion_counts", median, 2, exclusions, "Tumor Histology", "# Deletions")

plot_ordered_boxes(results_by_donor, "histology_tier3", "donor_deletion_counts", median, 2, exclusions, "Tumor Histology", "# Deletions")

plot_ordered_boxes(results_by_donor, "histology_tier2", "donor_deletion_counts", median, 2, exclusions, "Tumor Histology", "# Deletions")

plot_ordered_boxes(results_by_donor, "histology_tier1", "donor_deletion_counts", median, 2, exclusions, "Tumor Histology", "# Deletions")

plot_ordered_boxes(results_by_donor, "histology_abbreviation", "donor_deletion_counts", median, 2, exclusions, "Tumor Histology", "# Deletions")

plot_ordered_boxes(results_by_donor, "organ_system", "donor_deletion_counts", median, 2, exclusions, "Organ System", "# Deletions")

plot_ordered_boxes(results_by_donor, "dcc_project_code", "donor_deletion_counts", median, 2, exclusions, "Organ System", "# Deletions")

icd_lookup = cbind(unique(filtered_dtt[,.(tumor_histology, count)]), icd_explain(filtered_dtt$donor_diagnosis_icd10))
colnames(icd_lookup) = c("ICD10 Code", "# Patients", "ICD10 Description")
setcolorder(icd_lookup, c(1,3,2))
setorder(icd_lookup, -"# Patients")

#Deletions by Age
#deletions_by_age = deletions_by_age[!is.na(donor_age)]
ggplot(results_by_donor, aes(x=as.factor(age_bin), y=donor_deletion_counts)) + geom_boxplot() +
  xlab("Donor Age at Diagnosis") + ylab("# Deletions") + scale_x_discrete(labels=paste(seq(0,90,10), "-", seq(10,100,10)))
