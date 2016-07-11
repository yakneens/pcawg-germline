exclusions = c("", "Unknown", "Unknown (Periampullary?)")

plot_ordered_boxes(results_by_donor, "donor_diagnosis_icd10", "snv_counts", median, 2, exclusions, "Tumor Histology ICD10 Code", expression(paste(log[10]("# SNVs"))))

plot_ordered_boxes(results_by_donor, "histology_tier4", "snv_counts", median, 2, exclusions, "Tumor Histology", expression(paste(log[10]("# SNVs"))))

plot_ordered_boxes(results_by_donor, "histology_tier3", "snv_counts", median, 2, exclusions, "Tumor Histology", expression(paste(log[10]("# SNVs"))))

plot_ordered_boxes(results_by_donor, "histology_tier2", "snv_counts", median, 2, exclusions, "Tumor Histology", expression(paste(log[10]("# SNVs"))))

plot_ordered_boxes(results_by_donor, "histology_tier1", "snv_counts", median, 2, exclusions, "Tumor Histology", expression(paste(log[10]("# SNVs"))))

plot_ordered_boxes(results_by_donor, "histology_abbreviation", "snv_counts", median, 2, exclusions, "Tumor Histology", expression(paste(log[10]("# SNVs"))))

plot_ordered_boxes(results_by_donor, "organ_system", "snv_counts", median, 2, exclusions, "Organ System", expression(paste(log[10]("# SNVs"))))

plot_ordered_boxes(results_by_donor, "dcc_project_code", "snv_counts", median, 2, exclusions, "Organ System", expression(paste(log[10]("# SNVs"))))


#SNVs by Age
#snv_counts_by_age = snv_counts_by_age[!is.na(donor_age)]
ggplot(results_by_donor, aes(x=as.factor(age_bin), y=log10(snv_counts))) + geom_boxplot() +
  xlab("Donor Age at Diagnosis") + ylab(expression(paste(log[10]("# SNVs")))) + scale_x_discrete(labels=paste(seq(0,90,10), "-", seq(10,100,10)))


icd_lookup = cbind(unique(filtered_stt[,.(histology_tier4, count)]), icd_explain(filtered_stt$histology_tier4))
colnames(icd_lookup) = c("ICD10 Code", "# Patients", "ICD10 Description")
setcolorder(icd_lookup, c(1,3,2))
setorder(icd_lookup, -"# Patients")
