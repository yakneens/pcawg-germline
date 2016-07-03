#Deletions by ICD10 Disease Code
deletions_by_tumor_type = data.table(donor_deletion_counts, donor_meta$donor_diagnosis_icd10)
setnames(deletions_by_tumor_type, c("donor_deletion_counts", "tumor_histology"))

to_remap = cbind(c("25", "25.1", "25.0", "25.2", "25.0/25.9", "25.7", "25.1/25.2", "c92.0", "C92.00", "C26.8"), c("C25", "C25.1", "C25.0", "C25.2", "C25.0", "C25.7", "C25.1", "C92.0", "C92.0", "C26.9"))
apply(to_remap, 1, function(x) deletions_by_tumor_type[tumor_histology == x[1], tumor_histology := x[2]])


filtered_dtt = deletions_by_tumor_type[, `:=`(count = .N), by=tumor_histology][count > 2 & tumor_histology != "" & tumor_histology != "Unknown" & tumor_histology != "Unknown (Periampullary?)"]

ggplot(filtered_dtt, aes(x=reorder(tumor_histology, donor_deletion_counts, FUN=median), y=donor_deletion_counts)) + geom_boxplot() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  xlab("Tumor Histology ICD10 Code") + ylab("# Deletions")

icd_lookup = cbind(unique(filtered_dtt[,.(tumor_histology, count)]), icd_explain(filtered_dtt$tumor_histology))
colnames(icd_lookup) = c("ICD10 Code", "# Patients", "ICD10 Description")
setcolorder(icd_lookup, c(1,3,2))
setorder(icd_lookup, -"# Patients")

#Deletions by Histology Tier 4
deletions_by_tumor_type = data.table(donor_deletion_counts, donor_meta$histology_tier4)
setnames(deletions_by_tumor_type, c("donor_deletion_counts", "tumor_histology"))
filtered_dtt = deletions_by_tumor_type[, `:=`(count = .N), by=tumor_histology][count > 2 & tumor_histology != "" & tumor_histology != "Unknown" & tumor_histology != "Unknown (Periampullary?)"]

ggplot(filtered_dtt, aes(x=reorder(tumor_histology, donor_deletion_counts, FUN=median), y=donor_deletion_counts)) + geom_boxplot() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  xlab("Tumor Histology") + ylab("# Deletions")

#Deletions by Histology Tier 3
deletions_by_tumor_type = data.table(donor_deletion_counts, donor_meta$histology_tier3)
setnames(deletions_by_tumor_type, c("donor_deletion_counts", "tumor_histology"))
filtered_dtt = deletions_by_tumor_type[, `:=`(count = .N), by=tumor_histology][count > 2 & tumor_histology != "" & tumor_histology != "Unknown" & tumor_histology != "Unknown (Periampullary?)"]

ggplot(filtered_dtt, aes(x=reorder(tumor_histology, donor_deletion_counts, FUN=median), y=donor_deletion_counts)) + geom_boxplot() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  xlab("Tumor Histology") + ylab("# Deletions")

#Deletions by Histology Tier 2
deletions_by_tumor_type = data.table(donor_deletion_counts, donor_meta$histology_tier2)
setnames(deletions_by_tumor_type, c("donor_deletion_counts", "tumor_histology"))
filtered_dtt = deletions_by_tumor_type[, `:=`(count = .N), by=tumor_histology][count > 2 & tumor_histology != "" & tumor_histology != "Unknown" & tumor_histology != "Unknown (Periampullary?)"]

ggplot(filtered_dtt, aes(x=reorder(tumor_histology, donor_deletion_counts, FUN=median), y=donor_deletion_counts)) + geom_boxplot() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  xlab("Tumor Histology") + ylab("# Deletions")

#Deletions by Histology Tier 1
deletions_by_tumor_type = data.table(donor_deletion_counts, donor_meta$histology_tier1)
setnames(deletions_by_tumor_type, c("donor_deletion_counts", "tumor_histology"))
filtered_dtt = deletions_by_tumor_type[, `:=`(count = .N), by=tumor_histology][count > 2 & tumor_histology != "" & tumor_histology != "Unknown" & tumor_histology != "Unknown (Periampullary?)"]

ggplot(filtered_dtt, aes(x=reorder(tumor_histology, donor_deletion_counts, FUN=median), y=donor_deletion_counts)) + geom_boxplot() +
  xlab("Tumor Histology") + ylab("# Deletions")

#Deletions by Histology Abbreviation
deletions_by_tumor_type = data.table(donor_deletion_counts, donor_meta$histology_abbreviation)
setnames(deletions_by_tumor_type, c("donor_deletion_counts", "tumor_histology"))
filtered_dtt = deletions_by_tumor_type[, `:=`(count = .N), by=tumor_histology][count > 2 & tumor_histology != "" & tumor_histology != "Unknown" & tumor_histology != "Unknown (Periampullary?)"]

ggplot(filtered_dtt, aes(x=reorder(tumor_histology, donor_deletion_counts, FUN=median), y=donor_deletion_counts)) + geom_boxplot() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  xlab("Tumor Histology") + ylab("# Deletions")

#Deletions by Organ System
deletions_by_tumor_type = data.table(donor_deletion_counts, donor_meta$organ_system)
setnames(deletions_by_tumor_type, c("donor_deletion_counts", "tumor_histology"))
filtered_dtt = deletions_by_tumor_type[, `:=`(count = .N), by=tumor_histology][count > 2 & tumor_histology != "" & tumor_histology != "Unknown" & tumor_histology != "Unknown (Periampullary?)"]

ggplot(filtered_dtt, aes(x=reorder(tumor_histology, donor_deletion_counts, FUN=median), y=donor_deletion_counts)) + geom_boxplot() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  xlab("Organ System") + ylab("# Deletions")

#Deletions by Age
deletions_by_age = data.table(donor_deletion_counts, donor_meta$donor_age_at_diagnosis)
setnames(deletions_by_age, c("donor_deletion_counts", "donor_age"))
deletions_by_age[, age_bin := findInterval(donor_age, seq(0,100,10))]
deletions_by_age = deletions_by_age[!is.na(donor_age)]
ggplot(deletions_by_age, aes(x=as.factor(age_bin), y=donor_deletion_counts)) + geom_boxplot() +
  xlab("Donor Age at Diagnosis") + ylab("# Deletions") + scale_x_discrete(labels=paste(seq(0,90,10), "-", seq(10,100,10)))


