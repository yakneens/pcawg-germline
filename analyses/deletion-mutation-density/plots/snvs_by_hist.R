#SNVs by ICD10 Disease Code
snv_counts_by_tumor_type = data.table(log10(snv_counts), donor_meta$donor_diagnosis_icd10)
setnames(snv_counts_by_tumor_type, c("donor_snv_counts", "tumor_histology"))

to_remap = cbind(c("25", "25.1", "25.0", "25.2", "25.0/25.9", "25.7", "25.1/25.2", "c92.0", "C92.00", "C26.8"), c("C25", "C25.1", "C25.0", "C25.2", "C25.0", "C25.7", "C25.1", "C92.0", "C92.0", "C26.9"))
apply(to_remap, 1, function(x) snv_counts_by_tumor_type[tumor_histology == x[1], tumor_histology := x[2]])


filtered_stt = snv_counts_by_tumor_type[, `:=`(count = .N), by=tumor_histology][count > 2 & tumor_histology != "" & tumor_histology != "Unknown" & tumor_histology != "Unknown (Periampullary?)"]

ggplot(filtered_stt, aes(x=reorder(tumor_histology, donor_snv_counts, FUN=median), y=donor_snv_counts)) + geom_boxplot() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  xlab("Tumor Histology ICD10 Code") + ylab(expression(paste(log[10]("# SNVs"))))

icd_lookup = cbind(unique(filtered_stt[,.(tumor_histology, count)]), icd_explain(filtered_stt$tumor_histology))
colnames(icd_lookup) = c("ICD10 Code", "# Patients", "ICD10 Description")
setcolorder(icd_lookup, c(1,3,2))
setorder(icd_lookup, -"# Patients")

#SNVs by Histology Tier 4
snv_counts_by_tumor_type = data.table(log10(snv_counts), donor_meta$histology_tier4)
setnames(snv_counts_by_tumor_type, c("donor_snv_counts", "tumor_histology"))
filtered_stt = snv_counts_by_tumor_type[, `:=`(count = .N), by=tumor_histology][count > 2 & tumor_histology != "" & tumor_histology != "Unknown" & tumor_histology != "Unknown (Periampullary?)"]

ggplot(filtered_stt, aes(x=reorder(tumor_histology, donor_snv_counts, FUN=median), y=donor_snv_counts)) + geom_boxplot() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  xlab("Tumor Histology") + ylab(expression(paste(log[10]("# SNVs"))))

#SNVs by Histology Tier 3
snv_counts_by_tumor_type = data.table(log10(snv_counts), donor_meta$histology_tier3)
setnames(snv_counts_by_tumor_type, c("donor_snv_counts", "tumor_histology"))
filtered_stt = snv_counts_by_tumor_type[, `:=`(count = .N), by=tumor_histology][count > 2 & tumor_histology != "" & tumor_histology != "Unknown" & tumor_histology != "Unknown (Periampullary?)"]

ggplot(filtered_stt, aes(x=reorder(tumor_histology, donor_snv_counts, FUN=median), y=donor_snv_counts)) + geom_boxplot() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  xlab("Tumor Histology") + ylab(expression(paste(log[10]("# SNVs"))))

#SNVs by Histology Tier 2
snv_counts_by_tumor_type = data.table(log10(snv_counts), donor_meta$histology_tier2)
setnames(snv_counts_by_tumor_type, c("donor_snv_counts", "tumor_histology"))
filtered_stt = snv_counts_by_tumor_type[, `:=`(count = .N), by=tumor_histology][count > 2 & tumor_histology != "" & tumor_histology != "Unknown" & tumor_histology != "Unknown (Periampullary?)"]

ggplot(filtered_stt, aes(x=reorder(tumor_histology, donor_snv_counts, FUN=median), y=donor_snv_counts)) + geom_boxplot() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  xlab("Tumor Histology") + ylab(expression(paste(log[10]("# SNVs"))))

#SNVs by Histology Tier 1
snv_counts_by_tumor_type = data.table(log10(snv_counts), donor_meta$histology_tier1)
setnames(snv_counts_by_tumor_type, c("donor_snv_counts", "tumor_histology"))
filtered_stt = snv_counts_by_tumor_type[, `:=`(count = .N), by=tumor_histology][count > 2 & tumor_histology != "" & tumor_histology != "Unknown" & tumor_histology != "Unknown (Periampullary?)"]

ggplot(filtered_stt, aes(x=reorder(tumor_histology, donor_snv_counts, FUN=median), y=donor_snv_counts)) + geom_boxplot() +
  xlab("Tumor Histology") + ylab(expression(paste(log[10]("# SNVs"))))

#SNVs by Histology Abbreviation
snv_counts_by_tumor_type = data.table(log10(snv_counts), donor_meta$histology_abbreviation)
setnames(snv_counts_by_tumor_type, c("donor_snv_counts", "tumor_histology"))
filtered_stt = snv_counts_by_tumor_type[, `:=`(count = .N), by=tumor_histology][count > 2 & tumor_histology != "" & tumor_histology != "Unknown" & tumor_histology != "Unknown (Periampullary?)"]

ggplot(filtered_stt, aes(x=reorder(tumor_histology, donor_snv_counts, FUN=median), y=donor_snv_counts)) + geom_boxplot() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  xlab("Tumor Histology") + ylab(expression(paste(log[10]("# SNVs"))))

#SNVs by Organ System
snv_counts_by_tumor_type = data.table(log10(snv_counts), donor_meta$organ_system)
setnames(snv_counts_by_tumor_type, c("donor_snv_counts", "tumor_histology"))
filtered_stt = snv_counts_by_tumor_type[, `:=`(count = .N), by=tumor_histology][count > 2 & tumor_histology != "" & tumor_histology != "Unknown" & tumor_histology != "Unknown (Periampullary?)"]

ggplot(filtered_stt, aes(x=reorder(tumor_histology, donor_snv_counts, FUN=median), y=donor_snv_counts)) + geom_boxplot() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  xlab("Organ System") + ylab(expression(paste(log[10]("# SNVs"))))

#SNVs by Age
snv_counts_by_age = data.table(log10(snv_counts), donor_meta$donor_age_at_diagnosis)
setnames(snv_counts_by_age, c("donor_snv_counts", "donor_age"))
snv_counts_by_age[, age_bin := findInterval(donor_age, seq(0,100,10))]
snv_counts_by_age = snv_counts_by_age[!is.na(donor_age)]
ggplot(snv_counts_by_age, aes(x=as.factor(age_bin), y=donor_snv_counts)) + geom_boxplot() +
  xlab("Donor Age at Diagnosis") + ylab(expression(paste(log[10]("# SNVs")))) + scale_x_discrete(labels=paste(seq(0,90,10), "-", seq(10,100,10)))


