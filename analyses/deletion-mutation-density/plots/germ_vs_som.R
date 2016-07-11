plot_ranks <- function(data, col_1, col_2, comp_col, exclusions, count_cutoff, order_func){
  filtered_data = data[, count := .N, by=get(comp_col)][count > count_cutoff & !get(comp_col) %in% exclusions]
  
  col_1_ranked = filtered_data[,med_count := order_func(as.double(get(col_1))), by=get(comp_col)][order(-med_count), unique(get(comp_col))]
  
  col_2_ranked = filtered_data[,med_count := order_func(as.double(get(col_2))), by=comp_col][order(-med_count), unique(get(comp_col))]
  
  ranks = as.data.table(list(seq(1, length(col_1_ranked)), match(col_1_ranked, col_2_ranked), col_1_ranked))
  colnames(ranks) = c("col_1", "col_2", "labels")
  print(
    ggplot(data=ranks) + geom_point(aes(x=col_1, y=col_2)) + 
      xlab(col_1) + ylab(col_2) +
      ggtitle(label=paste("Rank of", comp_col, "by", order_func, col_1, "vs", col_2, 
                  ", r =", round(cor(ranks$col_1, 
                                                ranks$col_2), 2)))
  )
}

exclusions = c("", "Unknown", "Unknown (Periampullary?)")

plot_ranks(results_by_donor, "snv_counts", "donor_deletion_counts", "donor_diagnosis_icd10", exclusions, 2, median)
plot_ranks(results_by_donor, "snv_counts", "donor_deletion_counts", "histology_tier4", exclusions, 2, median)
plot_ranks(results_by_donor, "snv_counts", "donor_deletion_counts", "histology_tier3", exclusions, 2, median)
plot_ranks(results_by_donor, "snv_counts", "donor_deletion_counts", "histology_tier2", exclusions, 2, median)
plot_ranks(results_by_donor, "snv_counts", "donor_deletion_counts", "histology_tier1", exclusions, 2, median)
plot_ranks(results_by_donor, "snv_counts", "donor_deletion_counts", "histology_abbreviation", exclusions, 2, median)
plot_ranks(results_by_donor, "snv_counts", "donor_deletion_counts", "organ_system", exclusions, 2, median)
plot_ranks(results_by_donor, "snv_counts", "donor_deletion_counts", "age_bin", exclusions, 2, median)
plot_ranks(results_by_donor, "snv_counts", "donor_deletion_counts", "dcc_project_code", exclusions, 2, median)

