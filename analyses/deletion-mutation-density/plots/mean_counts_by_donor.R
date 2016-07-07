in_list = c("READ-US", "READ-US", "MELA-AU", "SKCM-US", "UCEC-US")

# Ratio of deletion vs flank mean SNV counts per donor against p-value, colored by project
ggplot(results_by_donor) + 
  aes(x=log10(mean_ratios), y=-log10(wilcox_pvals), col=dcc_project_code) + 
  geom_point(size=0.5) + 
  xlab(expression(paste(log[10](frac(deletion_density,flank_density))))) +
  ylab(expression(paste(-log[10](p-value)))) +
  geom_abline(aes(intercept=-log10(significance_threshold_donors), slope=0), col="red") +
  scale_fill_discrete(guide = FALSE) + theme(legend.position="none") +
  geom_label_repel(data=results_by_donor[wilcox_pvals < significance_threshold_donors], 
                   aes(label=dcc_project_code, fill=dcc_project_code), 
                   fontface = 'bold', color = 'white',
                   box.padding = unit(0.25, "lines"),
                   point.padding = unit(0.5, "lines"))

# Ratio of deletion vs flank mean SNV counts per donor against p-value, colored by Histology Tier 2
ggplot(results_by_donor) + 
  aes(x=log10(mean_ratios), y=-log10(wilcox_pvals), col=histology_tier2) + 
  geom_point(size=0.5) + 
  xlab(expression(paste(log[10](frac(deletion_density,flank_density))))) +
  ylab(expression(paste(-log[10](p-value)))) +
  geom_abline(aes(intercept=-log10(significance_threshold_donors), slope=0), col="red") +
  scale_fill_discrete(guide = FALSE) + theme(legend.position="none") +
  geom_label_repel(data=results_by_donor[wilcox_pvals < significance_threshold_donors], 
                   aes(label=histology_tier2, fill=histology_tier2), 
                   fontface = 'bold', color = 'white',
                   box.padding = unit(0.25, "lines"),
                   point.padding = unit(0.5, "lines"))

# Mean SNV counts per donor against deletion/flank p-value
ggplot(results_by_donor) +
  geom_point(aes(x=log10(snv_counts), y=-log10(wilcox_pvals))) +
  geom_abline(aes(intercept=-log10(significance_threshold_donors), slope=0), col="red") +
  geom_point(data=results_by_donor[wilcox_pvals < significance_threshold_donors], aes(x=log10(snv_counts), y=-log10(wilcox_pvals)), col="red")
  
  # geom_point(data=plot_data_to_label, aes(x=x, y=y), size=3, col="red") +
  # geom_label_repel(data=plot_data_to_label, 
  #                  aes(x, y, label=gene_symbol), 
  #                  fontface = 'bold', color = 'white',
  #                  box.padding = unit(0.25, "lines"),
  #                  point.padding = unit(0.5, "lines"), fill="red")

