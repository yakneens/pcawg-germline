### Carrier VS Non-Carrier Plots

# Number of carriers vs number of deletions
ggplot(log(carrier_summary), aes(x=carrier_counts, y=N)) + geom_point() + labs(x="Log(# Carriers)", y="Log(# Deletions)")

#Plot test p-values over ratio of mean densities between carriers and non-carriers. 
plot_data = results_by_del[,list(log10(c_vs_nc_mean_ratios), -log10(c_vs_nc_wilcox_pvals))]
setnames(plot_data, c("x", "y"))

# Add labels for COSMIC Cancer Census genes.
plot_data_to_label = plot_data[sig_idx[queryHits(del_gene_overlaps)][cosmic_hit_idx],]
plot_data_to_label$gene_symbol = cosmic_gene_hits$SYMBOL
additional_gene_rows = cbind(plot_data[sig_idx[queryHits(del_gene_overlaps)][additional_interesting_gene_idx],], master_hit_df[additional_interesting_gene_idx,.(SYMBOL)])
setnames(additional_gene_rows, c("x", "y", "gene_symbol"))

plot_data_to_label = rbind(plot_data_to_label, additional_gene_rows)

ggplot(plot_data) + aes(x=x, y=y) + geom_point(size=0.5) + 
  xlab(expression(paste(log[10](frac(carrier_density,non_carrier_density))))) +
  ylab(expression(paste(-log[10](p-value)))) +
  geom_abline(aes(intercept=-log10(significance_threshold_dels), slope=0), col="red") +
  geom_point(data=plot_data_to_label, aes(x=x, y=y), size=3, col="red") +
  geom_label_repel(data=plot_data_to_label, 
                   aes(x, y, label=gene_symbol), 
                   fontface = 'bold', color = 'white',
                   box.padding = unit(0.25, "lines"),
                   point.padding = unit(0.5, "lines"), fill="red")



results_by_del[, `:=`(chr=as.vector(seqnames(deletion_ranges)), width_bin=findInterval(del_widths, c(500,1000,2000,5000,10000,100000,250000,1000000)), carrier_count_bin=findInterval(carrier_counts, seq(0,300,25)))]

# Carrier vs Non-Carrier mean SNV counts
ggplot(results_by_del) + geom_point(aes(x=-log10(carriers), y=-log10(non_carriers)))  +  
  geom_point(aes(x=-log10(mean(carriers)), y=-log10(mean(non_carriers))), col="red", size=5)

# Carrier vs Non-Carrier SNV density plot
ggplot(melt(results_by_del, measure.vars = c("carriers", "non_carriers")), aes(-log10(value), y=..scaled..,fill=variable)) + geom_density(alpha=0.6) + xlab("-log10(normalized_snv_count)") + ylab("Density")

# Histogram of Mann Whitney U Test p-values among carriers and non-carriers
ggplot(results_by_del, aes(c_vs_nc_wilcox_pvals)) + geom_histogram(col=I("black"), fill=I("grey"))

# QQ Plot for the distribution of Mann Whitney p-values
ggplot(results_by_del) + stat_qq(aes(sample=c_vs_nc_wilcox_pvals), distribution=stats::qunif) + geom_abline(aes(slope=1, intercept=0))

# Ratio of mean SNV counts in carriers vs non-carriers for 100 smallest p-values
top_100_pval = order(results_by_del$c_vs_nc_wilcox_pvals)[1:100]
ggplot(results_by_del[top_100_pval,]) + geom_point(aes(x=seq(1,100),y=log10(carriers/non_carriers))) + geom_abline(aes(intercept=0,slope=0), col="red") + xlab("Index") + ylab(expression(paste(log[10](frac(carrier_density,non_carrier_density))))) + xlim(0,110) + ylim(-2,5)


### Deletions vs Flanks Plots

# Deletion vs Flank mean SNV counts
ggplot(results_by_del) + geom_point(aes(x=-log10(carriers), y=-log10(flanks))) +
  geom_point(aes(x=-log10(mean(carriers)), y=-log10(mean(flanks))), col="red", size=5)

# Density of SNV counts for deletions and flanks
ggplot(melt(results_by_del, measure.vars = c("carriers", "flanks")), aes(-log10(value), y=..scaled..,fill=variable)) + geom_density(alpha=0.6) + xlab("-log10(normalized_snv_count)") + ylab("Density")

# Histogram of Mann Whitney U Test among deletions and flanks
ggplot(results_by_del, aes(d_vs_f_wilcox_pvals)) + geom_histogram(col=I("black"), fill=I("grey"))

# QQ Plot for the p-value distribution 
ggplot(results_by_del) + stat_qq(aes(sample=d_vs_f_wilcox_pvals), distribution=stats::qunif) + geom_abline(aes(slope=1, intercept=0))
