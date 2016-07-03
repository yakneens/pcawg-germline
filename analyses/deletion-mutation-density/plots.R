ggplot(log(carrier_summary), aes(x=carrier_counts, y=Freq)) + geom_point() + labs(x="Log(# Carriers)", y="Log(# Deletions)")

#Plot test p-values over ratio of mean densities between carriers and non-carriers. 
plot_data = data.frame(log10(mean_ratios), -log10(count_df$wilcox_pvals))
colnames(plot_data) = c("x", "y")

plot_data_to_label = plot_data[sig_idx[queryHits(del_gene_overlaps)][cosmic_hit_idx],]
plot_data_to_label$gene_symbol = cosmic_gene_hits$SYMBOL
mapkapk5_row = c(plot_data[sig_idx[queryHits(del_gene_overlaps)][mapkapk5],], master_hit_df[mapkapk5,]$SYMBOL)
names(mapkapk5_row)[3] = "gene_symbol"
plot_data_to_label = rbind(plot_data_to_label, mapkapk5_row)

ggplot(plot_data) + aes(x=x, y=y) + geom_point(size=0.5) + 
  xlab(expression(paste(log[10](frac(carrier_density,non_carrier_density))))) +
  ylab(expression(paste(-log[10](p-value)))) +
  geom_abline(aes(intercept=-log10(significance_threshold), slope=0), col="red") +
  geom_point(data=plot_data_to_label, aes(x=x, y=y), size=3, col="red") +
  geom_label_repel(data=plot_data_to_label, 
                   aes(x, y, label=gene_symbol), 
                   fontface = 'bold', color = 'white',
                   box.padding = unit(0.25, "lines"),
                   point.padding = unit(0.5, "lines"), fill="red")






