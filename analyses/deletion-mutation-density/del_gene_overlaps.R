txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
ch = import.chain("~/Downloads/pcawg_data/external_data/hg19ToGRCh37.over.chain")
genes_grch37 = liftOver(genes(txdb), ch)

get_gene_overlaps <- function(dels_index){
  sel_dels = deletion_ranges[dels_index]
  del_gene_overlaps = findOverlaps(sel_dels, genes_grch37)
  gene_hits = genes_grch37[subjectHits(del_gene_overlaps)]
  gene_hit_details = select(org.Hs.eg.db, unlist(gene_hits)$gene_id, c("SYMBOL", "GENENAME"))
  
  deletion_hits = deletion_ranges[dels_index[queryHits(del_gene_overlaps)]]
  
  master_hit_df = cbind(results_by_del[dels_index[queryHits(del_gene_overlaps)],], gene_hit_details)
  master_hit_df[, `:=` (del_start = start(ranges(deletion_hits)), del_end = end(ranges(deletion_hits)), dels_index = dels_index[queryHits(del_gene_overlaps)])]
  
  return(master_hit_df)
  
}

get_cosmic_hits <- function(hits_dt) {
  
  data("cgc_67", package="COSMIC.67")
  if(any(hits_dt$ENTREZID %in% cgc_67$ENTREZID)){
    cosmic_gene_hits = hits_dt[na.omit(match(cgc_67$ENTREZID, hits_dt$ENTREZID)),] 
  }
  
  return(cosmic_gene_hits)
}

#sig_idx = intersect(which(count_df$wilcox_pvals < 0.01), which(flank_df$wilcox_pvals < 0.01))
sig_idx = which(results_by_del$c_vs_nc_wilcox_pvals < significance_threshold_dels)
#sig_idx = which(results_by_del$d_vs_f_wilcox_pvals < significance_threshold_dels)

c_vs_nc_hits = get_gene_overlaps(sig_idx)

#Hits with top 10 p-vals
top_10_pv = c_vs_nc_hits[order(c_vs_nc_wilcox_pvals)[1:10], .(c_vs_nc_wilcox_pvals, SYMBOL, GENENAME, c_vs_nc_mean_ratios, carrier_counts, chr, del_start, del_end)]
write.table(top_10_pv, sep=",", col.names=T)

for (el in c_vs_nc_hits[order(c_vs_nc_wilcox_pvals)[9], dels_index]){
  print(c_vs_nc_hits[dels_index == el, SYMBOL])
  print("Dels")
  write.table(normalized_snv_hits[which(deletion_carrier_mask[,el] == 1), el][which(normalized_snv_hits[which(deletion_carrier_mask[,el] == 1), el] > 0)], sep=",")
  print("Flanks")
  write.table(normalized_flank_hits[which(deletion_carrier_mask[,el] == 1), el][which(normalized_flank_hits[which(deletion_carrier_mask[,el] == 1), el] > 0)], sep=",")
  barplot(table(results_by_donor[which(deletion_carrier_mask[,el] == 1), histology_tier2]), 
          main=paste(c_vs_nc_hits[dels_index == el, SYMBOL], "-", c_vs_nc_hits[dels_index == el, carrier_counts], "carriers"), 
          las=2)
}


#Hits with top 10 mean-ratios
top_10_mr = c_vs_nc_hits[order(-c_vs_nc_mean_ratios)[1:10], .(c_vs_nc_wilcox_pvals, SYMBOL, GENENAME, c_vs_nc_mean_ratios, carrier_counts, chr, del_start, del_end)]
write.table(top_10_mr, sep=",", col.names=T)

for (el in c_vs_nc_hits[order(-c_vs_nc_mean_ratios)[10:10], dels_index]){
  print(c_vs_nc_hits[dels_index == el, SYMBOL])
  print("Dels")
  write.table(normalized_snv_hits[which(deletion_carrier_mask[,el] == 1), el][which(normalized_snv_hits[which(deletion_carrier_mask[,el] == 1), el] > 0)], sep=",")
  print("Flanks")
  write.table(normalized_flank_hits[which(deletion_carrier_mask[,el] == 1), el][which(normalized_flank_hits[which(deletion_carrier_mask[,el] == 1), el] > 0)], sep=",")
  barplot(table(results_by_donor[which(deletion_carrier_mask[,el] == 1), histology_tier2]), 
          main=paste(c_vs_nc_hits[dels_index == el, SYMBOL], "-", c_vs_nc_hits[dels_index == el, carrier_counts], "carriers"), 
          las=2)
}


sig_idx = which(results_by_del$d_vs_f_wilcox_pvals < significance_threshold_dels)

d_vs_f_hits = get_gene_overlaps(sig_idx)

#Hits with top 10 p-vals
top_10_pv = d_vs_f_hits[order(d_vs_f_wilcox_pvals)[1:10], .(d_vs_f_wilcox_pvals, SYMBOL, GENENAME, d_vs_f_mean_ratios, carrier_counts, chr, del_start, del_end)]
write.table(top_10_pv, sep=",", col.names=T)

for (el in d_vs_f_hits[order(d_vs_f_wilcox_pvals)[1:10], dels_index]){
  print(d_vs_f_hits[dels_index == el, SYMBOL])
  print("Dels")
  write.table(normalized_snv_hits[which(deletion_carrier_mask[,el] == 1), el][which(normalized_snv_hits[which(deletion_carrier_mask[,el] == 1), el] > 0)], sep=",")
  print("Flanks")
  write.table(normalized_flank_hits[which(deletion_carrier_mask[,el] == 1), el][which(normalized_flank_hits[which(deletion_carrier_mask[,el] == 1), el] > 0)], sep=",")
  barplot(table(results_by_donor[which(deletion_carrier_mask[,el] == 1), histology_tier2]), 
          main=paste(d_vs_f_hits[dels_index == el, SYMBOL], "-", d_vs_f_hits[dels_index == el, carrier_counts], "carriers"), 
          las=2)
}


#Hits with top 10 mean-ratios
top_10_mr = d_vs_f_hits[order(d_vs_f_mean_ratios)[1:10], .(d_vs_f_wilcox_pvals, SYMBOL, GENENAME, d_vs_f_mean_ratios, carrier_counts, chr, del_start, del_end)]
write.table(top_10_mr, sep=",", col.names=T)

for (el in d_vs_f_hits[order(d_vs_f_mean_ratios)[6], dels_index]){
  print(d_vs_f_hits[dels_index == el, SYMBOL])
  print("Dels")
  write.table(normalized_snv_hits[which(deletion_carrier_mask[,el] == 1), el][which(normalized_snv_hits[which(deletion_carrier_mask[,el] == 1), el] > 0)], sep=",")
  print("Flanks")
  write.table(normalized_flank_hits[which(deletion_carrier_mask[,el] == 1), el][which(normalized_flank_hits[which(deletion_carrier_mask[,el] == 1), el] > 0)], sep=",")
  barplot(table(results_by_donor[which(deletion_carrier_mask[,el] == 1), histology_tier2]), 
          main=paste(d_vs_f_hits[dels_index == el, SYMBOL], "-", d_vs_f_hits[dels_index == el, carrier_counts], "carriers"), 
          las=2)
}

cosmic_hit_idx = get_cosmic_hits(c_vs_nc_hits)
additional_interesting_genes = c("MAPKAPK5", "AVL9", "ADAM3A", "POLE2", "NR5A2", "MCPH1", "CCNY")
additional_interesting_gene_idx = which(c_vs_nc_hits$SYMBOL %in% additional_interesting_genes)


for (el in c_vs_nc_hits[adam3a, dels_index]){
  print(c_vs_nc_hits[dels_index == el, SYMBOL])
  print("Dels")
  write.table(normalized_snv_hits[which(deletion_carrier_mask[,el] == 1), el][which(normalized_snv_hits[which(deletion_carrier_mask[,el] == 1), el] > 0)], sep=",")
  print("Flanks")
  write.table(normalized_flank_hits[which(deletion_carrier_mask[,el] == 1), el][which(normalized_flank_hits[which(deletion_carrier_mask[,el] == 1), el] > 0)], sep=",")
  barplot(table(results_by_donor[which(deletion_carrier_mask[,el] == 1), histology_tier2]), 
          main=paste(c_vs_nc_hits[dels_index == el, SYMBOL], "-", c_vs_nc_hits[dels_index == el, carrier_counts], "carriers"), 
          las=2)
}


# hs_chrs = as(seqinfo(BSgenome.Hsapiens.1000genomes.hs37d5), "GRanges")
# hs_chrs = keepSeqlevels(hs_chrs, c(1:22, "X", "Y"))
# 
# 
# g_bins = tileGenome(seqlengths(hs_chrs), tilewidth=1000, cut.last.tile.in.chrom=T)
# g_bins_100kb = tileGenome(seqlengths(hs_chrs), tilewidth=100000, cut.last.tile.in.chrom=T)
# g_bins_1mb = tileGenome(seqlengths(hs_chrs), tilewidth=1000000, cut.last.tile.in.chrom=T)
# a = as.table(findOverlaps(g_bins_1mb, rowRanges(snv_samples[[1]])))
# a_idx = seq(a)[which(a >0)]
# a = a[which(a >0)]
# 
# b = snv_hits[1,] / del_widths * 1000
# b2 = b[which(deletion_carrier_mask[1,] == 1)]
# deletion_bins = findOverlaps(deletion_ranges, g_bins_1mb)
# b2_bins = subjectHits(deletion_bins[which(deletion_carrier_mask[1,] == 1),])
# plot.new()
# plot_snp_density = data.frame(seq(length(a)),a)
# colnames(plot_snp_density) = c("bin", "density")
# ggplot(data=plot_snp_density, aes(x=bin, y=density)) + geom_point() +  geom_smooth(span=0.00000000001)
# 
# smoothScatter(seq(length(a)),a)
# spl = smooth.spline(seq(length(a)), a, spar=0.1)
# lines(spl, col="green", lwd=2)
# points(b2_bins[which(b2 >0)], b2[which(b2>0)], col = "red")


