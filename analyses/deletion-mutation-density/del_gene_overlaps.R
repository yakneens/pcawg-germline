txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
ch = import.chain("~/Downloads/pcawg_data/external_data/hg19ToGRCh37.over.chain")
genes_grch37 = liftOver(genes(txdb), ch)

#sig_idx = intersect(which(count_df$wilcox_pvals < 0.01), which(flank_df$wilcox_pvals < 0.01))
sig_idx = which(count_df$wilcox_pvals < significance_threshold_dels)
#sig_idx = which(flank_df$wilcox_pvals < significance_threshold_dels)

sig_dels = deletion_ranges[sig_idx]
del_gene_overlaps = findOverlaps(sig_dels, genes_grch37)
gene_hits = genes_grch37[subjectHits(del_gene_overlaps)]

gene_hit_details = select(org.Hs.eg.db, unlist(gene_hits)$gene_id, c("SYMBOL", "GENENAME"))

data("cgc_67", package="COSMIC.67")
if(any(unlist(gene_hits)$gene_id %in% cgc_67$ENTREZID)){
  cosmic_gene_hits = gene_hit_details[na.omit(match(cgc_67$ENTREZID, unlist(gene_hits)$gene_id)),] 
}

master_hit_df = cbind(flank_df[sig_idx[queryHits(del_gene_overlaps)],], count_df[sig_idx[queryHits(del_gene_overlaps)],], gene_hit_details)

cosmic_hit_idx = match(cosmic_gene_hits$ENTREZID, master_hit_df$ENTREZID)
cosmic_hit_deletions = rownames(master_hit_df[cosmic_hit_idx,])

mapkapk5 = which(gene_hit_details$SYMBOL == "MAPKAPK5")
mapkapk5_del = rownames(master_hit_df[mapkapk5,])

interesting_deletions = c(cosmic_hit_deletions, mapkapk5_del)

for(i in interesting_deletions){
  print(table(na.omit(snv_hits[deletion_carrier_mask[, i] == 1,i])))
  print(table(na.omit(snv_hits[deletion_carrier_mask[, i] == 0,i])))
}

hs_chrs = as(seqinfo(BSgenome.Hsapiens.1000genomes.hs37d5), "GRanges")
hs_chrs = keepSeqlevels(hs_chrs, c(1:22, "X", "Y"))


g_bins = tileGenome(seqlengths(hs_chrs), tilewidth=1000, cut.last.tile.in.chrom=T)
g_bins_100kb = tileGenome(seqlengths(hs_chrs), tilewidth=100000, cut.last.tile.in.chrom=T)
g_bins_1mb = tileGenome(seqlengths(hs_chrs), tilewidth=1000000, cut.last.tile.in.chrom=T)
a = as.table(findOverlaps(g_bins_1mb, rowRanges(snv_samples[[1]])))
a_idx = seq(a)[which(a >0)]
a = a[which(a >0)]

b = snv_hits[1,] / del_widths * 1000
b2 = b[which(deletion_carrier_mask[1,] == 1)]
deletion_bins = findOverlaps(deletion_ranges, g_bins_1mb)
b2_bins = subjectHits(deletion_bins[which(deletion_carrier_mask[1,] == 1),])
plot.new()
plot_snp_density = data.frame(seq(length(a)),a)
colnames(plot_snp_density) = c("bin", "density")
ggplot(data=plot_snp_density, aes(x=bin, y=density)) + geom_point() +  geom_smooth(span=0.00000000001)

smoothScatter(seq(length(a)),a)
spl = smooth.spline(seq(length(a)), a, spar=0.1)
lines(spl, col="green", lwd=2)
points(b2_bins[which(b2 >0)], b2[which(b2>0)], col = "red")
