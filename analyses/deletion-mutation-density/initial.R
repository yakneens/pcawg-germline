library(VariantAnnotation)
library(data.table)
library(pcawg.common)
library(GenomicFeatures)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(rtracklayer)
library(org.Hs.eg.db)
library(COSMIC.67)


#Read in metadata
sample_meta = get_pcawg_metadata("~/Downloads/pcawg_data/sample_metadata/pcawg_summary.tsv", split_multi_tumors = F)

sample_meta = keep_first_of_multi_tumors(sample_meta)

clinical_meta = get_clinical_metadata("~/Downloads/pcawg_data/clinical_metadata/pcawg_donor_clinical_May2016_v2.tsv", sample_meta)

hist_meta = get_histology_metadata("~/Downloads/pcawg_data/clinical_metadata/pcawg_specimen_histology_May2016_v2.tsv", sample_meta)

#Load germline deletions data
load("~/Downloads/pcawg_data/germline_deletions/dels.Rdata")
#load("~/Downloads/pcawg_data/germline_deletions/dels_chr22.Rdata")

set_deletion_range_ends(deletions)

deletion_genotypes = geno(deletions)$GT[,match(sample_meta$normal_wgs_aliquot_id, colnames(geno(deletions)$GT))] 


#Load somatic SNV data
load("~/Downloads/pcawg_data/snv_samples.Rdata")
#load("~/Downloads/pcawg_data/snv_samples_chr22.Rdata")

snv_samples = snv_samples[match(sample_meta$tumor_wgs_aliquot_id, names(snv_samples))]

#Detect and remove samples that are missing deletions or SNV calls
missing_sample_indices = c(which(is.na(colnames(deletion_genotypes))), which(is.na(names(snv_samples))))

sample_meta = sample_meta[-missing_sample_indices,]
clinical_meta = clinical_meta[-missing_sample_indices,]
hist_meta = hist_meta[-missing_sample_indices,]
deletion_genotypes = deletion_genotypes[,-missing_sample_indices]
snv_samples = snv_samples[-missing_sample_indices]

colnames(deletion_genotypes) = sample_meta$donor_unique_id
names(snv_samples) = sample_meta$donor_unique_id

#Compute matrix with 1 for heterozygous carriers of a deletion, 0 for non-carriers (hom ref), and NA for others
deletion_carrier_mask = get_het_carrier_mask(deletion_genotypes)

#Number of het carriers per deletion
carrier_counts = rowSums(deletion_carrier_mask, na.rm = T)

min_carrier_count = 5

deletion_filter = which(carrier_counts < min_carrier_count)

deletion_ranges = rowRanges(deletions)[-deletion_filter]
deletion_genotypes = deletion_genotypes[-deletion_filter,]
deletion_carrier_mask = t(deletion_carrier_mask[-deletion_filter,])

#Compute a matrix with count of SNVs overlying each deletion
snv_hits = do.call(rbind, lapply(snv_samples, function(x) as.table(findOverlaps(deletion_ranges, rowRanges(x)))))

colnames(snv_hits) = rownames(deletion_genotypes)

# Normalize each sample by number of SNVs in that sample
# Normalize each deletion by deletion width
del_widths =  width(ranges(deletion_ranges))
snv_counts = unlist(lapply(snv_samples, length))
  
normalized_snv_hits = snv_hits / snv_counts / del_widths[col(snv_hits)]

wilcox_pvals = c()

for(i in 1 : dim(normalized_snv_hits)[2]) {
  wilcox_pvals[i] = wilcox.test(normalized_snv_hits[deletion_carrier_mask[,i] == 1,i], normalized_snv_hits[deletion_carrier_mask[,i] == 0,i])$p.value
}

mean_counts = list()

#Calculate mean counts among carriers and non-carriers.
for(i in 1 : dim(normalized_snv_hits)[2]) {
  mean_counts[[i]] = c(mean(normalized_snv_hits[deletion_carrier_mask[,i] == 1,i], na.rm = T), mean(normalized_snv_hits[deletion_carrier_mask[,i] == 0,i], na.rm = T))
}

mean_counts = do.call(rbind, mean_counts)

count_df = as.data.frame(mean_counts)
colnames(count_df) = c("carriers", "non_carriers")
count_df$wilcox_pvals = pvals
count_df$widths = width(deletion_ranges)
count_df$carrier_counts = carrier_counts[-deletion_filter]

mean_difs = mean_counts[,1] - mean_counts[,2]

pre_flanks = flank(deletion_ranges, width(deletion_ranges) / 2, start = T)
pre_hits = do.call(rbind, lapply(snv_samples, function(x) as.table(findOverlaps(pre_flanks, rowRanges(x)))))

post_flanks = flank(deletion_ranges, width(deletion_ranges) / 2, start = F)
post_hits = do.call(rbind, lapply(snv_samples, function(x) as.table(findOverlaps(post_flanks, rowRanges(x)))))

flank_hits = pre_hits + post_hits
normalized_flank_hits = flank_hits / snv_counts / del_widths[col(flank_hits)]

mean_flank_counts = list()
for(i in 1 : dim(normalized_snv_hits)[2]) {
  mean_flank_counts[[i]] = c(mean(normalized_snv_hits[deletion_carrier_mask[,i] == 1,i], na.rm = T), mean(normalized_flank_hits[deletion_carrier_mask[,i] == 1,i], na.rm = T))
}

flank_wilcox_pvals = c()
flank_ttest_pvals = c()
flank_cors = c()
for(i in 1 : dim(normalized_snv_hits)[2]) {
  #flank_pvals[i] = t.test(na.omit(normalized_snv_hits[deletion_carrier_mask[,i] == 1,i]), na.omit(normalized_flank_hits[deletion_carrier_mask[,i] == 1,i]))$p.value
  flank_wilcox_pvals[i] = wilcox.test(na.omit(normalized_snv_hits[deletion_carrier_mask[,i] == 1,i]), na.omit(normalized_flank_hits[deletion_carrier_mask[,i] == 1,i]), paired = T)$p.value
  flank_ttest_pvals[i] = t.test(na.omit(normalized_snv_hits[deletion_carrier_mask[,i] == 1,i]), na.omit(normalized_flank_hits[deletion_carrier_mask[,i] == 1,i]), paired = T)$p.value
  flank_cors[i] = cor(na.omit(normalized_snv_hits[deletion_carrier_mask[,i] == 1,i]), na.omit(normalized_flank_hits[deletion_carrier_mask[,i] == 1,i]))
}

mean_flank_counts = do.call(rbind, mean_flank_counts)
mean_flank_difs = mean_flank_counts[,1] - mean_flank_counts[,2]

flank_df = as.data.frame(mean_flank_counts)
colnames(flank_df) = c("dels", "flanks")
flank_df$widths = width(deletion_ranges)
flank_df$carrier_counts = carrier_counts[-deletion_filter]
flank_df$wilcox_pvals = flank_wilcox_pvals
flank_df$ttest_pvals = flank_ttest_pvals
flank_df$chr = as.character(seqnames(deletion_ranges))


txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
import.chain("~/Downloads/pcawg_data/external_data/hg19ToGRCh37.over.chain")
genes_grch37 = liftOver(genes(txdb), ch)

num_hypotheses = dim(normalized_snv_hits)[2]
significance_threshold = 0.05 / num_hypotheses

#sig_idx = intersect(which(count_df$wilcox_pvals < 0.01), which(flank_df$wilcox_pvals < 0.01))
sig_idx = which(count_df$wilcox_pvals < significance_threshold)
#sig_idx = which(flank_df$wilcox_pvals < significance_threshold)

sig_dels = deletion_ranges[sig_idx]
del_gene_overlaps = findOverlaps(sig_dels, genes_grch37)
gene_hits = genes_grch37[subjectHits(del_gene_overlaps)]

gene_hit_details = select(org.Hs.eg.db, unlist(gene_hits)$gene_id, c("SYMBOL", "GENENAME"))

data("cgc_67", package="COSMIC.67")
if(any(unlist(gene_hits)$gene_id %in% cgc_67$ENTREZID)){
  cosmic_gene_hits = gene_hit_details[na.omit(match(cgc_67$ENTREZID, unlist(gene_hits)$gene_id)),] 
}

master_hit_df = cbind(flank_df[sig_idx[queryHits(del_gene_overlaps)],], count_df[sig_idx[queryHits(del_gene_overlaps)],], gene_hit_details)
