#!/usr/bin/env Rscript

library(VariantAnnotation)
library(data.table)
library(GenomicFeatures)
library(progress)
library(logging)
library(purrr)

selected_chrom = 14
donor_meta_path = "~/Downloads/pcawg_data/del_density/input_data/donor_meta.RData"
deletion_ranges_path = "~/Downloads/pcawg_data/del_density/input_data/deletion_ranges.RData"
snv_ranges_path = "~/Downloads/pcawg_data/del_density/input_data/snv_ranges.RData"
carrier_mask_path = "~/Downloads/pcawg_data/del_density/input_data/deletion_carrier_mask.RData"
result_path = "~/Downloads/pcawg_data/del_density/"
non_carriers = FALSE



load(donor_meta_path)
load(deletion_ranges_path)
load(snv_ranges_path)


loginfo("Loaded input data")


#deletion_filter = which(as.character(seqnames(deletion_ranges)) != selected_chrom)
num_bins = 250
margin_width = 50000
start_c = 106030296 - margin_width
end_c = 107011401 + margin_width
window_width = end_c - start_c
chrom = "14"

margin_width = 100000
start_c = 56742258 - margin_width
end_c = 57178108 + margin_width
window_width = end_c - start_c
chrom = "7"
num_bins = 200

region = GRanges(IRanges(start=start_c, end=end_c), seqnames=chrom)

region_tiles = tile(region, n = as.integer(num_bins))
pb = progress_bar$new(format=":current/:total [:bar] :percent :elapsed :eta",total = length(snv_ranges))
res = as.data.table(t(snv_ranges %>% map_df(~{pb$tick();overlaps = countOverlaps(unlist(region_tiles), .);})))[,donor_id := donor_meta$donor_unique_id]
res[,histo := donor_meta[donor_unique_id == res$donor_id]$histology_tier2]
res[,project := donor_meta[donor_unique_id == res$donor_id]$dcc_project_code]
res[histo != "Lymphoid", hist := "Other"]
res[histo != "Lymphoid", project := "Other"]

res[!histo %in% c("Colon/Rectum", "Skin"), histo := "Other"]
res[!histo %in% c("Colon/Rectum", "Skin"), project := "Other"]

res_m = melt(res, id.vars = c("donor_id", "histo", "project"), variable.name = "bin_index", value.name = "snv_count")
res_m$histo = factor(res_m$histo, levels=c("Colon/Rectum", "Skin", "Other"))
res_m$project = factor(res_m$project, levels=c("COAD-US", "READ-US", "MELA-AU", "SKCM-US", "Other"))
ggplot(res_m[,mean(snv_count), by=c("bin_index", "histo")], aes(x=as.integer(gsub("V", "", bin_index)), y=V1, fill=histo)) + 
  geom_bar(stat="identity", position="stack") +
  geom_vline(xintercept=32, linetype="dashed") + 
  geom_vline(xintercept=168, linetype="dashed") + 
  facet_grid(~histo) +
  xlab("Bin") + 
  ylab("Mean SNV Load")

ggplot(res_m[,mean(snv_count), by=c("bin_index", "project")], aes(x=as.integer(gsub("V", "", bin_index)), y=V1, fill=project)) + 
  geom_bar(stat="identity", position="stack") +
  geom_vline(xintercept=32, linetype="dashed") + 
  geom_vline(xintercept=168, linetype="dashed") + 
  facet_wrap(~project) +
  xlab("Bin") + 
  ylab("Mean SNV Load")

