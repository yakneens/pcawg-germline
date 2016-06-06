# v = readVcf(vcf.files[i], genome = "hs37d5")
# vmat = data.frame("pid" = meta$donor_unique_id[i], 
#                   "chr" = as.character(seqnames(rowRanges(v))),
#                   "pos" = start(rowRanges(v)),
#                   "ref" = as.character(rowRanges(v)$REF),
#                   "alt" = as.character(unlist(rowRanges(v)$ALT)),
#                   "ncaller" = info(v)$NumCallers,
#                   "caller" = unlist(lapply(info(v)$Caller, function(x) paste0(x, collapse = "-"))),
#                   "maf" = unlist(info(v)$"1000genomes_AF"), stringsAsFactors = FALSE)
# vmat$maf[is.na(vmat$maf)] = 0
library(VariantAnnotation)
v = readVcf("~/Downloads/preliminary_final_release/snv_mnv/0009b464-b376-4fbc-8a56-da538269a02f.annotated.snv_mnv.vcf.gz", genome = "hs37d5")

vmat = data.frame("pid" = "0009b464-b376-4fbc-8a56-da538269a02f", 
  "chr" = as.character(seqnames(rowRanges(v))),
  "pos" = start(rowRanges(v)),
  "ref" = as.character(rowRanges(v)$REF),
  "alt" = as.character(unlist(rowRanges(v)$ALT)),
  "ncaller" = info(v)$NumCallers,
  "caller" = unlist(lapply(info(v)$Caller, function(x) paste0(x, collapse = "-"))),
  "maf" = unlist(info(v)$"1000genomes_AF"), stringsAsFactors = FALSE)
vmat$maf[is.na(vmat$maf)] = 0

#Read in sample metadata
sample_meta = read.table("~/Downloads/pcawg_data/sample_metadata/pcawg_summary.tsv", header=TRUE, sep="\t")

#Read in deletions from Chromosome 22 only 
rng = GRanges(seqnames="22", ranges=IRanges(start=0, end=51304566))
tab = TabixFile("~/Downloads/pcawg_data/DEL.pcawg.rf.vcf.gz")
dels = readVcf(tab, "hs37d5", param=rng)

sub_dels = data.frame(start(dels), info(dels)$END, geno(dels)$GT)

normal_wgs_aliquots = colnames(dels)
tumor_aliquots = sample_meta[which(sample_meta$normal_wgs_aliquot_id %in% normal_wgs_aliquots),]$tumor_wgs_aliquot_id

snp_samples = list()

for(aliquot in tumor_aliquots){
  tab = TabixFile(paste("~/Downloads/preliminary_final_release/snv_mnv/",aliquot,".annotated.snv_mnv.vcf.gz"))
  snp_samples[aliquot] = readVcf(tab, "hs37d5", param=rng)
}
