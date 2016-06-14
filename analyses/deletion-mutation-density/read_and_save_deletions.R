library(VariantAnnotation)
library(pcawg.common)

#Read in deletions from Chromosome 22 only 
chr2_rng = GRanges(seqnames="22", ranges=IRanges(start=0, end=51304566))
germline_deletions_chr22 = read_germline_deletions_from_vcf("~/Downloads/pcawg_data/germline_deletions/DEL.pcawg.rf.vcf.gz", chr2_rng)

save_deletions_to_r_file(germline_deletions_chr22, "~/Downloads/pcawg_data/germline_deletions/dels_chr22.Rdata")


germline_deletions = read_germline_deletions_from_vcf("~/Downloads/pcawg_data/germline_deletions/DEL.pcawg.rf.vcf.gz")

save_deletions_to_r_file(germline_deletions, "~/Downloads/pcawg_data/germline_deletions/dels.Rdata")