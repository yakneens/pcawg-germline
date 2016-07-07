library(VariantAnnotation)
library(pcawg.common)

#Read in sample metadata
sample_meta = get_pcawg_metadata("~/Downloads/pcawg_data/sample_metadata/pcawg_summary.tsv")


#Read in snvs from Chromosome 22 only 
chr22_rng = GRanges(seqnames="22", ranges=IRanges(start=0, end=51304566))
snv_samples_chr22 = read_snv_samples_from_vcf("~/Downloads/pcawg_data/somatic_snv_mnv/", sample_meta, chr22_rng)

#Save SNV samples in R format
save_snv_samples(snv_samples_chr22, "~/Downloads/pcawg_data/snv_samples_chr22.Rdata")

#Read in SNV samples
snv_samples = read_snv_samples_from_vcf("~/Downloads/pcawg_data/somatic_snv_mnv/", sample_meta)

#Save SNV samples in R format
save_snv_samples(snv_samples, "~/Downloads/pcawg_data/snv_samples.Rdata")

