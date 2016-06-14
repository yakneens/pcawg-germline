library(VariantAnnotation)
library(pcawg.common)

#Read in sample metadata
sample_meta = get_pcawg_metadata("~/Downloads/pcawg_data/sample_metadata/pcawg_summary.tsv")

#Read in SNV samples
snv_samples = read_snv_samples_from_vcf("~/Downloads/pcawg_data/somatic_snv_mnv/", sample_meta)

#Save SNV samples in R format
save_snv_samples(snv_samples, "~/Downloads/pcawg_data/snv_samples.Rdata")

