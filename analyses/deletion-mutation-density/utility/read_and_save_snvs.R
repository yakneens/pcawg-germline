library(VariantAnnotation)
library(pcawg.common)

#Read in sample metadata
sample_meta = get_pcawg_metadata("~/Downloads/pcawg_data/sample_metadata/pcawg_summary.tsv", "~/Downloads/pcawg_data/sample_metadata/PCAWG Excluded Donors%2FSamples - Excluded_donors_2016_08_30.tsv")


#Read in snvs from Chromosome 22 only 
chr22_rng = GRanges(seqnames="22", ranges=IRanges(start=0, end=51304566))
snv_samples_chr22 = read_snv_samples_from_vcf("~/Downloads/pcawg_data/final_consensus_somatic_snv/final_consensus_12oct_passonly/snv_mnv/", sample_meta, chr22_rng)

#Save SNV samples in R format
save_snv_samples(snv_samples_chr22, "~/Downloads/pcawg_data/snv_samples_chr22.Rdata")

#Read in SNV samples

snv_samples = read_snv_samples_from_vcf("~/Downloads/pcawg_data/final_consensus_somatic_snv/final_consensus_12oct_passonly/snv_mnv/", sample_meta)

#Save SNV samples in R format
save_snv_samples(snv_samples, "~/Downloads/pcawg_data/snv_samples.Rdata")

germline_snv_samples = read_germline_snvs_from_vcf("~/Downloads/pcawg_data/germline_snv/pcawg2727.chr1_22.1000gp_phase3_shapeit2_mvncall_integrated_v5b.gmaf1.genotypes.release_180216.vcf.gz",chr22_rng)
save_snv_samples(germline_snv_samples, "~/Downloads/pcawg_data/germline_snv_samples.Rdata")
