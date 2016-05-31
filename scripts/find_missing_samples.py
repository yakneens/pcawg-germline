import sys
import os
from os import listdir

results_path = "/shared/data/samples/vcf/annai_rtg/"
sample_lists = ["annai_germline_icgc_dkfz_ids_stripped.txt", "annai_germline_tcga_osdc_ids_stripped.txt", "annai_germline_icgc_osdc_ids_stripped.txt"]
results = listdir(results_path)
full_sample_list = []
for sample_list in sample_lists:
    f = open(sample_list, "r")
    full_sample_list = full_sample_list + f.read().splitlines()
    f.close()

full_sample_set = set(full_sample_list)

print len(full_sample_list)
print len(full_sample_set)

for el in full_sample_set:
    if el not in results:
        print el