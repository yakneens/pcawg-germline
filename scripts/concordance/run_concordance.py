import os, fnmatch
from subprocess import call
import gzip, re

rootdir = '/icgc/pcawg/analysis/iakhnin/germline_genotype_concordance'
raw_dir = rootdir + '/raw/'
processed_dir = rootdir + '/processed/'
results_dir = rootdir + '/results/'
header_file_name = "header.txt"
snp6_vcf_dir = '/icgc/pcawg/analysis/waszak/tcga_snp6/vcf_files/'
snp6_vcf_suffix = '.affy-snp6.cqc_0.4.qcr_0.86.birdseed-dev_2.6.cr_0.95.chr1_22.bi.snv.polymorphic.vcf.gz'
#rtg_vcf_dir = "/icgc/pcawg/project/results/germline_variants/annai/tcga/"
rtg_vcf_dir = "/icgc/pcawg/analysis/iakhnin/test/"
snp6_to_rtg_mapping_file = "/icgc/pcawg/analysis/iakhnin/PCAWG_SNP6_uuid_to_annai_vcf_uuid.txt"
header_size = 122

def find_sample_files(pattern, path):
    print "Looking for RTG samples"
    result = {}
    for root, dirs, files in os.walk(path):
        for name in files:
            if fnmatch.fnmatch(name, pattern):
                current_sample_filename = os.path.join(root, name)
                file_uuid = os.path.basename(root)
                result[file_uuid] = current_sample_filename
    print str(len(result)) + " samples found"
    return result

def get_snp6_to_rtg_mappings(mapping_file_path):
    snp6_to_rtg = {}
    rtg_to_snp6 = {}
    f = open(mapping_file_path, "r")
    for line in f:
        rec = line.strip().split(",")
        snp6_to_rtg[rec[0]] = rec[1]
        rtg_to_snp6[rec[1]] = rec[0]
    return snp6_to_rtg, rtg_to_snp6

def get_results_file_name(snp6_uuid, rtg_uuid):
    return snp6_uuid + "_" + rtg_uuid + "_snp6_rtg_concordance_test.txt"


def generate_pbs():
    result_list = find_sample_files('*.vcf.gz', rtg_vcf_dir)
    snp6_to_rtg, rtg_to_snp6 = get_snp6_to_rtg_mappings(snp6_to_rtg_mapping_file)
    for rtg_file_uuid in result_list.keys():
        if rtg_to_snp6.has_key(rtg_file_uuid):              
            sample_file = result_list[rtg_file_uuid]
            
            p = open("/icgc/pcawg/analysis/iakhnin/germline_genotype_concordance/submit_concordance.pbs", "w")
            p.write("#!/bin/bash\n")
            p.write("#PBS -o " + raw_dir + rtg_file_uuid + "/pbs.out\n")
            p.write("#PBS -e " + raw_dir + rtg_file_uuid + "/pbs.err\n")
            p.write("#PBS -l nodes=1:ppn=8,mem=8gb\n")
            p.write("#PBS -l walltime=00::15:00\n")
            p.write("#PBS -N siakhnin_pcawg_snp6\n")
            p.write("#PBS -M iakhnin@embl.de\n")
            p.write("#PBS -m e\n")
            p.write("sample_uuid = $(python /icgc/pcawg/analysis/iakhnin/pcawg-germline/process_header.py " + sample_file + " " + rtg_file_uuid + ")\n")
            header_file_path = raw_dir + rtg_file_uuid + "/" + header_file_name 
            
            snp6_uuid = rtg_to_snp6[rtg_file_uuid]
            results_file_name = raw_dir + get_results_file_name(snp6_uuid, rtg_file_uuid) 
            
            snp6_file_path = snp6_vcf_dir + snp6_uuid + snp6_vcf_suffix
            p.write("/icgc/pcawg/analysis/iakhnin/bcftools-1.2/bcftools reheader -h " + header_file_path + " -o " + raw_dir + rtg_file_uuid + "/$sample_uuid.vcf.gz " + sample_file + "\n")
            p.write("tabix -f -p vcf " + raw_dir + rtg_file_uuid + "/$sample_uuid.vcf.gz\n")
            p.write("vcf-compare -a -g -m " + snp6_uuid + ":$sample_uuid " + snp6_file_path + " " + raw_dir + rtg_file_uuid + "/$sample_uuid.vcf.gz > " + raw_dir + rtg_file_uuid  + "/" + results_file_name + "\n")
            p.write("mv " + raw_dir + rtg_file_uuid + "/*.vcf.gz* " + processed_dir + "\n")
            p.write("mv " + raw_dir + rtg_file_uuid + "/" + results_file_name + " " + results_dir + "\n")
            p.write("rm -rf" + raw_dir + rtg_file_uuid + "\n")
            
            #call("qsub /icgc/pcawg/analysis/iakhnin/germline_genotype_concordance/submit_concordance.pbs", shell=True)         
    
generate_pbs()   
