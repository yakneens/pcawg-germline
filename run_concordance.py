import os, fnmatch
from subprocess import call
import gzip, file, re

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
    result = []
    for root, dirs, files in os.walk(path):
        for name in files:
            if fnmatch.fnmatch(name, pattern):
                current_sample_filename = os.path.join(root, name)
                result.append(current_sample_filename)
    print len(result) + " samples found"
    return result

def process_header(sample_filename):
    print "Processing header for " + sample_filename
    f = gzip.open(sample_filename, 'r')
    h = file.open(raw_dir + header_file_name, 'w')
    current_sample_uuid = ''
    for line in f:
        line_to_write = line
        if re.match('##reference=', line):
            line_to_write = '##reference=<ID=hs37d5>\n'
        if re.match('##SAMPLE=', line):
            res = re.search('[0-9a-f]{12}4[0-9a-f]{3}[89ab][0-9a-f]{15}\Z',line)
            if res != None:
                curret_sample_uuid = res.group(0)
                line_to_write = '##SAMPLE=<ID=' + current_sample_uuid + '>'
        h.write(line_to_write)
        
        if re.match('#CHROM', line):
            break
        
    return current_sample_uuid
    
 
def reheader_sample(sample_filename, sample_uuid):
    print "Reheadering " + sample_filename
    call("bcftools reheader -h " + raw_dir + header_file_name + " -o " + raw_dir + sample_uuid + ".vcf.gz " + sample_filename, shell=True)
    return raw_dir + sample_uuid + ".vcf.gz"   

def regenerate_tabix(sample_filename):
    print "Generating index for " + sample_filename
    os.chdir(raw_dir)
    call("tabix -f -p vcf " + sample_filename, shell=True)
    
def get_snp6_to_rtg_mappings(mapping_file_path):
    snp6_to_rtg = {}
    rtg_to_snp6 = {}
    f = file.open(mapping_file_path, "r")
    for line in f.readline():
        rec = line.split(",")
        snp6_to_rtg[rec[0]] = rec[1]
        rtg_to_snp6[rec[1]] = rec[0]
    return snp6_to_rtg, rtg_to_snp6

def get_results_file_name(snp6_uuid, rtg_uuid):
    return snp6_uuid + "_" + rtg_uuid + "_snp6_rtg_concordance_test.txt"

def calculate_genotype_concordance(snp6_uuid, rtg_uuid, rtg_file_path):
    print "Calculating genotype concordance between SNP6 sample - " + snp6_uuid + " and RTG sample - " + rtg_uuid 
    snp6_file_path = snp6_vcf_dir + snp6_uuid + snp6_vcf_suffix
    results_file_name = raw_dir + get_results_file_name(snp6_uuid, rtg_uuid) 
    os.chdir(raw_dir)
    call("vcf-compare -a -g -m " + snp6_uuid + ":" + rtg_uuid + snp6_file_path + " " + rtg_file_path + "  > " + results_file_name, shell=True)
    return results_file_name

def clean_up(results_file_name):
    print "Cleaning up for " + results_file_name
    os.chdir(raw_dir)
    call("mv *.vcf.gz* " + processed_dir, shell=True)
    call("mv " + results_file_name + " " + results_dir, shell=True)
    call("rm " + header_file_name, shell=True)

def snp6_to_rtg_concordance_test():
    result_list = find_sample_files('*.vcf.gz', rtg_vcf_dir)
    snp6_to_rtg, rtg_to_snp6 = get_snp6_to_rtg_mappings(snp6_to_rtg_mapping_file)
    
    for sample_file in result_list:
        snp6_uuid = process_header(sample_file)
        reheadered_sample_filename = reheader_sample(sample_file, snp6_uuid)
        regenerate_tabix(reheadered_sample_filename)
        results_file_name = calculate_genotype_concordance(snp6_uuid, snp6_to_rtg[snp6_uuid], reheadered_sample_filename)

snp6_to_rtg_concordance_test()    
    
    