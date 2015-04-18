import os, fnmatch
from subprocess import call
import gzip, file, re

rootdir = '/icgc/pcawg/analysis/iakhnin/germline_genotype_concordance'
raw_dir = rootdir + '/raw/'
processed_dir = rootdir + '/processed/'
results_dir = rootdir + '/results/'
header_size = 122

def find_sample_files(pattern, path):
    result = []
    for root, dirs, files in os.walk(path):
        for name in files:
            if fnmatch.fnmatch(name, pattern):
                current_sample_filename = os.path.join(root, name)
                result.append(current_sample_filename)
    return result

def process_header(sample_filename):
    f = gzip.open(sample_filename, 'r')
    h = file.open(raw_dir + "header.txt", 'w')
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
    call("bcftools reheader -h " + raw_dir + "header.txt" + " -o " + raw_dir + sample_uuid + ".vcf.gz " + sample_filename, shell=True)
    return raw_dir + sample_uuid + ".vcf.gz"   

def regenerate_tabix(sample_filename):
    os.chdir(raw_dir)
    call("tabix -f -p vcf " + sample_filename, shell=True)

result_list = find_sample_files('*.vcf.gz', '/icgc/pcawg/project/results/germline_variants/annai/tcga')

for sample_file in result_list:
    wgs_uuid = process_header(sample_file)
    reheadered_sample_filename = reheader_sample(sample_file, wgs_uuid)
    regenerate_tabix(reheadered_sample_filename)
    
    
    