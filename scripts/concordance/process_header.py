import gzip, re, os
import sys
rootdir = '/icgc/pcawg/analysis/iakhnin/germline_genotype_concordance'
raw_dir = rootdir + '/raw/'
processed_dir = rootdir + '/processed/'
results_dir = rootdir + '/results/'
header_file_name = "header.txt"

args = sys.argv
sample_filename = args[1]
sample_uuid = args[2]

f = gzip.open(sample_filename, 'r')
os.makedirs(raw_dir + sample_uuid)
h = open(raw_dir + sample_uuid + "/" + header_file_name, 'w')
current_sample_uuid = ''
for line in f:
    line_to_write = line
    if re.match('##reference=', line):
        line_to_write = '##reference=<ID=hs37d5>\n'
    if re.match('##SAMPLE=', line):
        res = re.search('[a-f0-9]{8}-?[a-f0-9]{4}-?4[a-f0-9]{3}-?[89ab][a-f0-9]{3}-?[a-f0-9]{12}',line)
        if res != None:
            current_sample_uuid = res.group(0)
    line_to_write = '##SAMPLE=<ID=' + current_sample_uuid + '>'
    h.write(line_to_write)
    
    if re.match('#CHROM', line):
        break
    
print current_sample_uuid