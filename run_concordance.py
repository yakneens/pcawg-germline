import os, fnmatch
rootdir = '/icgc/pcawg/project/results/germline_variants/annai/tcga'

def find(pattern, path):
    result = []
    for root, dirs, files in os.walk(path):
        for name in files:
            if fnmatch.fnmatch(name, pattern):
                result.append(os.path.join(root, name))
    return result

find('*.vcf.gz', '/icgc/pcawg/project/results/germline_variants/annai/tcga')