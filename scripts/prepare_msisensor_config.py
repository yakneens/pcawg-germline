#!/usr/bin/env python2
"""
Generate json config files for a PCAWG MSIsensor run.

Tumor-normal pairing is read from the PCAWG sample metainfo table.
"""

import sys, os, argparse, json
import pandas as pd
from generate_tumor_normal_pairs import pair_tumor_normal, PCAWG_METAINFO

PAIRED_INFO = os.path.abspath(os.path.dirname(sys.argv[0]) + '/..') + '/configuration/release_may2016.v1.4.tumor_normal_pairs.tsv'

def generate_config(paired_info_fn, out_dir):
    pairs = pd.read_csv(paired_info_fn, sep = '\t',
        names = ['normal_wgs_aliquot_id', 'Status', 'normal_bam_path',
                'tumor_wgs_aliquot_id', 'tumor_bam_path'])
    pairs = pairs.loc[pairs['Status'] == 'Success'][['normal_bam_path', \
        'tumor_wgs_aliquot_id', 'tumor_bam_path']]
    try:
        os.makedirs(out_dir)
    except:
        pass
    for d in pairs.to_dict(orient = 'records'):
        with open("{}/{}.json".format(out_dir, d['tumor_wgs_aliquot_id']), 'w') as o:
            o.write(json.dumps({"sample" : d}))

if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('-i', '--input', help = 'PCAWG metainfo tsv file [%(default)s]', default = PCAWG_METAINFO)
    p.add_argument('-p', '--pairs', help = 'Paired tumor-normal info (intermediate output) [%(default)s]', default = PAIRED_INFO)
    p.add_argument('-o', '--output', help = 'Output directory', required = True)
    args = p.parse_args()
    if os.path.exists(args.pairs) == False:
        print 'Generating tumor-normal pairing...'
        pcawg = pd.read_csv(args.input, sep = '\t')
        pair_tumor_normal(pcawg, open(args.pairs, 'w'))
    else:
        print "Existing tumor-normal pairing found in " + args.pairs
    generate_config(args.pairs, args.output)
