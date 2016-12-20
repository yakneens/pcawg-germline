#!/usr/bin/env python2
"""

"""

import sys, argparse, os, os.path
import pandas as pd

PCAWG_METAINFO = os.path.abspath(os.path.dirname(sys.argv[0]) + '/..') + '/configuration/release_may2016.v1.4.tsv'
GNOS_BAMS = '/gnosdata/data'
GNOS_TCGA_BAMS = '/gnosdata/tcga'

def locate(bamfn):
    bam_path = GNOS_BAMS + '/' + bamfn
    if os.path.exists(bam_path) == False:
        bam_path = GNOS_TCGA_BAMS + '/' + bamfn
        if os.path.exists(bam_path) == False:
            return None
    return bam_path

def pair_tumor_normal(pcawg, of):
    p = pcawg[['normal_wgs_aliquot_id',
               'normal_wgs_bwa_alignment_gnos_id',
               'normal_wgs_bwa_alignment_bam_file_name',
               'tumor_wgs_aliquot_id',
               'tumor_wgs_bwa_alignment_gnos_id',
               'tumor_wgs_bwa_alignment_bam_file_name']]
    of.write('normal_wgs_aliquot_id\tStatus\tnormal_bam_path\ttumor_wgs_aliquot_id\ttumor_bam_path\n')
    for i, row in p.iterrows():
        n_id = row['normal_wgs_aliquot_id'].split(',')
        if len(n_id) > 1:
            sys.stderr.write('Warning: multiple normal IDs - skipping: %s\n' % (n_id))
            for nid in n_id:
                of.write('%s\tMultipleNormalIDs\tNone\tNone\tNone\n' % (nid))
            continue

        n_id = n_id[0]
        n_path = row['normal_wgs_bwa_alignment_gnos_id']
        n_bam = row['normal_wgs_bwa_alignment_bam_file_name']
        nbam_path = locate(n_path + '/' + n_bam)
        if nbam_path is None:
            sys.stderr.write('Warning: unable to locate normal bam - skipping: %s\n' % (n_path))
            of.write('%s\tNormalBamNotFound\tNone\tNone\tNone\n' % (n_id))
            continue

        t_id = row['tumor_wgs_aliquot_id'].split(',')
        t_path = row['tumor_wgs_bwa_alignment_gnos_id'].split(',')
        t_bam = row['tumor_wgs_bwa_alignment_bam_file_name'].split(',')
        for tix, tid in enumerate(t_id):
            tbam_path = locate(t_path[tix] + '/' + t_bam[tix])
            if tbam_path is None:
                sys.stderr.write('Warning: unable to locate tumor bam - skipping: %s\n' % (t_path[ix]))
                of.write('%s\tTumorBamNotFound\t%s\t%s\tNone\n' % (n_id, nbam_path, tid))
                continue
            of.write('%s\n' % ('\t'.join([n_id, 'Success', nbam_path, tid, tbam_path])))
    #print pcawg.columns

if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('-i', '--input', help = 'PCAWG metainfo tsv file [%(default)s]', default = PCAWG_METAINFO)
    p.add_argument('-o', '--output', help = 'A list of tumor-normal pairs to output')
    args = p.parse_args()
    pcawg = pd.read_csv(args.input, sep = '\t')
    pair_tumor_normal(pcawg, open(args.output, 'w'))
