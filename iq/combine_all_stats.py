#!/usr/bin/env python
'''
####################################################################################
#
# Copyright Melbourne Genomics Health Alliance members. All rights reserved.
#
####################################################################################
#
# Purpose:
#     combines previously calculated overlap and coverage stats and combines into a single file
#
# Usage:
#     python combine_all_stats.py --overlap combined-overlap.txt --cds final_stats_cds.txt --exons final_stats_raw.txt
#
####################################################################################
'''

import sys

def combine(capture, exon_coverage, cds_coverage, target):
    '''
        combine capture stats with coverage stats
    '''
    target.write('original_gene\thgnc\tnextera_cds_bases\tnextera_exon_bases\trefseq_cds_bases\trefseq_exon_bases\tcds_coverage\texon_coverage\talternative_names\texon_mean_of_mean_coverage\texon_mean_coverage_sd\texon_mean_of_median_coverage\texon_mean_of_percent_20x\texon_mean_percent_sd\tcds_mean_of_mean_coverage\tcds_mean_coverage_sd\tcds_mean_of_median_coverage\tcds_mean_of_percent_20x\tcds_mean_percent_sd\tvep_annotation\n')
    exon_info = {}
    for line in exon_coverage:
        gene = line.split('\t')[0]
        exon_info[gene] = line

    cds_info = {}
    for line in cds_coverage:
        gene = line.split('\t')[0]
        cds_info[gene] = line

    for line in capture:
        gene = line.split('\t')[0]
        if gene == 'Gene':
            continue
        capture_fields = line.strip('\n').split('\t')
        target.write('\t'.join([capture_fields[-1]] + capture_fields[:-1])) # move original gene name to start

        if gene in exon_info:
            target.write('\t')
            target.write('\t'.join(exon_info[gene].strip('\n').split('\t')[1:6]))
        else:
            target.write('\t')
            target.write('0\t0\t0\t0\t0')

        if gene in cds_info:
            target.write('\t')
            target.write('\t'.join(cds_info[gene].strip('\n').split('\t')[1:6]))
        else:
            target.write('\t')
            target.write('0\t0\t0\t0\t0')

        target.write('\t1') # vep_annotation
        target.write('\n')

def main():
    '''
        execute via command line
    '''
    import argparse
    parser = argparse.ArgumentParser(description='Combine all stats')
    parser.add_argument('--capture_overlap', help='stats on capture overlap')
    parser.add_argument('--exons', help='exon coverage stats')
    parser.add_argument('--cds', help='cds coverage stats')
    args = parser.parse_args()

    combine(open(args.capture_overlap, 'r'), open(args.exons, 'r'), open(args.cds, 'r'), sys.stdout)

if __name__ == '__main__':
    main()

