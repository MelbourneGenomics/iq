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
    target.write('Gene\tNextera CDS Bases\tNextera Exon Bases\tRefSeq CDS Bases\tRefSeq Exon Bases\tCDS Coverage\tExon Coverage\tAlternative Names\tExon Mean of Mean Coverage\tExon Mean Coverage SD\tExon Mean of Median Coverage\tExon Mean of Percent>20x\tExon Mean Percent SD\tCDS Mean of Mean Coverage\tCDS Mean Coverage SD\tCDS Mean of Median Coverage\tCDS Mean of Percent>20x\tCDS Mean Percent SD\n')
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
        target.write(line.strip('\n'))

        if gene in exon_info:
            target.write('\t')
            target.write('\t'.join(exon_info[gene].strip('\n').split('\t')[1:6]))
        else:
            target.write('\t')
            target.write('0\t0\t0\t0\t0')

        if gene in cds_info:
            target.write('\t')
            target.write('\t'.join(cds_info[gene].strip('\n').split('\t')[1:6]))
            target.write('\n')
        else:
            target.write('\t')
            target.write('0\t0\t0\t0\t0\n')

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

