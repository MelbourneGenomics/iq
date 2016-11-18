#!/usr/bin/env python
'''
####################################################################################
#
# Melbourne Genomics Master Exon Generation Script
#
# Copyright Melbourne Genomics Health Alliance members. All rights reserved.
#
####################################################################################
#
# Purpose:
#     Generates exons.bed, the master bed file used for gene lists
#
# Usage:
#     python combine_overlap_stats.py --exons exon_file --cds cds_file  > combined.tsv
#
####################################################################################
'''

import sys

def combine(exons, cds, target, log):
    '''
        combine stats from exons and cs files and write to target
    '''
    target.write('Gene\tCapture CDS Bases\tCapture Exon Bases\tRefSeq CDS Bases\tRefSeq Exon Bases\tCDS Coverage\tExon Coverage\tAlternative Names\tOriginal Gene\n')
    cds_set = {}
    for line in cds:
        cds_set[line.split('\t')[0].upper()] = line

    stats = {'max': 0, 'min': 0, 'other': 0, 'notfound': 0}
    for line in exons:
        exon_fields = line.strip('\n').split('\t')
        gene_original = exon_fields[0]
        gene = gene_original.upper()
        if gene == 'GENE':
            continue
        if gene in cds_set:
            cds_fields = cds_set[gene].strip('\n').split('\t')
            target.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\n'.format(gene_original, cds_fields[1], exon_fields[1], cds_fields[2], exon_fields[2], cds_fields[3], exon_fields[3], cds_fields[4], cds_fields[5]))
            if cds_fields[3] == '100.00':
                stats['max'] += 1
            elif cds_fields[3] == '0.00':
                stats['min'] += 1
            else:
                stats['other'] += 1
        else:
            target.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\n'.format(gene_original, 0, exon_fields[1], 0, exon_fields[2], 0, exon_fields[3], exon_fields[4], exon_fields[5]))
            stats['notfound'] += 1

    log.write('{0}\n'.format(stats))

def main():
    '''
        combine stats using command line
    '''
    import argparse
    parser = argparse.ArgumentParser(description='Compare BAMs')
    parser.add_argument('--exons', required=True, help='exon overlap file')
    parser.add_argument('--cds', required=True, help='CDS overlap file')
    args = parser.parse_args()

    combine(open(args.exons, 'r'), open(args.cds, 'r'), target=sys.stdout, log=sys.stderr)

if __name__ == '__main__':
    main()
