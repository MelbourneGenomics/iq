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
#         Generates gene coverage stats from a bam file
#
# Usage:
#         add_bam_to_iq --bam bam --exon exon > stats
#
####################################################################################
'''

import argparse
import collections
import os
import random
import sys

COVERAGE_THRESHOLD = 20

def mean(vals):
    '''
        returns the mean of vals
    '''
    return 1. * sum(vals) / len(vals)

def median(vals):
    '''
        returns the median of vals
    '''
    return sorted(vals)[len(vals) / 2]

def percent_above_threshold(vals, threshold):
    '''
        what percentage of vals are above threshold?
    '''
    return 100. * sum([x > threshold for x in vals]) / len(vals)

def run(cmd, log):
    '''
        execute a system command
    '''
    log.write("executing: {0}...\n".format(cmd))
    os.system(cmd)
    log.write("executing: done\n")

def write_gene(name, vals, fh_out, threshold):
    '''
        write the stats for name to fh_out
    '''
    fh_out.write('{0}\t{1:.2f}\t{2}\t{3:.2f}\n'.format(name, mean(vals), median(vals), percent_above_threshold(vals, threshold)))

def summarize_bed(fh_in, fh_out, log):
    '''
        calculate coverage stats from the incoming bed coverage data
    '''
    fh_out.write('Gene\tMean Coverage\tMedian Coverage\tPercent>{0}\n'.format(COVERAGE_THRESHOLD))
    log.write("generating stats...\n")
    stats = collections.defaultdict(list)
    for line in fh_in:
        fields = line.strip().split()
        gene = fields[3].strip()
        coverage = int(fields[5].strip())
        stats[gene].append(coverage)
    log.write("sorting and writing stats...\n")
    for key in sorted(stats):
        write_gene(key, stats[key], fh_out, COVERAGE_THRESHOLD)
    log.write("done\n")

def main():
    '''
        execute from command line
    '''
    parser = argparse.ArgumentParser(description='Generate coverage data')
    parser.add_argument('--bam', help='bam file to process')
    parser.add_argument('--exon', help='list of exons to filter against')
    args = parser.parse_args()
    bam_file = args.bam
    exon_file = args.exon
    #run('bedtools bamtobed -i {0} | bedtools coverage -a - -b {1} -d'.format(bam_file, exon_file)) # write directly (use for health hack)
    # generate output for merge_iq
    tmp_file = 'tmp.{0}'.format(random.randint(1, 1000000))
    run('bedtools bamtobed -i {0} | bedtools coverage -a - -b {1} -d > {2}'.format(bam_file, exon_file, tmp_file), sys.stderr)
    summarize_bed(open(tmp_file, 'r'), sys.stdout, sys.stderr)

if __name__ == '__main__':
    main()
