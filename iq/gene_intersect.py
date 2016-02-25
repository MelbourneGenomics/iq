#!/usr/bin/env python
'''
####################################################################################
#
# Copyright Melbourne Genomics Health Alliance members. All rights reserved.
#
####################################################################################
#
# Purpose:
#     takes genes from one set, looks for hgnc matches, and calculates overlap
#
# Usage:
#     python gene_intersect.py --capture capture_file --hgnc hgnc_names < exons.bed > stats.tsv
#
####################################################################################
'''

import collections
import datetime
import sys

def write_log(log, msg):
    '''
        write date stamped log message
    '''
    now = datetime.datetime.now().strftime('%y%m%d-%H%M%S')
    log.write('{0}: {1}\n'.format(now, msg))

def find_intersect(refseq, capture_handle, hgnc_handle, target, log):
    '''
        measure overlapping bases between refseq and capture
        also find alternative hgnc gene names
    '''
    write_log(log, 'Parsing capture...')
    # capture
    capture = collections.defaultdict(set)
    for idx, line in enumerate(capture_handle):
        if line.startswith('#'):
            continue
        fields = line.strip().split()
        if len(fields) > 2:
            chrom = fields[0].upper()
            for pos in xrange(int(fields[1]), int(fields[2])):
                capture[chrom].add(pos)
        if idx % 100000 == 0:
            write_log(log, '{0} lines processed...'.format(idx))
    write_log(log, 'Capture: done: {0}'.format(' '.join(capture.keys())))

    # refseq
    write_log(log, 'Processing refseq...')
    refseq_total = collections.defaultdict(set)
    refseq_match = collections.defaultdict(set)
    for idx, line in enumerate(refseq):
        if line.startswith('#'):
            continue
        fields = line.strip().split()
        if len(fields) > 3:
            chrom = fields[0].upper()
            gene = fields[3].upper()
            start = int(fields[1])
            finish = int(fields[2])
            if gene not in refseq_total:
                refseq_total[gene] = set() # so that zero length genes are included
            for pos in xrange(start, finish):
                refseq_total[gene].add(pos)
                if pos in capture[chrom]:
                    refseq_match[gene].add(pos)
            if idx % 100000 == 0:
                write_log(log, 'Checked {0}... last was {1}'.format(idx, chrom))
    write_log(log, 'Processing refseq: done')

    # refseq other names
    # refseq_alt = collections.defaultdict(set)
    write_log(log, 'Parsing hgnc...')
    hgnc = collections.defaultdict(set)
    skipped = 0
    idx = 0
    for idx, line in enumerate(hgnc_handle):
        fields = line.strip().split('\t')
        if len(fields) > 10:
            correct = fields[1].upper()
            old = set()
            alias = fields[8]
            previous = fields[10]
            if len(alias) > 0:
                for i in alias.strip('"').split('|'):
                    old.add(i.upper())
            if len(previous) > 0:
                for i in previous.strip('"').split('|'):
                    old.add(i.upper())
            hgnc[correct] = old
            if idx < 3:
                write_log(log, hgnc)
        else:
            skipped += 1
    write_log(log, 'Parsing hgnc: done. Processed {0} skipped {1}'.format(idx, skipped))
    write_log(log, hgnc['SDHAF3'])

    # write results
    target.write('{0}\t{1}\t{2}\t{3}\t{4}\n'.format('Gene', 'Bases Covered', 'Total Bases', '% Covered', 'Alternate Names'))
    for gene in sorted(refseq_total):
        if len(refseq_total[gene]) == 0:
            target.write('{0}\t{1}\t{2}\t{3:.2f}\t{4}\n'.format(gene, len(refseq_match[gene]), len(refseq_total[gene]), 0, ', '.join(sorted(list(hgnc[gene])))))
        else:
            target.write('{0}\t{1}\t{2}\t{3:.2f}\t{4}\n'.format(gene, len(refseq_match[gene]), len(refseq_total[gene]), 100. * len(refseq_match[gene]) / len(refseq_total[gene]), ', '.join(sorted(list(hgnc[gene])))))

def main():
    '''
        find coverage using command line
    '''
    import argparse
    parser = argparse.ArgumentParser(description='Find intersect stats between capture and stdin')
    parser.add_argument('--capture', required=True, help='capture bed')
    parser.add_argument('--hgnc', required=True, help='hgnc db file')
    args = parser.parse_args()
    find_intersect(sys.stdin, open(args.capture, 'r'), open(args.hgnc, 'r'), sys.stdout, sys.stderr)

if __name__ == '__main__':
    main()

