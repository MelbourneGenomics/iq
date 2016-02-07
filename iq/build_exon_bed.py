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
#     python build_exon_list.py > final.bed
#
# Dependencies:
#     * bedtools
#     * curl
#     * mysql
#     * sort
#
# Notes
#     * downloads all required files
#     * rm *717* to get rid of tmp files
####################################################################################
'''

import datetime
import os
import sys

def write_log(log, msg):
    '''
        write message to specified log target
    '''
    log.write('{0}\n'.format(msg))

def run(cmd, log):
    '''
        executes a system command
    '''
    write_log(log, 'executing {0}\n'.format(cmd))
    os.system(cmd)

def write_row(target, chrom, start, end, gene):
    '''
        write bed compatible line to target
    '''
    if start != end:
        target.write('{0}\t{1}\t{2}\t{3}\n'.format(chrom, start, end, gene))

def download_refseq(log):
    '''
        download required data from ucsc
        e.g. format:
        chr1    A3GALT2 33772366,33777652,33778101,33778407,33786676, 33773054,33777790,33778191,33778491,33786699,
        e.g. format:
        chr1    11873 12227 DDX11L1_exon_0_0_chr1_11874_f 0 +
    '''
    idx = 717
    # get the data
    run("mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -e 'select chrom, name2, cdsStart, cdsEnd, exonStarts, exonEnds from refGene' hg19 > mysql{0}".format(idx), log)
    return open('mysql{0}'.format(idx), 'r')

def process_exon(tmp_fh, cds, start_in, end_in, cds_start, cds_end, fields, log):
    '''
        process one exon in the context of cds and write to tmp_fh
    '''
    clipped = 0
    clipped_count = 0
    if len(start_in) > 0:
        start = int(start_in)
        end = int(end_in)
        if cds:
            if end < cds_start or start > cds_end: # exon is to the left or to the right of coding -> clip it
                write_log(log, 'clipped entire gene: {0} from {1}'.format(end - start, fields[1]))
                clipped = end - start
                clipped_count = 1
            elif start >= cds_start and end <= cds_end: # exon is contained inside coding -> write entire exon
                write_row(tmp_fh, fields[0], start, end, fields[1])
            elif start <= cds_start and end <= cds_end: # exon overlaps at start but doesn't extend past end -> start from cds
                write_row(tmp_fh, fields[0], cds_start, end, fields[1])
                clipped = cds_start - start
                clipped_count = 1
                write_log(log, 'clipped {0} at start from {1}'.format(cds_start - start, fields[1]))
            elif start >= cds_start and end > cds_end: # exon overlaps at end but not start -> stop at cds
                write_row(tmp_fh, fields[0], start, cds_end, fields[1])
                clipped = end - cds_end
                clipped_count = 1
                write_log(log, 'clipped {0} at end from {1}'.format(end - cds_end, fields[1]))
            elif cds_start >= start and cds_end <= end: # coding inside exon -> write coding
                write_row(tmp_fh, fields[0], cds_start, cds_end, fields[1])
                clipped = end - cds_end + cds_start - start
                clipped_count = 1
                write_log(log, 'clipped {0} at start and end from {1}'.format(end - cds_end + cds_start - start, fields[1]))
            else:
                write_log(log, 'unhandled exon {0}->{1} cds {2}->{3}'.format(start, end, cds_start, cds_end))
        else:
            write_row(tmp_fh, fields[0], start, end, fields[1])
    return (clipped, clipped_count)

def calculate(refseq, cds=False, padding=0, target=sys.stdout, log=sys.stderr):
    '''
        from refseq, generate a bed file of exons and write to stdout
    '''
    idx = 717

    clipped = 0
    clipped_count = 0
    with open('tmp{0}'.format(idx), 'w') as tmp_fh:
        first = True
        for line in refseq:
            if first:
                first = False
                continue
            fields = line.strip().split('\t') # chr, gene, starts, ends
            # remove chromosomes containing _
            if '_' in fields[0]:
                continue
            # remove anything following _ in gene
            if '_' in fields[1]:
                fields[1] = fields[1][:fields[1].find('_')]
            cds_start = int(fields[2])
            cds_end = int(fields[3])
            starts = fields[4].split(',')
            ends = fields[5].split(',')
            for start_in, end_in in zip(starts, ends):
                delta_clipped, delta_count = process_exon(tmp_fh, cds, start_in, end_in, cds_start, cds_end, fields, log)
                clipped += delta_clipped
                clipped_count += delta_count
    write_log(log, 'clipped: {0} clipped_count: {1}'.format(clipped, clipped_count))
    # sort it
    run('sort -k1,1 -k2,2n < tmp{0} | uniq > sorted{0}.bed'.format(idx), log)

    # add slop
    if padding > 0:
        run('curl http://support.illumina.com/content/dam/illumina-support/documents/documentation/chemistry_documentation/samplepreps_nextera/nexterarapidcapture/nexterarapidcapture_exome_targetedregions_v1.2.bed > nextera{0}.bed'.format(idx), log)
        run('mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -e "select chrom, size from hg19.chromInfo" > hg19.{0}.genome'.format(idx), log)

        run('bedtools slop -i nextera{0}.bed -g hg19.{0}.genome -b {1} > padding{0}.bed'.format(idx, padding), log)
        run('bedtools intersect -a sorted{0}.bed -b padding{0}.bed > padded{0}.bed'.format(idx), log)
        run('sort -k1,1 -k2,2n < padded{0}.bed > paddedsorted{0}.bed'.format(idx), log)
        run('bedtools merge -i paddedsorted{0}.bed -nms > merged{0}.bed'.format(idx), log)
    else:
        run('cp sorted{0}.bed merged{0}.bed'.format(idx, idx), log)

    # breaks multiple genes into separate lines
    target.write('#version %s\n' % datetime.datetime.now().strftime("%Y%m%d"))
    for line in open('merged{0}.bed'.format(idx), 'r'):
        if line.startswith('#'):
            target.write(line)
        else:
            fields = line.strip().split('\t')
            genes = set(fields[3].split(';'))
            for gene in genes:
                target.write('%s\t%s\t%s\t%s\n' % (fields[0], fields[1], fields[2], gene))

def main():
    '''
        execute using command line arguments
    '''
    import argparse
    parser = argparse.ArgumentParser(description='Generate exons')
    parser.add_argument('--cds', action='store_true', default=False, help='extract cds regions')
    parser.add_argument('--refseq', required=False, help='provide downloaded refseq file')
    args = parser.parse_args()
    if args.refseq is None:
        refseq = download_refseq(log=sys.stderr)
    else:
        refseq = open(args.refseq, 'r')

    calculate(refseq, cds=args.cds, padding=0, target=sys.stdout, log=sys.stderr)

if __name__ == '__main__':
    main()
