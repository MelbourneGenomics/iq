#!/usr/bin/env python

import collections
import sys

def add_gene_coords(tsv_in, tsv_out, bed, hgnc_handle, log):
    log.write('reading coords...\n')
    chroms = {}
    low = {}
    high = {}
    for idx, line in enumerate(bed):
        fields = line.strip().split('\t')
        if len(fields) < 4:
            continue
        chrom = fields[0]
        start = int(fields[1])
        finish = int(fields[1])
        gene = fields[3]

        chroms[gene] = chrom
        if gene not in low or low[gene] > start:
            low[gene] = start
        if gene not in high or high[gene] < finish:
            high[gene] = finish
        if idx % 100000 == 0:
            log.write('processed {}...\n'.format(idx))

    log.write('Parsing hgnc...\n')
    hgnc_alt = {} # maps old to correct
    hgnc = collections.defaultdict(set)
    skipped = 0
    idx = 0
    for idx, line in enumerate(hgnc_handle):
        fields = line.strip().split('\t')
        if len(fields) > 10:
            original = fields[1]
            correct = original.upper()
            alias = fields[8]
            previous = fields[10]
            if len(alias) > 0:
                for i in alias.strip('"').split('|'):
                    hgnc_alt[i.upper()] = correct
                    hgnc[correct].add(i.upper())
            if len(previous) > 0:
                for i in previous.strip('"').split('|'):
                    hgnc_alt[i.upper()] = correct
                    hgnc[correct].add(i.upper())
        else:
            skipped += 1

    log.write('Parsing hgnc: done. Processed {0} skipped {1}\n'.format(idx, skipped))

    log.write('processing tsv...\n')
    first = True
    for idx, line in enumerate(tsv_in):
        fields = line.strip().split('\t')
        if first:
            new_fields = [fields[0], fields[1], 'chromosome', 'start_pos', 'end_pos'] + fields[2:]
            first = False 
        elif fields[1] in chroms:
            new_fields = [fields[0], fields[1], chroms[fields[1]], str(low[fields[1]]), str(high[fields[1]])] + fields[2:]
        elif fields[1] in hgnc:
            for alt in hgnc[fields[1]]:
                if alt in chroms:
                    new_fields = [fields[0], fields[1], chroms[alt], str(low[alt]), str(high[alt])] + fields[2:]
                    break
        else:
            log.write('WARNING: gene {} not found.\n'.format(fields[1]))
            new_fields = [fields[0], fields[1], '?', '0', '0'] + fields[2:]
        tsv_out.write('\t'.join(new_fields))
        tsv_out.write('\n')
        if idx % 1000 == 0:
            log.write('processed {}...\n'.format(idx))

if __name__ == '__main__':
    add_gene_coords(sys.stdin, sys.stdout, open(sys.argv[1], 'r'), open(sys.argv[2], 'r'), sys.stderr)
