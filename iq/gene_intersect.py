
'''
    takes genes from one set, looks for hgnc matches, and calculates overlap
'''

import collections
import datetime
import sys

def write_log(log, msg):
    now = datetime.datetime.now().strftime('%y%m%d-%H%M%S')
    log.write('{0}: {1}\n'.format(now, msg))

def find_intersect(refseq, capture_handle, hgnc_handle, target, log):
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
    refseq_total = collections.defaultdict(int)
    refseq_match = collections.defaultdict(int)
    for idx, line in enumerate(refseq):
        if line.startswith( '#' ):
            continue
        fields = line.strip().split()
        if len(fields) > 3:
            chrom = fields[0].upper()
            gene = fields[3].upper()
            start = int(fields[1]) 
            finish = int(fields[2])  
            refseq_total[gene] += finish - start
            refseq_match[gene] += sum([1 if x in capture[chrom] else 0 for x in xrange(start, finish)])
            if idx % 100000 == 0:
                write_log(log, 'Checked {0}... last was {1}'.format(idx, chrom))
    write_log(log, 'Processing refseq: done')

    # refseq other names TODO
    refseq_alt = collections.defaultdict(set)
    write_log(log, 'Parsing hgnc...' )
    hgnc = collections.defaultdict(set)
    skipped = 0
    idx = 0
    for idx, l in enumerate(hgnc_handle):
        f = l.strip().split('\t')
        if len(f) > 10:
            correct = f[1].upper()
            old = set()
            alias = f[8]
            previous = f[10]
            if len(alias) > 0:
                for i in alias.strip('"').split('|'):
                    old.add( i.upper() )
            if len(previous) > 0:
                for i in previous.strip('"').split('|'):
                    old.add( i.upper() )
            hgnc[correct] = old
            if idx < 3:
                write_log(log, hgnc)
        else:
            skipped += 1
    write_log(log,'Parsing hgnc: done. Processed {0} skipped {1}'.format(idx, skipped))
    write_log(log, hgnc['SDHAF3'] )

    # write results
    target.write('{0}\t{1}\t{2}\t{3}\t{4}\n'.format('Gene', 'Bases Covered', 'Total Bases', '% Covered', 'Alternate Names'))
    for gene in sorted(refseq_total):
        target.write('{0}\t{1}\t{2}\t{3:.2f}\t{4}\n'.format(gene, refseq_match[gene], refseq_total[gene], 100. * refseq_match[gene] / refseq_total[gene], ','.join(sorted(list(hgnc[gene])))))

def main():
    '''
        find coverage using command line
    '''
    import argparse
    parser = argparse.ArgumentParser(description='Validate cpipe output')
    parser.add_argument('--capture', required=True, help='capture bed')
    parser.add_argument('--hgnc', required=True, help='hgnc db file')
    args = parser.parse_args()
    find_intersect(sys.stdin, open(args.capture, 'r'), open(args.hgnc, 'r'), sys.stdout, sys.stderr)

if __name__ == '__main__':
    main()

