
#chr1    161480623       161480746       FCGR2A  1       22
#chr1    161480623       161480746       FCGR2A  2       23

import sys

exon_stats = { 'min': 1e9, 'max': 0 }

def write_coverage( c ):
  sys.stdout.write( '{0}\t{1}\t{2}\t{3}\n'.format( c['chr'], c['gene'], c['start'], ','.join(c['coverage']) ) )
  exon_stats['min'] = min( exon_stats['min'], len(c['coverage']) )
  exon_stats['max'] = max( exon_stats['max'], len(c['coverage']) )

current = None
for i, line in enumerate( sys.stdin ):
  fields = line.strip().split()
  if fields[4] == '1':
    if current is not None:
      write_coverage( current )
    current = { 'chr': fields[0], 'start': fields[1], 'gene': fields[3], 'coverage': [ fields[5] ] }
  else:
    current['coverage'].append( fields[5] )
    # sanity
    expected = current['start']
    if expected != fields[1]:
      sys.stderr.write( 'ERROR on line {0}: expected {1} actual {2}\n'.format( i, expected, fields[1] ) )
  if i % 1000000 == 0:
    sys.stderr.write( 'Processed {0} lines...\n'.format( i ) )

write_coverage( current )

sys.stderr.write( 'Processed {0} lines. Min exon: {1}, Max exon: {2}\n'.format( i, exon_stats['min'], exon_stats['max'] ) )
