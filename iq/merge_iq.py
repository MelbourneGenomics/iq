
# merges the output from add_bam_to_iq
# usage: merge_stats bedfile...
import math
import numpy
import sys

def mean( vals ):
  return 1. * sum( vals ) / len( vals )

def sd( vals ):
  return numpy.std( vals )
  #valmean = mean( vals )
  #std = math.sqrt(mean(abs(x - x.mean())**2))

def merge_stats( sources, fh_out ):
  #fh_out.write( '# n={0}\n'.format( len( sources ) ) )
  fh_out.write( '{0}\t{1}\t{4}\t{2}\t{3}\t{5}\n'.format( "Gene", "Mean of Mean Coverage", "Mean of Median Coverage", "Mean of Percent>20x", "Mean Coverage SD", "Mean Percent SD" ) )
  while True:
    lines = [ source.readline() for source in sources ]
    summary = [ [], [], [] ]
    gene = None
    for i, line in enumerate(lines):
      fields = [ x.strip() for x in line.strip().split() ]
      if len(fields) == 0:
        return # done
      if fields[0] != 'Gene':
        if gene is None:
          gene = fields[0]
        elif gene is not None and gene != fields[0]:
          sys.stderr.write( 'Error: gene mismatch on line {0}'.format( i ) )
          return
        summary[0].append( float( fields[1] ) )
        summary[1].append( float( fields[2] ) )
        summary[2].append( float( fields[3] ) )
    if len(summary[0]) > 0:
      fh_out.write( '{0}\t{1:.2f}\t{4:.2f}\t{2:.2f}\t{3:.2f}\t{5:.2f}\n'.format( gene, mean( summary[0] ), mean( summary[1] ), mean( summary[2] ), sd( summary[0] ), sd( summary[2] ) ) )

if __name__ == '__main__':
  merge_stats( [ open(x, 'r') for x in sys.argv[1:] ], sys.stdout )
