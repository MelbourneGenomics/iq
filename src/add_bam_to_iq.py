
# generates gene coverage stats from a bam file
# usage: add_bam_to_iq bam exon > stats
import collections
import os
import random
import sys

coverage_threshold = 20


def mean( vals ):
  return 1. * sum( vals ) / len( vals )

def median( vals ):
  return sorted( vals )[ len(vals) / 2 ]

def percent_above_threshold( vals, threshold ):
  return 100. * sum( [ x > threshold for x in vals ] ) / len( vals )

def run( cmd ):
  sys.stderr.write( "executing: {0}...\n".format( cmd ) )
  os.system( cmd )
  sys.stderr.write( "executing: done\n" )

def write_gene( name, vals, fh_out, threshold ):
  fh_out.write( '{0}\t{1:.2f}\t{2}\t{3:.2f}\n'.format( name, mean( vals ), median( vals ), percent_above_threshold( vals, threshold ) ) )

def summarize_bed( fh_in, fh_out ):
  fh_out.write( 'Gene\tMean Coverage\tMedian Coverage\tPercent>{0}\n'.format( coverage_threshold ) )
  sys.stderr.write( "generating stats...\n" )
  stats = collections.defaultdict(list)
  for line in fh_in:
    fields = line.strip().split()
    gene = fields[3].strip()
    coverage = int(fields[5].strip())
    stats[gene].append( coverage )
  sys.stderr.write( "sorting and writing stats...\n" )
  for key in sorted(stats):
    write_gene( key, stats[key], fh_out, coverage_threshold )
  sys.stderr.write( "done\n" )

if __name__ == '__main__':
  bam_file = sys.argv[1]
  exon_file = sys.argv[2]
  tmp_file = 'tmp.{0}'.format( random.randint(1, 1000000) )
  run( 'bedtools bamtobed -i {0} | bedtools coverage -a - -b {1} -d > {2}'.format( bam_file, exon_file, tmp_file ) )
  summarize_bed( open( tmp_file, 'r' ), sys.stdout ) 
