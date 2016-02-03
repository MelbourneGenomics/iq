
# remove genes from a bed file
# usage filter_bed.py genes_file < exons.bed > filtered.bed
import sys

def filter_bed(exons, target, genes):
    exclude = set()
    for line in genes:
        exclude.add( line.strip().split()[0].upper() )

    sys.stderr.write( '{0} genes selected for exclusion\n'.format( len( exclude ) ) )

    skipped = 0
    i = 0
    for i, line in enumerate(exons):
        fields = line.strip().split()
        if len(fields) > 3 and fields[3].upper() in exclude:
            skipped += 1
        else:
            target.write(line)

    sys.stderr.write( 'Skipped {0} out of {1}\n'.format( skipped, i ) )

if __name__ == '__main__':
    filter_bed(sys.stdin, sys.stdout, open(sys.argv[1], 'r'))
