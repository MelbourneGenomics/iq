#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --account=VR0320
#SBATCH --time=12:00:00
#SBATCH --mem-per-cpu=16784

# this is an example script building exon lists and finally a tsv file of all statistics
# assumes the existence of 6 bam files in the data directory; replace with your own list
# assumes the existence of target_regions.bed; see the README for instructions on obtaining these

VER=160302

# build exons
#python iq/build_exon_bed.py > data/exons-ref.$VER
#python iq/build_exon_bed.py --cds --refseq mysql717 > data/exons-cds.$VER

# get latest data
#wget ftp://ftp.ebi.ac.uk/pub/databases/genenames/new/tsv/hgnc_complete_set.txt
#mv hgnc_complete_set.txt data/hgnc_complete_set.$VER.txt 

# measure overlap of genes with capture 
python iq/gene_intersect.py --capture ./data/target_regions.bed --hgnc ./data/hgnc_complete_set.$VER.txt < ./data/exons-cds.$VER > ./data/coverage-cds.$VER.txt
python iq/gene_intersect.py --capture ./data/target_regions.bed --hgnc ./data/hgnc_complete_set.$VER.txt < ./data/exons-raw.$VER > ./data/coverage-raw.$VER.txt
python iq/combine_overlap_stats.py --exons ./data/coverage-raw.$VER.txt --cds ./data/coverage-cds.$VER.txt > ./data/combined-coverage.$VER.txt

# cds 
python iq/add_bam_to_iq.py --bam ./data/sample1.bam --exon ./data/exons-cds.$VER > 005.out
python iq/add_bam_to_iq.py --bam ./data/sample2.bam --exon ./data/exons-cds.$VER > 007.out
python iq/add_bam_to_iq.py --bam ./data/sample3.bam --exon ./data/exons-cds.$VER > 009.out
python iq/add_bam_to_iq.py --bam ./data/sample4.bam --exon ./data/exons-cds.$VER > 010.out
python iq/add_bam_to_iq.py --bam ./data/sample5.bam --exon ./data/exons-cds.$VER > 012.out
python iq/add_bam_to_iq.py --bam ./data/sample6.bam --exon ./data/exons-cds.$VER > 014.out

# now merge
python iq/merge_iq.py 005.out 007.out 009.out 010.out 012.out 014.out > ./data/final_stats_cds.$VER.txt

# raw 
python iq/add_bam_to_iq.py --bam ./data/sample1.bam --exon ./data/exons-raw.$VER > 005.out
python iq/add_bam_to_iq.py --bam ./data/sample2.bam --exon ./data/exons-raw.$VER > 007.out
python iq/add_bam_to_iq.py --bam ./data/sample3.bam --exon ./data/exons-raw.$VER > 009.out
python iq/add_bam_to_iq.py --bam ./data/sample4.bam --exon ./data/exons-raw.$VER > 010.out
python iq/add_bam_to_iq.py --bam ./data/sample5.bam --exon ./data/exons-raw.$VER > 012.out
python iq/add_bam_to_iq.py --bam ./data/sample6.bam --exon ./data/exons-raw.$VER > 014.out

# now merge
python iq/merge_iq.py 005.out 007.out 009.out 010.out 012.out 014.out > ./data/final_stats_raw.$VER.txt

# final combo
python iq/combine_all_stats.py --capture_overlap ./data/combined-coverage.$VER.txt --exons ./data/final_stats_raw.$VER.txt --cds ./data/final_stats_cds.$VER.txt > ./data/all_stats.$VER.txt

# now add coords to output
python iq/add_gene_coords.py ./data/exons-raw.$VER ./data/hgnc_complete_set.$VER.txt < ./data/all_stats.$VER.txt > ../data/all_stats_coords.$VER.tsv
