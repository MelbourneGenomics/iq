
Analyzing, measuring, and predicting alignment quality
======================================================

Purpose
-------

Usage
-----

Generate a bed file of exons from UCSC refseq:
```
python build_exon_bed.py > exons-cds.160207
python build_exon_bed.py --cds --refseq mysql717 > exons-raw.160207
rm *717*
```

Calculate the intersection of genes in a bed file with a capture file, along with alternative gene names:
```
python gene_intersect.py --capture target_regions.bed --hgnc hgnc_complete_set.txt < exons-cds.160207 > overlap-cds.160207.txt
python gene_intersect.py --capture target_regions.bed --hgnc hgnc_complete_set.txt < exons-raw.160207 > overlap-raw.160207.txt
```

Obtaining Nextera 1.2 target regions for the previous step:
```
curl http://support.illumina.com/content/dam/illumina-support/documents/documentation/chemistry_documentation/samplepreps_nextera/nexterarapidcapture/nexterarapidcapture_exome_targetedregions_v1.2.bed > target_regions.bed
```

Obtaining HGNC alternate names for the previous step:
```
wget ftp://ftp.ebi.ac.uk/pub/databases/genenames/new/tsv/hgnc_complete_set.txt
```

Combine calculated overlap stats into a single file:
```
python combine_overlap_stats.py --exons overlap-raw.160207.txt --cds overlap-cds.160207.txt > combined-overlap.txt
```

Generate coverage stats using multiple bam inputs:
```
python add_bam_to_iq.py --bam sample1.bam --exon exons-cds.160207 > 1.out
python add_bam_to_iq.py --bam sample2.bam --exon exons-cds.160207 > 2.out
python add_bam_to_iq.py --bam sample3.bam --exon exons-cds.160207 > 3.out

python merge_iq.py 1.out 2.out 3.out > final_stats_cds.txt

python add_bam_to_iq.py --bam sample1.bam --exon exons-raw.160207 > 1.out
python add_bam_to_iq.py --bam sample2.bam --exon exons-raw.160207 > 2.out
python add_bam_to_iq.py --bam sample3.bam --exon exons-raw.160207 > 3.out

python merge_iq.py 1.out 2.out 3.out > final_stats_raw.txt
```

Combine all results:
```
python combine_all_stats.py --capture_overlap combined-overlap.txt --cds final_stats_cds.txt --exons final_stats_raw.txt
```

Running unit tests
------------------
To run the tests:
```
python -m unittest discover
```

Calculating unit test coverage
------------------------------
To calculate code coverage:
```
coverage run -m unittest discover
coverage report -m
```

Code style
----------
We ran:
```
pylint --disable line-too-long *.py
```

TODO
----
* more unit tests

Licence
-------
The MIT License (MIT)

Copyright (c) [2016] [Melbourne Genomics Health Alliance]

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
