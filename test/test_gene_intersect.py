
import StringIO
import unittest

import iq.gene_intersect

class TestGeneIntersect(unittest.TestCase):
    def test_refseq_correct(self):
        refseq = ['chrom\t100\t200\tmygene']
        hgnc_handle = ['0\tcorrected\t2\t3\t4\t5\t6\t7\taliases\t9\tmygene']
        capture_handle = ['chrom\t160\t220\tblah']
        target = StringIO.StringIO()
        log = StringIO.StringIO()
        iq.gene_intersect.find_intersect(refseq, capture_handle, hgnc_handle, target, log)
        fields = target.getvalue().split('\n')[1].split('\t')
        assert fields[0] == 'corrected'

    def test_refseq_no_correct(self):
        refseq = ['chrom\t100\t200\tmygene']
        hgnc_handle = ['0\tcorrected\t2\t3\t4\t5\t6\t7\taliases\t9\t10']
        capture_handle = ['chrom\t160\t220\tblah']
        target = StringIO.StringIO()
        log = StringIO.StringIO()
        iq.gene_intersect.find_intersect(refseq, capture_handle, hgnc_handle, target, log)
        fields = target.getvalue().split('\n')[1].split('\t')
        assert fields[0] == 'mygene'

if __name__ == '__main__':
    unittest.main()
