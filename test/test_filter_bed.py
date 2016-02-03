
import StringIO
import unittest

import iq.filter_bed

class TestFilterBed(unittest.TestCase):
    def test_simple(self):
        src = []
        genes = []
        target = StringIO.StringIO()
        iq.filter_bed.filter_bed(src, target, genes)

if __name__ == '__main__':
    unittest.main()
