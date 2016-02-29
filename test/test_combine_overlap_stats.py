import StringIO
import unittest

import iq.combine_overlap_stats

class TestCombineOverlapStats(unittest.TestCase):
    def test_simple(self):
        exons = ['A1CF\t1\t2\t50.00\tALT1,ALT2', 'A2M\t3\t4\t75.00\t']
        cds = ['A2M\t5\t6\t83.33\tALT3']
        target = StringIO.StringIO()
        log = StringIO.StringIO()
        iq.combine_overlap_stats.combine(exons, cds, target, log)
        lines = target.getvalue().split('\n')
        assert len(lines) == 4
        assert lines[1] == 'A1CF\t0\t1\t0\t2\t0\t50.00\tALT1,ALT2' # no cds data
        assert lines[2] == 'A2M\t5\t3\t6\t4\t83.33\t75.00\tALT3' # data for both
        assert lines[3] == ''

    def test_case(self):
        exons = ['a1CF\t1\t2\t50.00\tALT1,ALT2', 'A2m\t3\t4\t75.00\t']
        cds = ['A2M\t5\t6\t83.33\tALT3']
        target = StringIO.StringIO()
        log = StringIO.StringIO()
        iq.combine_overlap_stats.combine(exons, cds, target, log)
        lines = target.getvalue().split('\n')
        assert len(lines) == 4
        assert lines[1] == 'a1CF\t0\t1\t0\t2\t0\t50.00\tALT1,ALT2' # no cds data
        assert lines[2] == 'A2m\t5\t3\t6\t4\t83.33\t75.00\tALT3' # data for both
        assert lines[3] == ''

if __name__ == '__main__':
    unittest.main()
