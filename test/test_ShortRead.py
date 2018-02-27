import unittest
from ShortRead import ShortRead, IncorrectCigarException, CLIPzCigar

class CLIPzCigarTest(unittest.TestCase):

    def setUp(self):
        # plus strand
        self.cigar1 = CLIPzCigar(cigar="30", strand="+", chrom="chr1", position=1000)
        self.cigar2 = CLIPzCigar(cigar="20MTC5IG1MTG11DA9", strand="+", chrom="chr1", position=1000)
        self.cigar3 = CLIPzCigar(cigar="19MTC4DA15MTG3IT1", strand="+", chrom="chr1", position=1000)
        self.cigar4 = CLIPzCigar(cigar="5DADCDGDTITITIADC10", strand="+", chrom="chr1", position=1000)
        self.cigar5 = CLIPzCigar(cigar="MTCMTCITITMTCDADAITDA10", strand="+", chrom="chr1", position=1000)

        # minus strand
        self.cigar1_minus = CLIPzCigar(cigar="30", strand="-", chrom="chr1", position=1000)
        self.cigar2_minus = CLIPzCigar(cigar="20MTC5IG1MTG11DA9", strand="-", chrom="chr1", position=1000)
        self.cigar3_minus = CLIPzCigar(cigar="19MTC4DA15MTG3IT1", strand="-", chrom="chr1", position=1000)
        self.cigar4_minus = CLIPzCigar(cigar="5DADCDGDTITITIADC10", strand="-", chrom="chr1", position=1000)
        self.cigar5_minus = CLIPzCigar(cigar="MTCMTCITITMTCDADAITDA10", strand="-", chrom="chr1", position=1000)
        self.cigar6_minus = CLIPzCigar(cigar="MCTMCT10MAGMAG", strand="-", chrom="chr1", position=1000)

    def test_cigar_properly_split_the_string(self):
        self.assertEqual(self.cigar1.cigar_list, ['30'])
        self.assertEqual(self.cigar2.cigar_list, ['20', 'MTC', '5', 'IG', '1', 'MTG', '11', 'DA', '9'])
        self.assertEqual(self.cigar3.cigar_list, ['19', 'MTC', '4', 'DA', '15', 'MTG', '3', 'IT', '1'])
        self.assertEqual(self.cigar4.cigar_list, ['5', 'DA', 'DC', 'DG', 'DT', 'IT', 'IT', 'IA','DC', '10'])
        self.assertEqual(self.cigar5.cigar_list, ["MTC", "MTC", "IT", "IT", "MTC",'DA', 'DA', 'IT', 'DA', "10"])

        # the same should be for minus strand
        self.assertEqual(self.cigar1_minus.cigar_list, ['30'])
        self.assertEqual(self.cigar2_minus.cigar_list, ['20', 'MTC', '5', 'IG', '1', 'MTG', '11', 'DA', '9'])
        self.assertEqual(self.cigar3_minus.cigar_list, ['19', 'MTC', '4', 'DA', '15', 'MTG', '3', 'IT', '1'])
        self.assertEqual(self.cigar4_minus.cigar_list, ['5', 'DA', 'DC', 'DG', 'DT', 'IT', 'IT', 'IA','DC', '10'])
        self.assertEqual(self.cigar5_minus.cigar_list, ["MTC", "MTC", "IT", "IT", "MTC",'DA', 'DA', 'IT', 'DA', "10"])
        self.assertEqual(self.cigar6_minus.cigar_list, ["MCT", "MCT", "10", "MAG", "MAG"])

    def test_cigar_parsed_properly_for_plus_strand(self):
        # cigar2
        self.assertEqual(len(self.cigar2.features['MTC']), 1)
        self.assertEqual(len(self.cigar2.features['IG']), 1)
        self.assertEqual(len(self.cigar2.features['MTG']), 1)
        self.assertEqual(len(self.cigar2.features['DA']), 1)

        # cigar3
        self.assertEqual(len(self.cigar3.features['MTC']), 1)
        self.assertEqual(len(self.cigar3.features['DA']), 1)
        self.assertEqual(len(self.cigar3.features['MTG']), 1)
        self.assertEqual(len(self.cigar3.features['IT']), 1)

        # cigar4
        self.assertEqual(len(self.cigar4.features['DA']), 1)
        self.assertEqual(len(self.cigar4.features['DC']), 2)
        self.assertEqual(len(self.cigar4.features['DG']), 1)
        self.assertEqual(len(self.cigar4.features['DT']), 1)
        self.assertEqual(len(self.cigar4.features['IT']), 2)
        self.assertEqual(len(self.cigar4.features['IA']), 1)

        # cigar5
        self.assertEqual(len(self.cigar5.features['MTC']), 3)
        self.assertEqual(len(self.cigar5.features['IT']), 3)
        self.assertEqual(len(self.cigar5.features['DA']), 3)

    def test_cigar_parsed_properly_for_minus_strand(self):
        # cigar2_minus
        self.assertEqual(len(self.cigar2_minus.features['MAG']), 1)
        self.assertEqual(len(self.cigar2_minus.features['IC']), 1)
        self.assertEqual(len(self.cigar2_minus.features['MAG']), 1)
        self.assertEqual(len(self.cigar2_minus.features['DT']), 1)

        # cigar3_minus
        self.assertEqual(len(self.cigar3_minus.features['MAG']), 1)
        self.assertEqual(len(self.cigar3_minus.features['DT']), 1)
        self.assertEqual(len(self.cigar3_minus.features['MAC']), 1)
        self.assertEqual(len(self.cigar3_minus.features['IA']), 1)

        # cigar4_minus
        self.assertEqual(len(self.cigar4_minus.features['DT']), 1)
        self.assertEqual(len(self.cigar4_minus.features['DG']), 2)
        self.assertEqual(len(self.cigar4_minus.features['DC']), 1)
        self.assertEqual(len(self.cigar4_minus.features['DA']), 1)
        self.assertEqual(len(self.cigar4_minus.features['IA']), 2)
        self.assertEqual(len(self.cigar4_minus.features['IT']), 1)

        # cigar5_minus
        self.assertEqual(len(self.cigar5_minus.features['MAG']), 3)
        self.assertEqual(len(self.cigar5_minus.features['IA']), 3)
        self.assertEqual(len(self.cigar5_minus.features['DT']), 3)

        # cigar6_minus
        self.assertEqual(len(self.cigar6_minus.features['MTC']), 2)
        self.assertEqual(len(self.cigar6_minus.features['MGA']), 2)


    def test_features_on_plus_strand_have_correct_position(self):
        # cigar2
        self.assertEqual(self.cigar2.features['MTC'][0].iv.start, 1020)
        self.assertEqual(self.cigar2.features['IG'][0].iv.start, 1025)
        self.assertEqual(self.cigar2.features['MTG'][0].iv.start, 1027)
        self.assertEqual(self.cigar2.features['DA'][0].iv.start, 1039)

        # cigar3
        self.assertEqual(self.cigar3.features['MTC'][0].iv.start, 1019)
        self.assertEqual(self.cigar3.features['DA'][0].iv.start, 1024)
        self.assertEqual(self.cigar3.features['MTG'][0].iv.start, 1040)
        self.assertEqual(self.cigar3.features['IT'][0].iv.start, 1043)

        # cigar4
        self.assertEqual(self.cigar4.features['DA'][0].iv.start, 1005)
        self.assertEqual(self.cigar4.features['DC'][0].iv.start, 1006)
        self.assertEqual(self.cigar4.features['DG'][0].iv.start, 1007)
        self.assertEqual(self.cigar4.features['DT'][0].iv.start, 1008)
        self.assertEqual(self.cigar4.features['IT'][0].iv.start, 1008)
        self.assertEqual(self.cigar4.features['IT'][1].iv.start, 1008)
        self.assertEqual(self.cigar4.features['IA'][0].iv.start, 1008)
        self.assertEqual(self.cigar4.features['DC'][1].iv.start, 1009)

        # cigar5
        self.assertEqual(self.cigar5.features['MTC'][0].iv.start, 1000)
        self.assertEqual(self.cigar5.features['MTC'][1].iv.start, 1001)
        self.assertEqual(self.cigar5.features['IT'][0].iv.start, 1001)
        self.assertEqual(self.cigar5.features['IT'][1].iv.start, 1001)
        self.assertEqual(self.cigar5.features['MTC'][2].iv.start, 1002)
        self.assertEqual(self.cigar5.features['DA'][0].iv.start, 1003)
        self.assertEqual(self.cigar5.features['DA'][1].iv.start, 1004)
        self.assertEqual(self.cigar5.features['IT'][2].iv.start, 1004)
        self.assertEqual(self.cigar5.features['DA'][2].iv.start, 1005)

    def test_features_on_minus_strand_have_correct_position(self):
        # cigar2_minus
        self.assertEqual(self.cigar2_minus.features['MAG'][0].iv.start, 1020)
        self.assertEqual(self.cigar2_minus.features['IC'][0].iv.start, 1025)
        self.assertEqual(self.cigar2_minus.features['MAC'][0].iv.start, 1027)
        self.assertEqual(self.cigar2_minus.features['DT'][0].iv.start, 1039)

        # cigar3_minus
        self.assertEqual(self.cigar3_minus.features['MAG'][0].iv.start, 1019)
        self.assertEqual(self.cigar3_minus.features['DT'][0].iv.start, 1024)
        self.assertEqual(self.cigar3_minus.features['MAC'][0].iv.start, 1040)
        self.assertEqual(self.cigar3_minus.features['IA'][0].iv.start, 1043)

        # cigar4_minus
        self.assertEqual(self.cigar4_minus.features['DT'][0].iv.start, 1005)
        self.assertEqual(self.cigar4_minus.features['DG'][0].iv.start, 1006)
        self.assertEqual(self.cigar4_minus.features['DC'][0].iv.start, 1007)
        self.assertEqual(self.cigar4_minus.features['DA'][0].iv.start, 1008)
        self.assertEqual(self.cigar4_minus.features['IA'][0].iv.start, 1008)
        self.assertEqual(self.cigar4_minus.features['IA'][1].iv.start, 1008)
        self.assertEqual(self.cigar4_minus.features['IT'][0].iv.start, 1008)
        self.assertEqual(self.cigar4_minus.features['DG'][1].iv.start, 1009)

        # cigar5_minus
        self.assertEqual(self.cigar5_minus.features['MAG'][0].iv.start, 1000)
        self.assertEqual(self.cigar5_minus.features['MAG'][1].iv.start, 1001)
        self.assertEqual(self.cigar5_minus.features['IA'][0].iv.start, 1001)
        self.assertEqual(self.cigar5_minus.features['IA'][1].iv.start, 1001)
        self.assertEqual(self.cigar5_minus.features['MAG'][2].iv.start, 1002)
        self.assertEqual(self.cigar5_minus.features['DT'][0].iv.start, 1003)
        self.assertEqual(self.cigar5_minus.features['DT'][1].iv.start, 1004)
        self.assertEqual(self.cigar5_minus.features['IA'][2].iv.start, 1004)
        self.assertEqual(self.cigar5_minus.features['DT'][2].iv.start, 1005)

        # cigar6_minus
        self.assertEqual(self.cigar6_minus.features['MTC'][0].iv.start, 1012)
        self.assertEqual(self.cigar6_minus.features['MTC'][1].iv.start, 1013)
        self.assertEqual(self.cigar6_minus.features['MGA'][0].iv.start, 1000)
        self.assertEqual(self.cigar6_minus.features['MGA'][1].iv.start, 1001)

    def test_cigar_length_is_correct(self):
        self.assertEqual(self.cigar1.length, 30)
        self.assertEqual(self.cigar2.length, 49)
        self.assertEqual(self.cigar3.length, 45)
        self.assertEqual(self.cigar4.length, 20)
        self.assertEqual(self.cigar5.length, 16)

        # minus should be the same
        self.assertEqual(self.cigar1_minus.length, 30)
        self.assertEqual(self.cigar2_minus.length, 49)
        self.assertEqual(self.cigar3_minus.length, 45)
        self.assertEqual(self.cigar4_minus.length, 20)
        self.assertEqual(self.cigar5_minus.length, 16)
        self.assertEqual(self.cigar6_minus.length, 14)


class ShortReadTest(unittest.TestCase):

    def setUp(self):
        # plus strand
        self.sr1 = ShortRead(chrom="chr1", start=1000, end=1030, strand="+", clipz_cigar="30")
        self.sr2 = ShortRead(chrom="chr5", start=1000, end=1049, strand="+", clipz_cigar="20MTC5IG1MTG11DA9")
        self.sr3 = ShortRead(chrom="chr1", start=1000, end=1045, strand="+", clipz_cigar="19MTC4DA15MTG3IT1")

        # minus strand
        self.sr1_minus = ShortRead(chrom="chr1", start=1000, end=1030, strand="-", clipz_cigar="30")
        self.sr4_minus = ShortRead(chrom="chr22", start=1000, end=1049, strand="-", clipz_cigar="3MTG2DT4IC1MTC36")

    def test_raises_error_when_incorrect_cigar(self):
        with self.assertRaises(IncorrectCigarException):
            ShortRead(chrom="chr22", start=1000, end=1048, strand="-", clipz_cigar="3MTG2DT4IT")

    def test_returns_correct_truncation_position(self):
        self.assertEqual(self.sr1.get_truncation_position(), 1000)
        self.assertEqual(self.sr2.get_truncation_position(), 1000)
        self.assertEqual(self.sr3.get_truncation_position(), 1000)
        self.assertEqual(self.sr1_minus.get_truncation_position(), 1029)
        self.assertEqual(self.sr4_minus.get_truncation_position(), 1048)





if __name__ == '__main__':
        unittest.main()
