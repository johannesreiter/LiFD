import unittest

from varcode import Variant
from pyensembl import ensembl_grch37
from varcode.effects.effect_classes import *
from varcode.effects.effect_ordering import *

class VarcodeTest(unittest.TestCase):

    def setUp(self):
        return

    def tearDown(self):
        return
        
    def test_Varcode(self):
        variants = ( # chr, start_pos, reference allele, alternate allele, worst mutation effect
            (17, 7573996, 'A', 'G', 'Substitution'), 
            (2, 198283615, 'C', 'G', 'IntronicSpliceSite'), 
            (19, 47503648, 'G', 'A', 'PrematureStop'), 
            (14, 69256615, 'CGGTGGCAGCGG', '', 'Deletion'), 
            (5, 112175217, 'A', '', 'FrameShift')
        )
        
        for var in variants:
            var_poss = Variant(contig=var[0], start=var[1], ref=var[2], alt=var[3], ensembl=ensembl_grch37)
            self.assertEqual(var_poss.effects().top_priority_effect().__class__.__name__, var[4])