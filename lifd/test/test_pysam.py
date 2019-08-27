# Requires CGI database in order to run the test
import unittest
import os
import pandas as pd

import pysam
from lifd.settings import DB_DIR

class PysamTest(unittest.TestCase):

    def setUp(self):
        return

    def tearDown(self):
        return
        
    def test_FastaFileAvailability(self):
        self.assertTrue(os.path.isfile(os.path.join(DB_DIR, 'hg19.fa')))
    
    def test_Pysam(self):
        genome = pysam.Fastafile(os.path.join(DB_DIR, 'hg19.fa'))
        tmp = pd.read_csv(os.path.join(DB_DIR, 'catalog_of_validated_oncogenic_mutations.tsv'), delimiter='\t', encoding='Latin-1')
        for entries in tmp['gdna']:
            entry = entries.split('__')[0]
            if 'dup' not in entry:
                continue
                
            chrom, gen = entry.split(':')

            for i, char in enumerate(gen[2:]):
                if char.isalpha():
                    index = i
                    break
            gloc = None
            
            start = None
            end = None
            if '_' in gen:
                gcoords = gen[2:i+2].split('_')
                start = gcoords[0]
                end = gcoords[1]
            else:
                start = gen[2:i+2]
                end = gen[2:i+2]
            start = int(start)
            end = int(end)
            self.assertEqual(genome.fetch(chrom, start-1, end), gen[i+5:])
            self.assertTrue('n' not in genome.fetch(chrom, start-1, end))
            self.assertTrue('N' not in genome.fetch(chrom, start-1, end))