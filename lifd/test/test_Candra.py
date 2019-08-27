
import unittest
import logging
import os
import pandas as pd
import shutil

from lifd.predictors.candra import Candra
from lifd.utils import NT_VAR_COL

from lifd.test.utils_tests import VARIANTS, CANCER_TYPE

__author__ = 'Johannes REITER'
__date__ = 'Jan 19, 2019'

# get logger
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

INPUT_FP = '{}{}'.format(CANCER_TYPE, Candra.INPUT_SUFFIX)
OUTPUT_FP = '{}{}'.format(CANCER_TYPE, Candra.OUTPUT_SUFFIX)


class CandraTest(unittest.TestCase):

    def setUp(self):
        self._delete_tmp_files()

    def tearDown(self):
        self._delete_tmp_files()

    def _delete_tmp_files(self):

        if os.path.exists(INPUT_FP) and os.path.isfile(INPUT_FP):
            os.remove(INPUT_FP)

        if os.path.exists(OUTPUT_FP) and os.path.isfile(OUTPUT_FP):
            os.remove(OUTPUT_FP)

    def test_CandraAvailability(self):

        # is perl installed
        available = shutil.which('perl') is not None
        self.assertTrue(available, 'perl for Candra is not installed')

        # is candra script available?
        script_exists = os.path.exists(os.path.join(Candra.PATH, 'open_candra.pl'))
        self.assertTrue(script_exists,
                        'open_candra.pl is not available at expected path {}. '.format(Candra.PATH)
                        + 'Change location in {}.'.format(os.path.join('src', 'lifd', 'predictors', 'candra.py')))

    def test_Candra(self):

        # mirroring mini pandas dataframe with the essential information of variants
        vars_df = pd.DataFrame(VARIANTS)

        # generate input file
        self.assertFalse(os.path.exists(INPUT_FP), 'Candra input file with variants should not yet exist.')
        Candra.generate_input_file(
            INPUT_FP, vars_df, chromosome_col='Chromosome', position_col='StartPosition',
            reference_col='ReferenceAllele', alternate_col='AlternateAllele')
        self.assertTrue(os.path.exists(INPUT_FP), 'Candra input file with variants should exist.')

        # generate Candra output
        self.assertFalse(os.path.exists(OUTPUT_FP), 'Candra output file should not yet exist.')
        Candra.run(INPUT_FP, OUTPUT_FP, CANCER_TYPE)
        self.assertTrue(os.path.exists(OUTPUT_FP) and os.path.isfile(OUTPUT_FP), 'Candra output file should exist.')

        candra_df = Candra.read_results(OUTPUT_FP)

        self.assertEqual(candra_df.shape, (2, 13),
                         'Dimensions of the dataframe checking: {}'.format(candra_df.shape))

        self.assertEqual(candra_df.loc[VARIANTS[0][NT_VAR_COL]][Candra.SCORE_COL], 10.757)
        self.assertEqual(candra_df.loc[VARIANTS[0][NT_VAR_COL]][Candra.SIGNIFICANCE_COL], 0.085632)
        self.assertEqual(candra_df.loc[VARIANTS[0][NT_VAR_COL]][Candra.CATEGORY_COL], 'Driver')

    def test_Candra_get_results(self):

        # mirroring mini pandas dataframe with the essential information of variants
        var_df = pd.DataFrame(VARIANTS)
        n_cols = var_df.shape[1]

        # generate input file
        self.assertFalse(os.path.exists(INPUT_FP), 'Cravat input file with variants should not yet exist.')
        self.assertFalse(os.path.exists(OUTPUT_FP), 'Cravat output file should not yet exist.')

        # generate CanDrA output
        var_df = Candra.get_annotation(var_df, CANCER_TYPE, input_fp=INPUT_FP, output_fp=OUTPUT_FP)
        self.assertTrue(os.path.exists(INPUT_FP), 'Candra input file with variants should exist.')
        self.assertTrue(os.path.exists(OUTPUT_FP), 'Candra output file should exist.')
        self.assertEqual(var_df.shape, (len(VARIANTS), n_cols + 4),
                         'Dimensions of the dataframe checking: {}'.format(var_df.shape))

        self.assertEqual(var_df[var_df[NT_VAR_COL] == VARIANTS[0][NT_VAR_COL]][
                             Candra.LIFD_SCORE_COL].iloc[0], 10.757)
        self.assertEqual(var_df[var_df[NT_VAR_COL] == VARIANTS[0][NT_VAR_COL]][
                             Candra.LIFD_SIGNIFICANCE_COL].iloc[0],  0.085632)
        self.assertEqual(var_df[var_df[NT_VAR_COL] == VARIANTS[0][NT_VAR_COL]][Candra.LIFD_CATEGORY_COL].iloc[0],
                         'Driver')
