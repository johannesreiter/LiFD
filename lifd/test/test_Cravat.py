
import unittest
import os
import pandas as pd
import shutil

from lifd.predictors.cravat import Cravat
from lifd.utils import NT_VAR_COL
from lifd.settings import CHR_COL, POS_START_COL, REF_COL, ALT_COL

from lifd.test.utils_tests import CANCER_TYPE, VARIANTS, N_FUNC_VARS

__author__ = 'Johannes REITER'
__date__ = 'Jan 19, 2019'


INPUT_FP = '{}{}'.format(CANCER_TYPE, Cravat.INPUT_SUFFIX)
OUTPUT_DIR = './TMP'


class CravatTest(unittest.TestCase):

    def setUp(self):
        if os.path.exists(INPUT_FP) and os.path.isfile(INPUT_FP):
            os.remove(INPUT_FP)

        output_fp = os.path.join(OUTPUT_DIR, '{}{}'.format(CANCER_TYPE, Cravat.OUTPUT_SUFFIX))
        if os.path.exists(output_fp) and os.path.isfile(output_fp):
            os.remove(output_fp)

    def tearDown(self):
        if os.path.exists(INPUT_FP) and os.path.isfile(INPUT_FP):
            os.remove(INPUT_FP)
        assert not os.path.exists(INPUT_FP), 'Cravat input file has not been removed'

        output_fp = os.path.join(OUTPUT_DIR, '{}{}'.format(CANCER_TYPE, Cravat.OUTPUT_SUFFIX))
        if os.path.exists(output_fp) and os.path.isfile(output_fp):
            os.remove(output_fp)

        for f in os.listdir(OUTPUT_DIR):
            if f.startswith(f'{CANCER_TYPE}{Cravat.OUTPUT_SUFFIX}'):
                os.remove(os.path.join(OUTPUT_DIR, f))

    # @unittest.skip('Not yet implemented')
    def test_CravatAvailability(self):

        command = 'cravat'
        available = shutil.which(command) is not None
        self.assertTrue(available, 'cravat seems to be not installed at this location!')

    def test_Cravat(self):

        # mirroring mini pandas dataframe with the essential information of variants
        vars_df = pd.DataFrame(VARIANTS)

        # generate input file
        self.assertFalse(os.path.exists(INPUT_FP), 'Cravat input file with variants should not yet exist.')
        Cravat.generate_input_file(
            INPUT_FP, vars_df, chromosome_col=CHR_COL, position_col=POS_START_COL,
            reference_col=REF_COL, alternate_col=ALT_COL, subject_col='Subject')
        self.assertTrue(os.path.exists(INPUT_FP), 'Cravat input file with variants should exist.')

        # generate CRAVAT output
        output_fp = os.path.join(OUTPUT_DIR, '{}{}'.format(CANCER_TYPE, Cravat.OUTPUT_SUFFIX))
        self.assertFalse(os.path.exists(output_fp), 'Cravat output file should not yet exist.')
        Cravat.run(INPUT_FP, OUTPUT_DIR, CANCER_TYPE, reference_genome='hg19')
        self.assertTrue(os.path.exists(output_fp), 'Cravat output file should exist.')

        cravat_df = Cravat.read_results(output_fp, cancer_type=CANCER_TYPE)

        # Cravat produces results with 25 columns
        self.assertEqual((N_FUNC_VARS, 25), cravat_df.shape,
                         'Dimensions of the dataframe checking: {}'.format(cravat_df.shape))

        self.assertEqual(cravat_df.loc[VARIANTS[0][NT_VAR_COL]][Cravat.CP_SCORE_COL.replace('CT', CANCER_TYPE)],
                         0.936)
        self.assertEqual(cravat_df.loc[VARIANTS[0][NT_VAR_COL]][Cravat.CP_SCORE_CT_COL.replace('CT', CANCER_TYPE)],
                         0.913)
        self.assertEqual(cravat_df.loc[VARIANTS[0][NT_VAR_COL]][Cravat.CP_COL], 0)

    # @unittest.skip('Not yet implemented')
    def test_Cravat_get_results(self):

        # mirroring mini pandas dataframe with the essential information of variants
        var_df = pd.DataFrame(VARIANTS)
        n_cols = var_df.shape[1]

        # generate input file
        self.assertFalse(os.path.exists(INPUT_FP), 'Cravat input file with variants should not yet exist.')
        output_fp = os.path.join(OUTPUT_DIR, '{}{}'.format(CANCER_TYPE, Cravat.OUTPUT_SUFFIX))
        self.assertFalse(os.path.exists(output_fp), 'Cravat output file should not yet exist.')

        # generate CRAVAT output
        var_df = Cravat.get_annotation(var_df, CANCER_TYPE, input_fp=INPUT_FP, output_dir=OUTPUT_DIR)
        self.assertTrue(os.path.exists(INPUT_FP), 'Cravat input file with variants should exist.')
        self.assertTrue(os.path.exists(output_fp), 'Cravat output file should exist.')
        self.assertEqual(var_df.shape, (len(VARIANTS), n_cols + 4),
                         'Dimensions of the dataframe checking: {}'.format(var_df.shape))

        self.assertEqual(var_df[var_df[NT_VAR_COL] == VARIANTS[0][NT_VAR_COL]][
                             Cravat.CP_SCORE_COL].iloc[0], 0.936)
        self.assertEqual(var_df[var_df[NT_VAR_COL] == VARIANTS[0][NT_VAR_COL]][
                             Cravat.CP_SCORE_CT_COL].iloc[0],  0.913)
        self.assertEqual(var_df[var_df[NT_VAR_COL] == VARIANTS[0][NT_VAR_COL]][Cravat.CP_COL].iloc[0], 0)
