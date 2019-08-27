
import unittest
import logging
import os
import pandas as pd
import shutil

from lifd.predictors.fathmm import FatHMM
from lifd.utils import NT_VAR_COL, FATHMM_KEY_COL

from lifd.test.utils_tests import VARIANTS

__author__ = 'Johannes REITER'
__date__ = 'Jan 19, 2019'

# get logger
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


class FatHMMTest(unittest.TestCase):

    def setUp(self):
        # mirroring mini pandas dataframe with the essential information of variants
        self.var_df = pd.DataFrame(VARIANTS)
        self._delete_tmp_files()

    def tearDown(self):
        self._delete_tmp_files()

    def _delete_tmp_files(self):

        for sub in self.var_df.Subject.unique():
            input_fp = os.path.abspath('{}{}'.format(sub, FatHMM.INPUT_SUFFIX))
            output_fp = os.path.abspath('{}{}'.format(sub, FatHMM.OUTPUT_SUFFIX))

            if os.path.exists(input_fp) and os.path.isfile(input_fp):
                os.remove(input_fp)

            if os.path.exists(output_fp) and os.path.isfile(output_fp):
                os.remove(output_fp)

    def test_FatHMMAvailability(self):

        # is FatHMM installed?
        available = shutil.which('python2') is not None
        self.assertTrue(available, 'python2 is not available')

        # is candra script available?
        script_exists = os.path.exists(os.path.join(FatHMM.PATH, 'fathmm.py'))
        self.assertTrue(script_exists,
                        f'fathmm.py is not available at expected path {FatHMM.PATH}. '
                        + 'Change location in {}.'.format(os.path.join('src', 'lifd', 'predictors', 'fathmm.py')))

    # @unittest.skip('Not yet implemented')
    def test_FatHMM(self):

        for sub in self.var_df.Subject.unique():
            input_fp = os.path.abspath('{}{}'.format(sub, FatHMM.INPUT_SUFFIX))
            output_fp = os.path.abspath('{}{}'.format(sub, FatHMM.OUTPUT_SUFFIX))
            # generate input file
            self.assertFalse(os.path.exists(input_fp), 'FatHMM input file with variants should not yet exist.')
            FatHMM.generate_input_file(
                input_fp, self.var_df[self.var_df.Subject == sub])

            self.assertTrue(os.path.exists(input_fp), 'FatHMM input file with variants should exist.')

            # generate Candra output
            self.assertFalse(os.path.exists(output_fp), 'FATHMM output file should not yet exist.')
            FatHMM.run(input_fp, output_fp)
            self.assertTrue(os.path.exists(output_fp) and os.path.isfile(output_fp), 'FATHMM output file should exist.')

            if sub == 'Pat1':
                vep_df = FatHMM.read_results(output_fp)

                self.assertEqual(vep_df.shape, (1, 15),
                                 'Dimensions of the dataframe checking: {}'.format(vep_df.shape))

                self.assertEqual(vep_df.loc[VARIANTS[0][FATHMM_KEY_COL]][FatHMM.SCORE_COL], -2.32)

            elif sub == 'Pat2':
                vep_df = FatHMM.read_results(output_fp)

                self.assertEqual((2, 15), vep_df.shape,
                                 'Dimensions of the dataframe checking: {}'.format(vep_df.shape))

                self.assertEqual(vep_df.loc[VARIANTS[1][FATHMM_KEY_COL]][FatHMM.SCORE_COL], -5.74)

    # @unittest.skip('Not yet implemented')
    def test_FatHMM_get_results(self):

        # generate input file
        sub_name = 'Pat1'
        input_fp = '{}{}'.format(sub_name, FatHMM.INPUT_SUFFIX)
        output_fp = '{}{}'.format(sub_name, FatHMM.OUTPUT_SUFFIX)
        self.assertFalse(os.path.exists(input_fp), 'FatHMM input file with variants should not yet exist.')
        self.assertFalse(os.path.exists(output_fp), 'FatHMM output file should not yet exist.')
        n_rows, n_cols = self.var_df.shape

        # generate FatHMM output
        filter_cond = (self.var_df.Subject == sub_name)
        self.var_df = FatHMM.get_annotation(self.var_df, ds_name=sub_name, input_fp=input_fp, output_fp=output_fp,
                                            filter_condition=filter_cond)
        self.assertTrue(os.path.exists(input_fp), 'FATHMM input file with variants should exist.')
        self.assertTrue(os.path.exists(output_fp), 'FATHMM output file should exist.')
        # FATHMM adds one column
        self.assertEqual(self.var_df.shape, (n_rows, n_cols + 1),
                         'Dimensions of the dataframe checking: {}'.format(self.var_df.shape))

        self.assertEqual(
            self.var_df[self.var_df[NT_VAR_COL] == VARIANTS[0][NT_VAR_COL]][FatHMM.SCORE_COL].iloc[0], -2.32)
