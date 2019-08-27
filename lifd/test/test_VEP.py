
import unittest
import logging
import os
import pandas as pd
import shutil

from lifd.predictors.vep import Vep
from lifd.utils import NT_VAR_COL, FATHMM_KEY_COL

from lifd.test.utils_tests import VARIANTS

__author__ = 'Johannes REITER'
__date__ = 'Jan 19, 2019'

# get logger
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


class VepTest(unittest.TestCase):

    def setUp(self):
        # mirroring mini pandas dataframe with the essential information of variants
        self.var_df = pd.DataFrame(VARIANTS)
        self._delete_tmp_files()

    def tearDown(self):
        self._delete_tmp_files()

    def _delete_tmp_files(self):
        for sub in self.var_df.Subject.unique():
            input_fp = '{}{}'.format(sub, Vep.INPUT_SUFFIX)
            output_fp = '{}{}'.format(sub, Vep.OUTPUT_SUFFIX)

            if os.path.exists(input_fp) and os.path.isfile(input_fp):
                os.remove(input_fp)

            if os.path.exists(output_fp) and os.path.isfile(output_fp):
                os.remove(output_fp)

    def test_VEPAvailability(self):

        # is VEP installed?
        available = shutil.which(os.path.join(Vep.PATH, 'vep')) is not None
        self.assertTrue(available, 'VEP is not available')

    # @unittest.skip('Not yet implemented')
    def test_Vep(self):

        n_expected_cols = 18

        for sub in self.var_df.Subject.unique():
            input_fp = '{}{}'.format(sub, Vep.INPUT_SUFFIX)
            output_fp = '{}{}'.format(sub, Vep.OUTPUT_SUFFIX)
            # generate input file
            self.assertFalse(os.path.exists(input_fp), 'VEP input file with variants should not yet exist.')
            Vep.generate_input_file(
                input_fp, self.var_df[self.var_df.Subject == sub], var_key_col=NT_VAR_COL,
                chromosome_col='Chromosome', start_col='StartPosition',
                reference_col='ReferenceAllele', alternate_col='AlternateAllele')

            self.assertTrue(os.path.exists(input_fp), 'VEP input file with variants should exist.')

            # generate VEP output
            self.assertFalse(os.path.exists(output_fp), 'VEP output file should not yet exist.')
            Vep.run(input_fp, output_fp)
            self.assertTrue(os.path.exists(output_fp) and os.path.isfile(output_fp), 'VEP output file should exist.')

            if sub == 'Pat1':
                vep_df = Vep.read_results(output_fp)

                self.assertEqual(vep_df.shape, (1, n_expected_cols),
                                 'Dimensions of the dataframe checking: {}'.format(vep_df.shape))

                self.assertEqual(vep_df.loc[VARIANTS[0][NT_VAR_COL]][Vep.IMPACT_COL], 'MODERATE')
                self.assertEqual(vep_df.loc[VARIANTS[0][NT_VAR_COL]][FATHMM_KEY_COL], 'ENSP00000256078__G12V')

            elif sub == 'Pat2':
                vep_df = Vep.read_results(output_fp)

                self.assertEqual(vep_df.shape, (2, n_expected_cols),
                                 'Dimensions of the dataframe checking: {}'.format(vep_df.shape))

                self.assertEqual(vep_df.loc[VARIANTS[1][NT_VAR_COL]][Vep.IMPACT_COL], 'MODERATE')
                self.assertEqual(vep_df.loc[VARIANTS[1][NT_VAR_COL]][FATHMM_KEY_COL], 'ENSP00000341551__D351G')

    def test_Vep_get_results(self):

        # generate input file
        sub_name = 'Pat1'
        input_fp = '{}{}'.format(sub_name, Vep.INPUT_SUFFIX)
        output_fp = '{}{}'.format(sub_name, Vep.OUTPUT_SUFFIX)
        self.assertFalse(os.path.exists(input_fp), 'Cravat input file with variants should not yet exist.')
        self.assertFalse(os.path.exists(output_fp), 'Cravat output file should not yet exist.')

        n_rows, n_cols = self.var_df.shape

        # generate Vep output
        filter_cond = (self.var_df.Subject == sub_name)
        self.var_df = Vep.get_annotation(self.var_df, ds_name=sub_name, input_fp=input_fp, output_fp=output_fp,
                                         filter_condition=filter_cond)

        self.assertTrue(os.path.exists(input_fp), 'VEP input file with variants should exist.')
        self.assertTrue(os.path.exists(output_fp), 'VEP output file should exist.')
        # VEP is expected to add three new columns
        self.assertEqual((n_rows, n_cols + 3), self.var_df.shape,
                         'Dimensions of the dataframe checking: {}'.format(self.var_df.shape))

        self.assertEqual(self.var_df[self.var_df[NT_VAR_COL] == VARIANTS[0][NT_VAR_COL]][Vep.IMPACT_COL].iloc[0],
                         'MODERATE')
        self.assertEqual(self.var_df[self.var_df[NT_VAR_COL] == VARIANTS[0][NT_VAR_COL]][FATHMM_KEY_COL].iloc[0],
                         VARIANTS[0][FATHMM_KEY_COL])
