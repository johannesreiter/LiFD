
import unittest
import pandas as pd

from lifd.predictors.cravat import Cravat
from lifd.utils import NT_VAR_COL, PT_VAR_COL, FATHMM_KEY_COL, MUT_EFFECT_COL, FUNC_COL
from lifd.mutation_effects import annotate_effects

from lifd.test.utils_tests import CANCER_TYPE, VARS_MIN_INPUT

__author__ = 'Johannes REITER'
__date__ = 'Jan 19, 2019'


INPUT_FP = '{}{}'.format(CANCER_TYPE, Cravat.INPUT_SUFFIX)
OUTPUT_DIR = '.'


class MutationEffectTest(unittest.TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    # @unittest.skip('Not yet implemented')
    def test_annotate_effects(self):

        # mirroring pandas dataframe with the essential information of variants
        var_df = pd.DataFrame(VARS_MIN_INPUT)

        n_rows, n_cols = var_df.shape

        new_cols = [NT_VAR_COL, PT_VAR_COL, 'GeneSymbol', 'Transcript_ID', 'Protein_ID', FATHMM_KEY_COL,
                    MUT_EFFECT_COL, FUNC_COL]
        for c in new_cols:
            self.assertNotIn(c, var_df.columns)

        var_df = annotate_effects(var_df)

        self.assertEqual((n_rows, n_cols + len(new_cols)), var_df.shape)

        for c in new_cols:
            self.assertIn(c, var_df.columns)

        nt_var_key = '12__25398284__C__A'
        self.assertEqual(var_df[var_df['StartPosition'] == VARS_MIN_INPUT[0]['StartPosition']][NT_VAR_COL].iloc[0],
                         nt_var_key)
        self.assertEqual(var_df[var_df['StartPosition'] == VARS_MIN_INPUT[0]['StartPosition']]['GeneSymbol'].iloc[0],
                         'KRAS')

        nt_var_key = '19__48994757__-__G'
        self.assertEqual(var_df[var_df['StartPosition'] == VARS_MIN_INPUT[1]['StartPosition']][NT_VAR_COL].iloc[0],
                         nt_var_key)
        self.assertEqual(var_df[var_df['StartPosition'] == VARS_MIN_INPUT[1]['StartPosition']]['GeneSymbol'].iloc[0],
                         'LMTK3')

        nt_var_key = '17__7577568__C__T'
        self.assertEqual(var_df[var_df['StartPosition'] == VARS_MIN_INPUT[2]['StartPosition']][NT_VAR_COL].iloc[0],
                         nt_var_key)
        self.assertEqual(var_df[var_df['StartPosition'] == VARS_MIN_INPUT[2]['StartPosition']]['GeneSymbol'].iloc[0],
                         'TP53')

        nt_var_key = '18__48591889__A__G'
        self.assertEqual(var_df[var_df['StartPosition'] == VARS_MIN_INPUT[4]['StartPosition']][NT_VAR_COL].iloc[0],
                         nt_var_key)
        self.assertEqual(var_df[var_df['StartPosition'] == VARS_MIN_INPUT[4]['StartPosition']]['GeneSymbol'].iloc[0],
                         'SMAD4')
        self.assertEqual(var_df[var_df['StartPosition'] == VARS_MIN_INPUT[4]['StartPosition']]['Transcript_ID'].iloc[0],
                         'ENST00000342988')
        self.assertEqual(var_df[var_df['StartPosition'] == VARS_MIN_INPUT[4]['StartPosition']]['Protein_ID'].iloc[0],
                         'ENSP00000341551')
        self.assertEqual(var_df[var_df['StartPosition'] == VARS_MIN_INPUT[4]['StartPosition']][FATHMM_KEY_COL].iloc[0],
                         'ENSP00000341551__D351G')
        self.assertEqual(var_df[var_df['StartPosition'] == VARS_MIN_INPUT[4]['StartPosition']][MUT_EFFECT_COL].iloc[0],
                         'Substitution')
        self.assertEqual(var_df[var_df['StartPosition'] == VARS_MIN_INPUT[4]['StartPosition']][FUNC_COL].iloc[0], True)
