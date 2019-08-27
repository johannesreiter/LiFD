
import unittest
import logging
import os
import shutil
import pandas as pd

from lifd.lifd import LiFD, LIFD_COL, LIFD_SUP_COL
from lifd.mutation_effects import annotate_effects
from lifd.utils import FUNC_COL, NT_VAR_COL, MUT_EFFECT_COL
from lifd.settings import HOTSPOTS_FP, COSMIC_VARS_FP, ONCOKB_ALLVARS_FP, ONCOGENIC_VARS_FP

from lifd.databases.hotspots_database import HotspotsDB
from lifd.databases.cosmic_db import CosmicDB
from lifd.databases.oncokb_database import OncoKBDB
from lifd.databases.cgi_database import CgiDB
from lifd.cgi_settings import CGI_USER_ID, CGI_TOKEN

from lifd.predictors.vep import Vep
from lifd.predictors.candra import Candra
from lifd.predictors.cravat import Cravat
from lifd.predictors.cgi import Cgi
from lifd.predictors.fathmm import FatHMM

from lifd.test.utils_tests import VARS_MIN_INPUT

# get logger
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

OUTPUT_DIR = os.path.join('..', 'output_integration_testing')


class IntegrationTest(unittest.TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        if os.path.exists(OUTPUT_DIR) and os.path.isdir(OUTPUT_DIR):
            shutil.rmtree(OUTPUT_DIR)

    def test_database_availability(self):

        # testing COSMIC integration
        cosmic_db = CosmicDB(COSMIC_VARS_FP, lazy_loading=False)
        n_cosmic_vars = 6210551
        # n_cosmic_samples = 28103
        self.assertGreaterEqual(n_cosmic_vars, sum(cosmic_db.db_df.Occurrences),
                                'Cosmic Genome Screen v89 should contain {} variants instead of {}.'.format(
                                n_cosmic_vars, sum(cosmic_db.db_df.Occurrences)))

        # testing cancer hotspots integration
        hs_db = HotspotsDB(HOTSPOTS_FP, lazy_loading=False)
        self.assertLessEqual(2547, hs_db.db_df.shape[0])
        self.assertEqual(41, hs_db.db_df.shape[1])

        # testing OncoKB integration
        oncokb_db = OncoKBDB(ONCOKB_ALLVARS_FP, lazy_loading=False)
        self.assertGreaterEqual(3787, oncokb_db.db_df.shape[0])
        self.assertEqual(10, oncokb_db.db_df.shape[1])

        # testing CGI integration
        cgi_db = CgiDB(ONCOGENIC_VARS_FP, lazy_loading=False)
        self.assertGreaterEqual(5601, len(cgi_db.db_df))
        self.assertEqual(9, cgi_db.db_df.shape[1])

        nt_var_key = '12__25398284__G__A'
        pt_var_key = 'KRAS__G12V'
        self.assertGreaterEqual(1411, cosmic_db.in_database(nt_var_key, None))
        self.assertTrue(hs_db.in_database(nt_var_key, pt_var_key))
        self.assertTrue(oncokb_db.in_database(nt_var_key, pt_var_key))
        self.assertTrue(cgi_db.in_database(nt_var_key, pt_var_key))

        nt_var_key = '12__25398284__G__@'
        pt_var_key = 'KRAS__G12@'
        self.assertFalse(cosmic_db.in_database(nt_var_key, pt_var_key))
        self.assertFalse(hs_db.in_database(nt_var_key, pt_var_key))
        self.assertFalse(oncokb_db.in_database(nt_var_key, pt_var_key))
        self.assertFalse(cgi_db.in_database(nt_var_key, pt_var_key))

        # 1111 cancer hotspots
        # dbs = [hs_db]
        # prds = []
        # l = LiFD(databases=dbs, predictors=prds)

    # @unittest.skip('Not yet implemented')
    def test_lifd_predictions(self):
        # set up databases for first phase of LiFD
        hs_db = HotspotsDB(HOTSPOTS_FP)
        cgi_db = CgiDB(ONCOGENIC_VARS_FP)
        Cgi.set_login(CGI_USER_ID, CGI_TOKEN)
        oncokb_db = OncoKBDB(ONCOKB_ALLVARS_FP)
        cosmic_db = CosmicDB(COSMIC_VARS_FP)
        dbs = [
            hs_db, cgi_db, oncokb_db, cosmic_db
               ]

        # predictors for second phase of LiFD
        prds = [Vep, Candra, Cravat, FatHMM, Cgi]

        lifd = LiFD(databases=dbs, predictors=prds)

        # load some test data
        # mirroring pandas dataframe with the essential information of variants
        var_df = pd.DataFrame(VARS_MIN_INPUT)

        var_df = annotate_effects(var_df)
        export_fn = 'integration_testing_results.xlsx'
        export_fp = os.path.join(OUTPUT_DIR, export_fn)
        # self.assertFalse(os.path.exists(export_fp) and os.path.isfile(export_fp))
        var_df = lifd.run_lifd(var_df, output_dir=OUTPUT_DIR, export_fn=export_fn)

        self.assertTrue(os.path.exists(export_fp) and os.path.isfile(export_fp))
        self.assertEqual(6, len(var_df[var_df[LIFD_COL]]))

        self.assertTrue(var_df[var_df['StartPosition'] == VARS_MIN_INPUT[0]['StartPosition']][LIFD_COL].iloc[0])
        self.assertEqual(4.4,
                         var_df[var_df['StartPosition'] == VARS_MIN_INPUT[0]['StartPosition']][LIFD_SUP_COL].iloc[0])
        self.assertEqual(1121, var_df[var_df['StartPosition'] == VARS_MIN_INPUT[0]['StartPosition']][
            cosmic_db.pred_col].iloc[0])

        nt_var_key = '19__48994757__-__G'
        self.assertEqual(var_df[var_df['StartPosition'] == VARS_MIN_INPUT[1]['StartPosition']][NT_VAR_COL].iloc[0],
                         nt_var_key)
        self.assertEqual(var_df[var_df['StartPosition'] == VARS_MIN_INPUT[1]['StartPosition']]['GeneSymbol'].iloc[0],
                         'LMTK3')
        # high functionality prediction but still no LIFD because mutation is not in a known driver gene
        self.assertEqual(1.5,
                         var_df[var_df['StartPosition'] == VARS_MIN_INPUT[1]['StartPosition']][LIFD_SUP_COL].iloc[0])
        self.assertFalse(var_df[var_df['StartPosition'] == VARS_MIN_INPUT[1]['StartPosition']][LIFD_COL].iloc[0])
        self.assertTrue(var_df[var_df['StartPosition'] == VARS_MIN_INPUT[1]['StartPosition']][FUNC_COL].iloc[0])

        nt_var_key = '17__7577568__C__T'
        self.assertEqual(var_df[var_df['StartPosition'] == VARS_MIN_INPUT[2]['StartPosition']][NT_VAR_COL].iloc[0],
                         nt_var_key)
        self.assertEqual(var_df[var_df['StartPosition'] == VARS_MIN_INPUT[2]['StartPosition']]['GeneSymbol'].iloc[0],
                         'TP53')
        self.assertTrue(var_df[var_df['StartPosition'] == VARS_MIN_INPUT[2]['StartPosition']][LIFD_COL].iloc[0])
        self.assertEqual(3.75,
                         var_df[var_df['StartPosition'] == VARS_MIN_INPUT[2]['StartPosition']][LIFD_SUP_COL].iloc[0])
        self.assertEqual(266, var_df[var_df['StartPosition'] == VARS_MIN_INPUT[2]['StartPosition']][
            cosmic_db.pred_col].iloc[0])

        nt_var_key = '18__48591889__A__G'
        self.assertEqual(var_df[var_df['StartPosition'] == VARS_MIN_INPUT[4]['StartPosition']][NT_VAR_COL].iloc[0],
                         nt_var_key)
        self.assertEqual(var_df[var_df['StartPosition'] == VARS_MIN_INPUT[4]['StartPosition']]['GeneSymbol'].iloc[0],
                         'SMAD4')
        self.assertEqual(var_df[var_df['StartPosition'] == VARS_MIN_INPUT[4]['StartPosition']]['Transcript_ID'].iloc[0],
                         'ENST00000342988')
        self.assertEqual(var_df[var_df['StartPosition'] == VARS_MIN_INPUT[4]['StartPosition']]['Protein_ID'].iloc[0],
                         'ENSP00000341551')
        self.assertEqual(var_df[var_df['StartPosition'] == VARS_MIN_INPUT[4]['StartPosition']][MUT_EFFECT_COL].iloc[0],
                         'Substitution')

        self.assertFalse(var_df[var_df['StartPosition'] == VARS_MIN_INPUT[7]['StartPosition']][FUNC_COL].iloc[0])
        self.assertFalse(var_df[var_df['StartPosition'] == VARS_MIN_INPUT[7]['StartPosition']][LIFD_COL].iloc[0])
        self.assertEqual(var_df[var_df['StartPosition'] == VARS_MIN_INPUT[7]['StartPosition']][LIFD_SUP_COL].iloc[0], 0)

    @unittest.skip('Not yet implemented')
    def test_error_handling(self):
        raise NotImplementedError


if __name__ == '__main__':
    unittest.main(warnings='ignore')
