
import unittest
import logging
import os
import pandas as pd
import requests

from lifd.predictors.cgi import Cgi
from lifd.utils import NT_VAR_COL
from lifd.cgi_settings import CGI_USER_ID, CGI_TOKEN

from lifd.test.utils_tests import VARIANTS, CANCER_TYPE, N_FUNC_VARS

__author__ = 'Johannes REITER'
__date__ = 'Jan 19, 2019'

# get logger
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


DS_ID_SUBMISSION = 'UNITTEST_SUBMISSION'   # dataset identifier
DS_ID_READING = 'UNITTEST_READING'   # dataset identifier


class CgiTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        Cgi.set_login(CGI_USER_ID, CGI_TOKEN)

    def setUp(self):
        # mirroring mini pandas dataframe with the essential information of variants
        self.var_df = pd.DataFrame(VARIANTS)
        # self._delete_tmp_files()

        self.n_cols = self.var_df.shape[1]

        # delete job with which submission is tested
        existing_jobs = Cgi._get_existing_jobs()
        if DS_ID_SUBMISSION in existing_jobs.keys():
            most_recent_date = sorted(existing_jobs[DS_ID_SUBMISSION].keys())[-1]
            job_id = existing_jobs[DS_ID_SUBMISSION][most_recent_date]['id']
            r = requests.delete('https://www.cancergenomeinterpreter.org/api/v1/{}'.format(job_id), headers=Cgi.HEADERS)
            print(r.json())

            existing_jobs = Cgi._get_existing_jobs()
            print('Following {} CGI jobs still exist: {}'.format(len(existing_jobs),
                                                                 ', '.join(job_id for job_id in existing_jobs)))

        # if DS_ID_READING not in existing_jobs.keys():
        #     logger.error('ERROR: Reading from CGI cannot be tested because adequate job does not exist.')

    def tearDown(self):
        files = ['{}{}'.format(DS_ID_READING, Cgi.INPUT_SUFFIX), '{}{}'.format(DS_ID_READING, Cgi.OUTPUT_SUFFIX),
                 '{}{}'.format(DS_ID_SUBMISSION, Cgi.INPUT_SUFFIX), '{}{}'.format(DS_ID_SUBMISSION, Cgi.OUTPUT_SUFFIX)]
        for file in files:
            if os.path.isfile(file):
                os.unlink(file)

    def test_CGI_connection(self):

        # connect to CGI
        # see https://www.cancergenomeinterpreter.org/rest_api for more details
        self.assertIsNotNone(Cgi.USER_ID, 'Valid username needed')
        self.assertIsNotNone(Cgi.TOKEN, 'Valid token needed')
        r = requests.get(Cgi.URL, headers=Cgi.HEADERS)
        job_ids = r.json()
        self.assertTrue(isinstance(job_ids, list))

    # @unittest.skip('Not yet implemented')
    def test_CGI_submission(self):
        """
        Tests generation of input TSV file and submission to CGI
        :return:
        """

        input_fp = '{}{}'.format(DS_ID_SUBMISSION, Cgi.INPUT_SUFFIX)
        output_fp = '{}{}'.format(DS_ID_SUBMISSION, Cgi.OUTPUT_SUFFIX)
        # generate input file
        self.assertFalse(os.path.exists(input_fp), 'CGI input file with variants should not yet exist.')
        Cgi.generate_input_file(input_fp, self.var_df)

        self.assertTrue(os.path.exists(input_fp), 'CGI input file with variants should exist.')

        # submit variants to CGI and obtain results
        n_existing_jobs = len(Cgi._get_existing_jobs().keys())
        self.assertFalse(os.path.exists(output_fp), 'CGI output file should not yet exist.')
        Cgi.run(input_fp, output_fp, CANCER_TYPE, DS_ID_SUBMISSION)
        self.assertEqual(len(Cgi._get_existing_jobs().keys()), n_existing_jobs+1, 'New job was created.')

        # downloading/reading results is separately tested

    # @unittest.skip('Not yet implemented')
    def test_CGI_reading(self):
        """
        Test downloading and reading results
        :return:
        """

        input_fp = '{}{}'.format(DS_ID_READING, Cgi.INPUT_SUFFIX)  # DS_ID_SUBMISSION
        output_fp = '{}{}'.format(DS_ID_READING, Cgi.OUTPUT_SUFFIX)  # DS_ID_SUBMISSION
        existing_jobs = Cgi._get_existing_jobs()

        self.assertFalse(os.path.exists(output_fp) and os.path.isfile(output_fp),
                         'CGI output file should not exist.')

        if DS_ID_READING not in existing_jobs.keys():  # DS_ID_READING
            logger.warning('Adequate job to test CGI reading is missing.')
            # generate input file
            Cgi.generate_input_file(input_fp, self.var_df)
            self.assertTrue(os.path.isfile(input_fp),
                            f'Input file has to exist to test CGI reading: {os.path.abspath(input_fp)}')

            # submit and download results file from CGI
            Cgi.run(input_fp, output_fp, CANCER_TYPE, DS_ID_READING)
            logger.info(f'Submitted missing job {DS_ID_READING} such that reading can be tested.')
        else:
            # download results file from CGI
            Cgi.run(input_fp, output_fp, CANCER_TYPE, DS_ID_READING)  # DS_ID_READING

        self.assertTrue(os.path.exists(output_fp) and os.path.isfile(output_fp), 'CGI output file should exist.')
        cgi_df = Cgi.read_results(output_fp)

        # test if duplicates got merged
        self.assertEqual(cgi_df.shape, (3, 50),
                         'Dimensions of the dataframe checking: {}'.format(cgi_df.shape))

        self.assertEqual(cgi_df.loc[VARIANTS[0][NT_VAR_COL]]['driver'], 'known')
        self.assertEqual(cgi_df.loc[VARIANTS[0][NT_VAR_COL]]['driver_gene'], 'tumor_driver')
        self.assertEqual(cgi_df.loc[VARIANTS[0][NT_VAR_COL]]['cadd_phred'], 29.9)
        self.assertEqual(cgi_df.loc[VARIANTS[0][NT_VAR_COL]]['cancer'], CANCER_TYPE)

        self.assertEqual(cgi_df.loc[VARIANTS[1][NT_VAR_COL]]['driver'], 'predicted')
        self.assertEqual(cgi_df.loc[VARIANTS[1][NT_VAR_COL]]['driver_gene'], 'tumor_driver')
        self.assertEqual(cgi_df.loc[VARIANTS[1][NT_VAR_COL]]['cadd_phred'], 29.4)

    # @unittest.skip('Not yet implemented')
    def test_Cgi_get_results(self):

        # generate input file
        input_fp = '{}{}'.format(DS_ID_READING, Cgi.INPUT_SUFFIX)    # DS_ID_READING
        output_fp = '{}{}'.format(DS_ID_READING, Cgi.OUTPUT_SUFFIX)  # DS_ID_READING
        # download results from CGI and obtain results
        self.assertFalse(os.path.exists(output_fp), 'CGI output file should not yet exist.')

        # generate CGI output
        sub_name = 'Pat1'
        filter_cond = (self.var_df.Subject == sub_name)
        self.var_df = Cgi.get_annotation(self.var_df, input_fp=input_fp, output_fp=output_fp,
                                         filter_condition=filter_cond, ds_name=sub_name)
        self.assertTrue(os.path.exists(output_fp) and os.path.isfile(output_fp), 'CGI output file should exist.')

        # CGI should add 8 columns to the dataframe
        self.assertEqual(self.var_df.shape, (len(VARIANTS), self.n_cols + 8),
                         'Dimensions of the dataframe checking: {}'.format(self.var_df.shape))

        self.assertEqual(self.var_df[Cgi.DRIVER_COL].count(), 1, f'Results were only added for subject {sub_name}')

        self.assertEqual(self.var_df[self.var_df[NT_VAR_COL] == VARIANTS[0][NT_VAR_COL]][Cgi.DRIVER_COL].iloc[0],
                         'known')
        self.assertEqual(self.var_df[self.var_df[NT_VAR_COL] == VARIANTS[0][NT_VAR_COL]][Cgi.DRIVER_GENE_COL].iloc[0],
                         'tumor_driver')
        self.assertEqual(self.var_df[self.var_df[NT_VAR_COL] == VARIANTS[0][NT_VAR_COL]][Cgi.CADD_PHRED_COL].iloc[0], 29.9)

        # add results for all subjects
        self.var_df = Cgi.get_annotation(self.var_df, CANCER_TYPE, input_fp=input_fp, output_fp=output_fp)
        self.assertEqual(self.var_df[Cgi.KNOWN_DRIVER_COL].count(), N_FUNC_VARS, 'Results for all subjects were added.')
        self.assertEqual(self.var_df[self.var_df[NT_VAR_COL] == VARIANTS[0][NT_VAR_COL]][Cgi.DRIVER_COL].iloc[0],
                         'known')
        self.assertEqual(self.var_df[self.var_df[NT_VAR_COL] == VARIANTS[0][NT_VAR_COL]][Cgi.DRIVER_COL].iloc[1],
                         'known')
