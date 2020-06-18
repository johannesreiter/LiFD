""" Helper functions to run Cravat/CHASM Carter et al, Cancer Research 2009 """
import logging
import os
import numpy as np
import pandas as pd
import csv
from subprocess import Popen, PIPE

from lifd.predictors.predictor import Predictor
from lifd.settings import HOME_DIR, PRD_DIR
from lifd.utils import add_column
from lifd.utils import FATHMM_KEY_COL, FUNC_COL

__author__ = 'Johannes Reiter'
__date__ = 'Jan 26, 2019'


# get logger
logger = logging.getLogger('lifd.{}'.format(__name__))


class FatHMM(Predictor):

    # FatHMM settings
    PATH = os.path.join(PRD_DIR, 'fathmm', 'cgi-bin')
    COMMAND = 'python2 fathmm.py -w Cancer {} {}' # 'python3 fathmm.py -w Cancer {} {}' if updating external fathmm to python 3
    INPUT_SUFFIX = '_fathmm_input.txt'
    OUTPUT_SUFFIX = '_fathmm_output.tsv'
    SCORE_COL = 'FATHMM_score'

    # default significance threshold (recommended -0.75)
    SIGNIFICANCE_TH = -0.75

    # which columns are required to run the predictor
    REQUIRED_COLS = [FATHMM_KEY_COL]

    # which columns are required to evaluate the prediction
    PREDICTION_COLS = [SCORE_COL]

    # reference genome not relevant for FatHMM
    SUPP_REF_GENOMES = None

    @staticmethod
    def predict_functionality(row):
        """
        Return tuple where first element indicates whether the tool predicted functionality
        and the second element indicates whether the tool produced a valid result
        :param row:
        :return: tuple with the No. tools predicting functionality and No. tools producing a result
        """

        if pd.isnull(row[FatHMM.SCORE_COL]):
            return 0.0, 0.0
        else:
            if row[FatHMM.SCORE_COL] <= FatHMM.SIGNIFICANCE_TH:
                return 1.0, 1.0
            else:
                return 0.0, 1.0

    @staticmethod
    def get_annotation(var_df, ds_name, input_fp=None, output_prefix='', output_dir=None, output_fp=None,
                       filter_condition=None, reference_genome=None):
        """
        Add or also run FatHMM for the variants in the given dataframe
        :param var_df: pandas dataframe with variants
        :param ds_name: name of dataset, possibly related to given filter;
                        used for default naming of input and output files
        :param input_fp: absolute path to input file
        :param output_prefix: should be None
        :param output_dir: should be None
        :param output_fp: absolute path to output file
        :param filter_condition: select only a subset of given dataframe, for example a cancer type or a subject
        :param reference_genome: no reference genome for FatHMM
        :return: extended pandas dataframe with FATHMM scores
        """
        if filter_condition is None:
            select_indices = np.where(var_df[FUNC_COL])[0]
        else:
            select_indices = np.where(filter_condition & var_df[FUNC_COL])[0]

        if len(select_indices) == 0:
            logger.warning('No variants selected for FATHMM analysis of case {}!'.format(ds_name))

        if output_dir is None:
            output_dir = os.path.join('.')
        if input_fp is None:
            input_fp = os.path.join(output_dir, '{}{}'.format(ds_name, FatHMM.INPUT_SUFFIX))
        if output_fp is None:
            output_fp = os.path.join(output_dir, '{}{}'.format(ds_name, FatHMM.OUTPUT_SUFFIX))

        if not os.path.isfile(input_fp):
            FatHMM.generate_input_file(input_fp, var_df.iloc[select_indices])

        if not os.path.isfile(output_fp):
            # run FatHMM
            FatHMM.run(os.path.abspath(input_fp), os.path.abspath(output_fp))

        fathmm_df = FatHMM.read_results(output_fp)

        # add FATHMM results to the dataframe
        add_column(var_df, select_indices, FATHMM_KEY_COL, FatHMM.SCORE_COL, fathmm_df,
                   as_value=True, sub_key=FatHMM.SCORE_COL, as_float=True)

        if (var_df.loc[select_indices][FatHMM.SCORE_COL].count()
                < 0.2 * var_df.loc[select_indices][FATHMM_KEY_COL].count()
                + min(3.0, var_df.loc[select_indices][FATHMM_KEY_COL].count() / 2.0)):
            logger.warning('Only {} FATHMM predictions found for case {} but {} IDs.'.format(
                    var_df.loc[select_indices][FatHMM.SCORE_COL].count(), ds_name,
                    var_df.loc[select_indices][FATHMM_KEY_COL].count()))

        return var_df

    @staticmethod
    def generate_input_file(fp, data_df, ensp_col='ENSP', protein_sub_col='AA_change'):
        """
        Generate input file for given nonsynonymous variants to run FATHMM based on pandas dataframe
        :param fp: path to file
        :param data_df: pandas dataframe with ENSP protein ID and amino acid change or FATHMM_ID
        :param ensp_col: ENSP protein identifier, e.g. ENSP00000351805
        :param protein_sub_col: amino acid change, e.g. R344Q
        """

        if ((ensp_col is None or ensp_col not in data_df.columns or protein_sub_col is None
             or protein_sub_col not in data_df.columns)
                and (FATHMM_KEY_COL is None or FATHMM_KEY_COL not in data_df.columns)):
            logger.error(
                'Either FATHMM ID column ({}) needs to be provided or protein identifier and amino acid change!')
            logger.error('Some of the give cols {}, {}, {} are not in given dataframe columns: {}'.format(
                ensp_col, protein_sub_col, FATHMM_KEY_COL, ', '.join(c for c in data_df.columns)))
            raise RuntimeError('FATHMM file generation needs protein and substitution identifiers columns!')

        if FATHMM_KEY_COL is not None and FATHMM_KEY_COL in data_df.columns:
            # extract information from identifier
            extract = True
            fathmm_id_col = data_df.columns.get_loc(FATHMM_KEY_COL)
            ensp_col = aa_change_col = None
        else:
            ensp_col = data_df.columns.get_loc(ensp_col)
            aa_change_col = data_df.columns.get_loc(protein_sub_col)
            fathmm_id_col = None
            extract = False

        with open(fp, 'w') as input_file:
            fathmm_writer = csv.writer(input_file, delimiter=' ')

            for index, row in data_df.iterrows():
                if extract:
                    if pd.isnull(row[fathmm_id_col]):
                        continue

                    fathmm_id = row[fathmm_id_col]
                    ensp, aa_sub = fathmm_id.split('__')

                else:
                    ensp = row[ensp_col]
                    aa_sub = row[aa_change_col]

                # write FATHMM input file
                fathmm_writer.writerow([ensp, aa_sub])

            logger.info('Generated FATHMM input file from dataframe with {} entries: {}'.format(
                data_df.shape[0], fp))

    @staticmethod
    def run(input_fp, output_fp):
        """
        Run tool FATHMM and start mysql server if not yet running
        :param input_fp:
        :param output_fp:
        :return:
        """
        cmd = FatHMM.COMMAND.format(input_fp, output_fp)
        logger.info('Call Fathmm with command {}'.format(cmd))
        #     exit_code = call(cmd, cwd=fathmm_path, shell=True)

        with Popen(cmd, stdout=PIPE, cwd=FatHMM.PATH, shell=True) as p:
            for line in p.stdout:
                logger.info(line)
            exit_code = p.poll()

        if exit_code == 0 or exit_code is None:
            logger.info('Successfully called FATHMM: {}'.format(output_fp))
        else:
            logger.info('Fathmm call was unsuccessful. Is Fathmm installed and running?')
            logger.info('Try starting mysql server with: mysql.server start or check with mysql.server status')

            # does not work with mariadb in cluster
            with Popen('mysql.server start', stdout=PIPE, shell=True) as p:
                # mysql.server restart/status/stop
                for line in p.stdout:
                    logger.info(line)
                exit_code = p.poll()

            if exit_code == 0:
                logger.info('Successfully started mysql server.')
                # try once again
                with Popen(cmd, stdout=PIPE, cwd=FatHMM.PATH, shell=True) as p:
                    for line in p.stdout:
                        logger.info(line)
                    exit_code = p.poll()

                if exit_code == 0:
                    logger.info('Successfully called Fathmm: {}'.format(output_fp))
                else:
                    logger.info('Fathmm call was unsuccessful. Is Fathmm installed?')
                    logger.info('Check status of mysql server with mysql.server status')
            else:
                logger.info('Could not start mysql server. Exit code: {}'.format(exit_code))
                raise RuntimeError('Could not run Fathmm or invoke mysql server! ')

            # ERROR! The server quit without updating PID file (/usr/local/var/mysql/<computer_name>.err
            # // Just empty macbook.local.err file then restart the server
            # rm /usr/local/var/mysql/<COMPUTER NAME>.local.err
            # mysql.server start

    @staticmethod
    def read_results(output_fp):
        """
        Read FATHMM results from given filepath
        :param output_fp:
        :return: dictionary with results; key=<ProteinID>__<Substitution>
                 values = ['Prediction', 'Score', 'Warning']
        """
        if os.path.isfile(output_fp):
            fathmm_df = pd.read_csv(output_fp, delimiter='\t', quoting=csv.QUOTE_NONE)
            fathmm_df[FATHMM_KEY_COL] = fathmm_df.apply(
                lambda row: '{}__{}'.format(row['Protein ID'], row.Substitution), axis=1)
            fathmm_df.drop_duplicates(subset=FATHMM_KEY_COL, keep='first', inplace=True)
            fathmm_df.set_index(FATHMM_KEY_COL, inplace=True)
            fathmm_df.rename(columns={'Score': FatHMM.SCORE_COL}, inplace=True)
            fathmm_df.dropna(subset=[FatHMM.SCORE_COL], inplace=True)
            logger.info('Loaded {} FATHMM annotated variants.'.format(len(fathmm_df)))
            return fathmm_df
        else:
            logger.error('Missing FATHMM output file {}'.format(output_fp))
