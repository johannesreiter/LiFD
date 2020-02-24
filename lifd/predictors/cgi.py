""" Helper functions to run CGI (www.cancergenomeinterpreter.org; Tamborero et al, Genome Medicine 2018) """
import logging
import os
from collections import defaultdict
import numpy as np
import pandas as pd
import csv
import requests
import zipfile
import time

from lifd.predictors.predictor import Predictor
from lifd.utils import add_column, NT_VAR_COL, FUNC_COL, CT_COL

__author__ = 'Johannes Reiter'
__date__ = 'Jan 26, 2019'


# get logger
logger = logging.getLogger('lifd.{}'.format(__name__))


# NOTE: CGI cannot be run without a valid username and token!
class Cgi(Predictor):

    INPUT_SUFFIX = '_cgi_input.tsv'
    OUTPUT_SUFFIX = '_cgi_output.zip'
    # URL of CGI's REST API
    URL = 'https://www.cancergenomeinterpreter.org/api/v1'
    USER_ID = None
    TOKEN = None
    HEADERS = {'Authorization': '{} {}'.format(USER_ID, TOKEN)}

    DRIVER_COL = 'CGI_driver'
    DRIVER_GENE_COL = 'CGI_driver_gene'
    KNOWN_DRIVER_COL = 'CGI_known_driver'
    PREDICTED_DRIVER_COL = 'CGI_predicted_driver'
    DRIVER_MUT_PREDICTION_COL = 'CGI_driver_mut_prediction'
    GENE_ROLE_COL = 'CGI_gene_role'
    SOURCE_COL = 'CGI_driver_gene_source'
    CADD_PHRED_COL = 'cadd_phred'

    # which columns are required to run the predictor
    REQUIRED_COLS = [NT_VAR_COL]

    # which columns are required to evaluate the prediction
    PREDICTION_COLS = [DRIVER_COL]

    # which reference genomes are supported
    SUPP_REF_GENOMES = ['hg19']

    @staticmethod
    def set_login(user_id, token):
        """
        Set Login credentials
        :param user_id:
        :param token:
        """
        Cgi.USER_ID = user_id
        Cgi.TOKEN = token
        Cgi.HEADERS = {'Authorization': '{} {}'.format(Cgi.USER_ID, Cgi.TOKEN)}
        logger.info('Successfully set new username {} and token.'.format(Cgi.USER_ID))

    @staticmethod
    def predict_functionality(row):
        """
        Return tuple where first element indicates whether the tool predicted functionality
        and the second element indicates whether the tool produced a valid result
        :param row:
        :return: tuple with the No. tools predicting functionality and No. tools producing a result
        """

        if pd.isnull(row[Cgi.DRIVER_COL]):
            return 0.0, 0.0
        else:
            if row[Cgi.DRIVER_COL] == 'predicted' or row[Cgi.DRIVER_COL] == 'known':
                return 1.0, 1.0
            else:
                return 0.0, 1.0

    @staticmethod
    def get_annotation(var_df, ds_name, input_fp=None, output_prefix='', output_dir=None, output_fp=None,
                       filter_condition=None, reference_genome='hg19'):
        """
        Add CGI analysis results for whole dataset or a given subject's data
        :param var_df: pandas dataframe
        :param ds_name: name of dataset, possibly related to given filter;
                        used for default naming of input and output files
        :param input_fp: path to input file
        :param output_prefix: output naming prefix
        :param output_dir: path to output directory
        :param output_fp: path to expected output file
        :param filter_condition: select only a subset of given dataframe
                for example CGI could be run separately for different cancer types or subjects
        :param reference_genome: human reference genome, e.g. hg19
        :return: extended pandas dataframe with columns 'cadd_phred', 'CGI_driver', 'CGI_known_driver',
                 'CGI_predicted_driver', 'CGI_driver_gene', 'CGI_driver_gene_source', 'CGI_gene_role',
                 'CGI_driver_mut_prediction'
        """
        
        if filter_condition is None:
            select_indices = np.where(var_df[FUNC_COL])[0]
        else:
            select_indices = np.where(filter_condition & var_df[FUNC_COL])[0]

        if output_dir is None:
            output_dir = os.path.join('.')
        if input_fp is None:
            input_fp = os.path.join(output_dir, '{}{}'.format(ds_name, Cgi.INPUT_SUFFIX))
        if output_fp is None:
            output_fp = os.path.join(output_dir, '{}{}'.format(ds_name, Cgi.OUTPUT_SUFFIX))
            
        if not os.path.isfile(os.path.abspath(input_fp)):
            Cgi.generate_input_file(input_fp, var_df.iloc[select_indices])

        job_id = None
        if not os.path.isfile(output_fp) and os.path.isfile(input_fp):
            ct = var_df.iloc[select_indices][CT_COL].unique()[0]
            if ct == 'PANCAN':
                ct = 'CANCER'
            job_id = Cgi.run(os.path.abspath(input_fp), os.path.abspath(output_fp), ct, ds_name)

        if os.path.exists(output_fp) and zipfile.is_zipfile(output_fp):
            cgi_vars_df = Cgi.read_results(os.path.abspath(output_fp))
            # cgi_vars_df['driver']
            # add CGI results to the dataframe
            # interpreting CADD scores
            # >30 very high, >25 high, >20 medium, >10 low, <=10 very low
            add_column(var_df, select_indices, NT_VAR_COL, 'cadd_phred', cgi_vars_df,
                       sub_key='cadd_phred', as_float=True)
            # add three columns for whether variant is a known or predicted driver, or other (passenger)
            add_column(var_df, select_indices, NT_VAR_COL, Cgi.DRIVER_COL, cgi_vars_df,
                       sub_key='driver', as_float=False)
            add_column(var_df, select_indices, NT_VAR_COL, Cgi.KNOWN_DRIVER_COL,
                       cgi_vars_df[cgi_vars_df.driver == 'known'],
                       sub_key='driver', as_value=False)
            add_column(var_df, select_indices, NT_VAR_COL, Cgi.PREDICTED_DRIVER_COL,
                       cgi_vars_df[cgi_vars_df.driver != 'other'],
                       sub_key='driver', as_value=False)

            add_column(var_df, select_indices, NT_VAR_COL, Cgi.DRIVER_GENE_COL, cgi_vars_df,
                       sub_key='driver_gene', as_float=False)
            add_column(var_df, select_indices, NT_VAR_COL, Cgi.SOURCE_COL, cgi_vars_df,
                       sub_key='driver_gene_source', as_float=False)
            add_column(var_df, select_indices, NT_VAR_COL, Cgi.GENE_ROLE_COL, cgi_vars_df,
                       sub_key='gene_role', as_float=False)
            add_column(var_df, select_indices, NT_VAR_COL, Cgi.DRIVER_MUT_PREDICTION_COL, cgi_vars_df,
                       sub_key='driver_mut_prediction', as_float=False)

            # assert var_df['CGI_driver'].count() > 0.3 * len(select_indices), \
            #     'Only {} CGI predictions found for case {} but {} variants.'.format(
            #         var_df['CGI_driver'].count(), ds_name, len(select_indices)) # UNDO THIS COMENT
            # delete CGI request
            if job_id is not None:
                r = requests.delete('{}/{}'.format(Cgi.URL, job_id), headers=Cgi.HEADERS)
                r.json()
        else:
            logger.warning(f'Missing CGI output for case {ds_name} with input file: {input_fp}')
            existing_jobs = Cgi._get_existing_jobs()
            if ds_name not in existing_jobs.keys():
                most_recent_date = sorted(existing_jobs[ds_name].keys())[-1]
                job_id = existing_jobs[ds_name][most_recent_date]['id']
                job_info = Cgi.get_job_info(job_id)
                logger.warning(job_info)

        return var_df

    @staticmethod
    def generate_input_file(fp, data_df, chromosome_col='Chromosome', position_col='StartPosition',
                            reference_col='ReferenceAllele', alternate_col='AlternateAllele', var_key_col=NT_VAR_COL):
        """
        Generate input file for given nonsynonymous variants to run CGI (cancer genome interpreter) based
        on pandas dataframe
        :param fp: path to file
        :param data_df: pandas dataframe with Chromosome, Position, and Alt columns
        :param chromosome_col:
        :param position_col:
        :param alternate_col: column name where alternate allele is stored
        :param reference_col: column name where reference allele is stored
        :param var_key_col:
        """

        if ((chromosome_col is None or position_col is None or reference_col is None or alternate_col is None)
                and var_key_col is not None):
            # extract information from variant key
            extract = True
            var_key_col = data_df.columns.get_loc(var_key_col)
            chr_col = pos_col = ref_col = alt_col = None
        else:
            chr_col = data_df.columns.get_loc(chromosome_col)
            pos_col = data_df.columns.get_loc(position_col)
            alt_col = data_df.columns.get_loc(alternate_col)
            ref_col = data_df.columns.get_loc(reference_col)
            var_key_col = None
            extract = False

        with open(fp, 'w') as input_file:
            cgi_writer = csv.writer(input_file, delimiter='\t')
            cgi_writer.writerow(['chr', 'pos', 'ref', 'alt'])

            for index, row in data_df.iterrows():
                if extract:
                    var_key = row[var_key_col]
                    chrom, pos, ref, alt = var_key.split('__')
                else:
                    chrom = row[chr_col]
                    pos = row[pos_col]
                    ref = row[ref_col]
                    alt = row[alt_col]

                # write CGI input file
                cgi_writer.writerow(['chr' + str(chrom), pos, ref, alt])

                #         logger.info('Generated CGI input file from dataframe with {} entries: {}'.format(
                #             data_df.shape[0], fp))
            logger.info('Generated CGI input file from dataframe with {} entries: {}'.format(
                data_df.shape[0], fp))

    @staticmethod
    def get_job_info(job_id):
        """
        Load information about CGI analysis job with a given job ID through the CGI REST API
        :param job_id: job ID
        :return: dictionary with information
        """
        if Cgi.USER_ID is None or Cgi.TOKEN is None:
            raise RuntimeError('CGI requires a valid username {} and token.'.format(Cgi.USER_ID))

        r = requests.get('{}/{}'.format(Cgi.URL, job_id), headers=Cgi.HEADERS)
        #         {'metadata': {'user': 'reiter.j@gmail.com', 'title': 'Makohon-Pam01', 'id': '937c26fa23b0b1adacab',
        #                       'dataset': 'input.tsv', 'date': '2018-07-25 18:11:53', 'cancertype': 'PAAD'},
        #                       'status': 'Done'}
        job_info = r.json()
        return job_info
    
    @staticmethod
    def _get_existing_jobs():
        if Cgi.USER_ID is None or Cgi.TOKEN is None:
            raise RuntimeError('CGI requires a valid username {} and token.'.format(Cgi.USER_ID))
        # obtain information about jobs which were run already
        r = requests.get(Cgi.URL, headers=Cgi.HEADERS)
        job_ids = r.json()

        existing_jobs = defaultdict(dict)
        for job_id in job_ids:
            job_info = Cgi.get_job_info(job_id)
            if 'metadata' in job_info:
                existing_jobs[job_info['metadata']['title']][job_info['metadata']['date']] = job_info['metadata']
            else:
                logger.error('Key \"metadata\" not present in job_info!')
                logger.error(job_info)

        return existing_jobs

    @staticmethod
    def run(input_fp, output_fp, cancer_type, ds_id):
        """
        Run the tool CGI (cancer genome interpreter) if it has not been run yet and download results
        :param input_fp: path to input file
        :param output_fp: path to output file
        :param cancer_type: cancer type
        :param ds_id: dataset identifier (e.g., name of subject)
        """
        if Cgi.USER_ID is None or Cgi.TOKEN is None:
            raise RuntimeError('CGI requires a valid username {} and token.'.format(Cgi.USER_ID, Cgi.TOKEN))

        existing_jobs = Cgi._get_existing_jobs()
        
        # # subject has been previously analyzed with CGI
        # if ds_id in existing_jobs.keys():
        #     most_recent_date = sorted(existing_jobs[ds_id].keys())[-1]
        #     job_id = existing_jobs[ds_id][most_recent_date]['id']
        #     logger.info('CGI has been run previously for subject {} with ID {}. Use results from {}.'.format(
        #         ds_id, job_id, most_recent_date))

        # # no previous results found
        # else:
        
        # threading.Thread(
        #     target=requests.post, 
        #     args=(
        #         Cgi.URL,headers=Cgi.HEADERS, 
        #         files={'mutations': open(input_fp, 'rb')}, 
        #         data=Cgi.get_payload_run(cancer_type, ds_id)
        #     )
        # )
        
        r = requests.post(Cgi.URL,
                          headers=Cgi.HEADERS,
                          files={
                              'mutations': open(input_fp, 'rb'),
                          },
                          data=Cgi.get_payload_run(cancer_type, ds_id))
        job_id = r.json()
        if 'error_code' in job_id:
            logger.error('CGI job submission failed with error code {}: {}'.format(job_id['error_code'], job_id))
            raise RuntimeError('CGI job submission failed with response: {}'.format(job_id))
        else:
            logger.info('Successfully submitted job with ID {} to CGI for subject {}.'.format(job_id, ds_id))

        # download the results and write them to a zip file
        # check whether or not CGI completed the analysis
        job_info = Cgi.get_job_info(job_id)
        logger.debug(job_info)

        if 'error_code' in job_info:
            logger.warning('CGI reported an error code.')
            if 'message' in job_info:
                logger.error(job_info['message'])
            else:
                logger.error(job_info)

            logger.error('CGI analysis for case {} (id: {}) failed with error code {}'.format(
                ds_id, job_id, job_info['error_code']))
        
        else: 
            time_total = 0
            while Cgi.get_job_info(job_id)['status'] != 'Done':
                logger.info(
                    'CGI analysis for case {} (id: {}) is not yet completed. Time elapsed: {}'.format(
                    ds_id, job_id, time_total))
                time.sleep(20)
                time_total += 20
            
            logger.info(
                'CGI analysis for case {} (id: {}) is complete. Time elapsed: {}'.format(
                ds_id, job_id, time_total))
            payload = {'action': 'download'}
            r = requests.get('{}/{}'.format(Cgi.URL, job_id), headers=Cgi.HEADERS, params=payload)
            with open(output_fp, 'wb') as fd:
                fd.write(r._content)

        if os.path.isfile(output_fp):
            logger.info('Successfully downloaded results for subject {}.'.format(ds_id))
        else:
            raise RuntimeError(
                'CGI results were not successfully downloaded. Missing file {}.'.format(output_fp))
            
        # elif 'status' in job_info and job_info['status'] != 'Done':
        #     logger.info(
        #         'CGI analysis for case {} (id: {}) is not yet completed. Check again later! Status: {}'.format(
        #             ds_id, job_id, job_info['status']))

        # else:
        #     payload = {'action': 'download'}
        #     r = requests.get('{}/{}'.format(Cgi.URL, job_id), headers=Cgi.HEADERS, params=payload)
        #     with open(output_fp, 'wb') as fd:
        #         fd.write(r._content)

        #     if os.path.isfile(output_fp):
        #         logger.info('Successfully downloaded results for subject {}.'.format(ds_id))
        #     else:
        #         raise RuntimeError(
        #             'CGI results were not successfully downloaded. Missing file {}.'.format(output_fp))
        return job_id

    @staticmethod
    def read_results(output_fp):
        """
        Read results from the downloaded CGI (cancer genome interpreter) zip file into a pandas dataframe
        :param output_fp: path to a TSV (tab separated values) file generated by CGI
        :return: pandas dataframe
        """
        if os.path.isfile(output_fp):
            cgi_output_fn = 'mutation_analysis.tsv'

            archive = zipfile.ZipFile(output_fp, 'r')
            cgi_mut_file = archive.open(cgi_output_fn)

            cgi_df = pd.read_csv(cgi_mut_file, delimiter='\t')

            # create var_key column and use it as index column
            # var key is <chromosome>__<start position>__<reference allele>__<alternate allele>
            cgi_df[NT_VAR_COL] = cgi_df.apply(lambda row: '{}__{}__{}__{}'.format(
                row.chr, row.pos, row.ref, row.alt), axis=1)
            cgi_df.drop_duplicates(subset=NT_VAR_COL, keep='first', inplace=True)
            cgi_df.set_index(NT_VAR_COL, inplace=True)
            logger.info('Loaded {} CGI annotated variants from file {}.'.format(len(cgi_df), output_fp))
            return cgi_df
        else:
            logger.error('Missing CGI output file {}'.format(output_fp))
            raise RuntimeError('Missing CGI output file {}'.format(output_fp))

    @staticmethod
    def get_payload_run(cancer_type, title):
        """
        Format cancer type and data set identifier (e.g., subject name) for submission to CGI
        :param cancer_type: cancer type
        :param title: data set identifier (e.g., subject name)
        :return: dictionary formatted for CGI submission
        """
        # convert TCGA cancer type names to CGI names
        if cancer_type == 'COAD':
            ct = 'COC'
        elif cancer_type == 'BRCA':
            ct = 'BRCA'
        elif cancer_type == 'HGSC':  # High-grade serous ovarian cancer (HGSC)
            ct = 'OVSE'
        elif cancer_type == 'READ':
            ct = 'REC'
        elif cancer_type == 'CESC':  # Cervical Squamous Cell Carcinoma
            ct = 'CESC'
        elif cancer_type == 'KIRC':  # (RCCC) Renal clear cell
            ct = 'RCCC'
        elif cancer_type == 'KICH':  # Kidney Chromophobe
            ct = 'RCH'  # Renal chromophobe cell
        elif cancer_type == 'KIRP':  # Kidney renal papillary cell carcinoma
            ct = 'RPC'
        elif cancer_type == 'LIHC':  # Liver hepatocellular carcinoma
            ct = 'HC'  # Hepatic carcinoma
        elif cancer_type == 'LUAD':  # Lung Adenocarcinoma
            ct = 'LUAD'
        elif cancer_type == 'LUSC':  # Lung Squamous Cell Carcinoma
            ct = 'LUSC'
        elif cancer_type == 'UCEC':  # Uterine Corpus Endometrial Carcinoma
            ct = 'UCEC'
        elif cancer_type == 'SKCM':  # Skin Cutaneous Melanoma
            ct = 'CM'
        elif cancer_type == 'OV':  # Ovarian serous cystadenocarcinoma
            ct = 'OV'
        elif cancer_type == 'PAAD':  # pancreas adenocarcinoma
            ct = 'PAAD'
        elif cancer_type == 'PRAD':  # prostate adenocarcinoma
            ct = 'PRAD'
        elif cancer_type == 'STAD':  # stomach adenocarcinoma
            ct = 'STAD'
        else:
            logger.warning('Cancer type {} not yet implemented!'.format(cancer_type))
            ct = cancer_type

        return {'cancer_type': ct, 'title': title}
