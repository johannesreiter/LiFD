""" Helper functions to run CanDrA (http://bioinformatics.mdanderson.org/main/CanDrA; Mao et al, PLoS ONE 2013) """
import logging
import os
import numpy as np
import pandas as pd
import csv
from subprocess import Popen, PIPE

from lifd.predictors.predictor import Predictor
from lifd.settings import HOME_DIR, PRD_DIR
from lifd.utils import add_column, NAN, NT_VAR_COL, FUNC_COL, CT_COL
from lifd.settings import CHR_COL, POS_START_COL, REF_COL, ALT_COL

__author__ = 'Johannes Reiter'
__date__ = 'Jan 26, 2019'


# get logger
logger = logging.getLogger('lifd.{}'.format(__name__))


class Candra(Predictor):

    # CanDrA path: only support hg19 (NCI Build 37)
    PATH = os.path.join(PRD_DIR, 'CanDrA.v+')
    DATABASE = os.path.join(PATH, 'database')
    INPUT_SUFFIX = '_candra_input.tsv'
    OUTPUT_SUFFIX = '_candra_output.tsv'
    CANDRA_COMMAND = 'perl open_candra.pl {} {} > {}'

    # columns in CanDrA output file
    SCORE_COL = 'CanDrA_Score'
    SIGNIFICANCE_COL = 'CanDrA_Significance'
    CATEGORY_COL = 'CanDrA_Category'
    # columns in LiFD output file
    LIFD_SCORE_COL = 'CanDrA_score'
    LIFD_SIGNIFICANCE_COL = 'CanDrA_sig'   # candra significance value
    LIFD_CATEGORY_COL = 'CanDrA_clf'

    # which columns are required to evaluate the predictor
    REQUIRED_COLS = [NT_VAR_COL, CT_COL]

    # which columns are required to evaluate the prediction
    PREDICTION_COLS = [LIFD_CATEGORY_COL, LIFD_SIGNIFICANCE_COL]

    # which reference genomes are supported
    SUPP_REF_GENOMES = ['hg19']

    # default significance threshold, if None, then only the category column is considered for the prediction
    SIGNIFICANCE_TH = 0.05

    @staticmethod
    def predict_functionality(row):
        """
        Return tuple where first element indicates whether the tool predicted functionality
        and the second element indicates whether the tool produced a valid result
        :param row:
        :return: tuple with the No. tools predicting functionality and No. tools producing a result
        """

        if pd.isnull(row[Candra.LIFD_CATEGORY_COL]) or pd.isnull(row[Candra.LIFD_SIGNIFICANCE_COL]):
            return 0.0, 0.0
        else:
            if row[Candra.LIFD_CATEGORY_COL] == 'Driver' \
                    and (Candra.SIGNIFICANCE_TH is None or row[Candra.LIFD_SIGNIFICANCE_COL] <= Candra.SIGNIFICANCE_TH):
                return 1.0, 1.0
            else:
                return 0.0, 1.0

    @staticmethod
    def get_annotation(var_df, ds_name, input_fp=None, output_prefix='', output_dir=None, output_fp=None,
                       filter_condition=None, reference_genome='hg19'):
        """
        Add or also run CanDrA for the given dataframe
        :param var_df: pandas dataframe
        :param ds_name: name of dataset, possibly related to given filter;
                        used for default naming of input and output files
        :param input_fp: path to input file
        :param output_prefix: output naming prefix
        :param output_dir: path to output directory
        :param output_fp: if output_fp is given then CanDrA will not be run and it is assumed to exist already
        :param filter_condition: select only a subset of given dataframe, for example a cancer type or a subject
        :param reference_genome: human reference genome, e.g. hg19 or hg20
        :return: extended pandas dataframe with CanDrA scores
        """

        if filter_condition is None:
            select_indices = np.where(var_df[FUNC_COL])[0]
        else:
            select_indices = np.where(filter_condition & var_df[FUNC_COL])[0]

        if len(select_indices) == 0:
            logger.warning('No variants selected for CanDrA analysis of case {}!'.format(ds_name))

        if output_dir is None:
            output_dir = os.path.join('.')
        if input_fp is None:
            input_fp = os.path.join(output_dir, '{}{}'.format(ds_name, Candra.INPUT_SUFFIX))
        if output_fp is None:
            output_fp = os.path.join(output_dir, '{}{}'.format(ds_name, Candra.OUTPUT_SUFFIX))

        if not os.path.isfile(input_fp):
            Candra.generate_input_file(
                input_fp, var_df.iloc[select_indices], chromosome_col='Chromosome',
                position_col='StartPosition', reference_col='ReferenceAllele',
                alternate_col='AlternateAllele')

        if not os.path.isfile(output_fp):
            # run CanDrA
            if len(var_df.iloc[select_indices][CT_COL].unique()) != 1:
                # var_df[var_df.Subject == sub_name][CT_COL].unique()
                raise RuntimeError('No cancer type could be inferred for subject {}: {} {}'.format(
                    ds_name, var_df.iloc[select_indices][CT_COL].unique(),
                    len(var_df.iloc[select_indices][CT_COL].unique())))

            ct = var_df.iloc[select_indices][CT_COL].unique()[0]
            Candra.run(os.path.abspath(input_fp), os.path.abspath(output_fp), ct)

        candra_df = Candra.read_results(output_fp)

        # add CanDrA results to the dataframe
        Candra.add_amino_acid_change_column(var_df, candra_df)
        add_column(var_df, select_indices, NT_VAR_COL, Candra.LIFD_SCORE_COL, candra_df,
                   as_value=True, sub_key=Candra.SCORE_COL, as_float=True)
        add_column(var_df, select_indices, NT_VAR_COL, Candra.LIFD_CATEGORY_COL, candra_df, as_value=True,
                   sub_key=Candra.CATEGORY_COL)
        var_df['CanDrA_clf'].replace('nan', NAN, inplace=True)
        add_column(var_df, select_indices, NT_VAR_COL, Candra.LIFD_SIGNIFICANCE_COL, candra_df,
                   as_value=True, sub_key=Candra.SIGNIFICANCE_COL, as_float=True)

        return var_df

    @staticmethod
    def generate_input_file(input_fp, data_df, chromosome_col=CHR_COL, position_col=POS_START_COL,
                            reference_col=REF_COL, alternate_col=ALT_COL, var_key_col=NT_VAR_COL):
        """
        Generate input file for given nonsynonymous variants to run CanDrA based on pandas dataframe
        :param input_fp: path to file
        :param data_df: pandas dataframe with Chromosome, Position, and Alt columns
        :param chromosome_col:
        :param position_col:
        :param reference_col: column name where reference allele is stored
        :param alternate_col: column name where alternate allele is stored
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
            ref_col = data_df.columns.get_loc(reference_col)
            alt_col = data_df.columns.get_loc(alternate_col)
            var_key_col = None
            extract = False

        with open(input_fp, 'w') as input_file:
            candra_writer = csv.writer(input_file, delimiter='\t')

            for index, row in data_df.iterrows():
                if extract:
                    var_key = row[var_key_col]
                    chrom, pos, ref, alt = var_key.split('__')
                else:
                    chrom = row[chr_col]
                    pos = row[pos_col]
                    ref = row[ref_col]
                    alt = row[alt_col]

                # write CanDrA input file
                # is there strand information available?
                candra_writer.writerow([chrom, pos, ref, alt, '+'])

            logger.info('Generated CanDrA input file from dataframe with {} entries: {}'.format(
                data_df.shape[0], input_fp))

    @staticmethod
    def run(input_fp, output_fp, cancer_type):
        """
        Run the tool CanDrA
        :param input_fp: path to input file
        :param output_fp: path to output file
        :param cancer_type: Supported cancer types: BLCA, BRCA, CRC, CSCC, EC, GBM, GENERAL, KIRC,
                            LUAD, LUSMACC, LUSQUCC, MDB, MEL, OVC, PRCA, SCSC
        """

        downloaded_candra_cancer_types = [d for d in os.listdir(Candra.DATABASE) if
                                          os.path.isdir(os.path.join(Candra.DATABASE, d))]
        logger.debug('CanDrA database for the following cancer types found: {}'.format(
            ', '.join(ct for ct in downloaded_candra_cancer_types)))

        # convert TCGA cancer type names to CanDrA names
        if cancer_type == 'COAD' or cancer_type == 'READ':
            ct = 'CRC'
        elif cancer_type == 'CESC':  # Cervical Squamous Cell Carcinoma
            ct = 'CSCC'
        elif cancer_type == 'LUSC':
            ct = 'LUSQUCC'  # Lung Squamous Cell Carcinoma
        elif cancer_type == 'UCEC':  # Uterine Corpus Endometrial Carcinoma
            ct = 'EC'
        elif cancer_type == 'SKCM':  # Skin Cutaneous Melanoma
            ct = 'MEL'
        elif cancer_type == 'OV':  # Ovarian serous cystadenocarcinoma
            ct = 'OVC'
        elif cancer_type == 'PRAD':  # prostate
            ct = 'PRCA'
        elif ((isinstance(cancer_type, float) and np.isnan(cancer_type))
              or cancer_type == 'ACC'   # Adrenocortical carcinoma
              or cancer_type == 'THCA'  # thyroid carcinoma
              or cancer_type == 'PAAD' or cancer_type == 'STAD'):  # pancreatic or gastric/stomach
            logger.warning('Cancer type {} is not supported by CanDrA. Using GENERAL instead.'.format(cancer_type))
            ct = 'GENERAL'
        else:
            # BLCA - bladder urothelial carcinoma
            # GBM - Glioblastoma multiforme
            ct = cancer_type

        if ct not in downloaded_candra_cancer_types:
            logger.warning('Given cancer type {} is not found in CanDrA database: {}'.format(
                ct, ', '.join(ct for ct in downloaded_candra_cancer_types)))
            ct = 'GENERAL'
            # raise RuntimeError('Given cancer type {} is not (yet) available in CanDrA. Check '.format(ct)
            #                    + 'http://bioinformatics.mdanderson.org/main/CanDrA')

        cmd = Candra.CANDRA_COMMAND.format(ct, os.path.abspath(input_fp), os.path.abspath(output_fp))
        logger.info('Call CanDrA for cancer type {} with command {}'.format(ct, cmd))

        with Popen(cmd, stdout=PIPE, cwd=Candra.PATH, shell=True) as p:
            for line in p.stdout:
                logger.info(line)
            exit_code = p.poll()

        if exit_code == 0:
            logger.info('Successfully called CanDrA: {}'.format(output_fp))
        else:
            logger.error('CanDrA call was unsuccessful. Is CanDrA installed?')
            raise RuntimeError('CanDrA quit with an error: {}'.format(exit_code))

    @staticmethod
    def read_results(output_fp):
        """
        Read CanDrA results from given filepath
        :param output_fp:
        :return: dictionary with results; key=<CHR>__<POSITION>__<ALT_ALLELE>
                 values = ['CONSEQUENCE', 'AAS', 'AAS_Location', 'CanDrA_Score',
                 'CanDrA_Category', 'CanDrA_Significance']
        """
        if os.path.isfile(output_fp):
            candra_df = pd.read_csv(output_fp, delimiter='\t')
            logger.info('Read CanDrA output file with {} lines: {}'.format(len(candra_df), output_fp))
            # create var_key column and use it as index column
            # var key is <chromosome>__<start position>__<reference allele>__<alternate allele>
            candra_df[NT_VAR_COL] = candra_df.apply(
                lambda row: '{}__{}__{}__{}'.format(
                    row['Chrom'], row['Coordinate'], row['Ref_Allele'], row['Mut_Allele']), axis=1)

            # # check for duplicates
            # if len(candra_df[candra_df.duplicated(NT_VAR_COL, keep=False)]) > 0:
            #     df = candra_df[candra_df.duplicated(NT_VAR_COL, keep=False)]
            #     # logger.error(df)
            #     logger.debug('Found {} duplicates in CRAVAT output according to the {}!'.format(
            #         len(df), NT_VAR_COL))

            candra_df.drop_duplicates(subset=NT_VAR_COL, keep='first', inplace=True)
            candra_df.set_index(NT_VAR_COL, inplace=True)

            logger.info('Loaded {} CanDrA annotated variants from file {}.'.format(len(candra_df), output_fp))
            return candra_df
        else:
            logger.error('Missing CanDrA output file {}'.format(output_fp))
            raise RuntimeError('Missing CanDrA output file {}'.format(output_fp))

    @staticmethod
    def add_amino_acid_change_column(var_df, candra_vars):
        """
        Add amino acid change column ('AA_change') to a dataframe from a given CanDrA output dictionary
        :param var_df:
        :param candra_vars:
        """
        # export mutation in protein formatting first
        # http://www.hgmd.cf.ac.uk/docs/mut_nom.html
        aa_changes = list()

        # if key is not in given dictionary with data, use preexisting value
        col_name = 'AA_change'
        if col_name in var_df.columns:
            preexisting_value = True
        else:
            preexisting_value = False

        for index, row in var_df.iterrows():
            var_key = row[NT_VAR_COL]
            if var_key in candra_vars and len(candra_vars[var_key][1]) == 3 and candra_vars[var_key][2] != '':
                aas = candra_vars[var_key][1].split('/')
                aa_change = aas[0] + candra_vars[var_key][2] + aas[1]
                aa_changes.append(aa_change)

            # do not replace preexisting value if key was not provided in the dictionary candra_vars
            elif preexisting_value:
                aa_changes.append(row[col_name])

            else:
                aa_changes.append(NAN)

        var_df[col_name] = pd.Series(aa_changes, index=var_df.index)
