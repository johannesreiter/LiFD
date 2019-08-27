""" Helper functions to run Cravat/CHASM Carter et al, Cancer Research 2009; Tokheim & Karchin, 2018 """
import logging
import os
import numpy as np
import pandas as pd
import csv
from subprocess import Popen, PIPE

from lifd.predictors.predictor import Predictor
from lifd.utils import add_column, NAN, NT_VAR_COL, FUNC_COL, CT_COL

__author__ = 'Johannes Reiter'
__date__ = 'Jan 26, 2019'


# get logger
logger = logging.getLogger(__name__)


class Cravat(Predictor):

    INPUT_SUFFIX = '_cravat_input.tsv'
    OUTPUT_SUFFIX = '_cravat_output.xlsx'
    CRAVAT_COMMAND = 'cravat -n {} -l {} -t excel -a chasmplus chasmplus_{} -d {} {}'
    CP_COL = 'CHASMplus_Pvalue'
    CP_SCORE_COL = 'CHASMplus_Score'
    CP_SCORE_CT_COL = 'CHASMplus_CT_Score'
    CP_CT_COL = 'CHASMplus_CT_Pvalue'
    CP_CT_COL_COR = CP_CT_COL+'_corr'   # p-value corrected for multi-hypotheses testing

    # default q-value significance threshold
    CHASMPLUS_Q_TH = 0.1

    # which columns are required to run the predictor
    REQUIRED_COLS = [NT_VAR_COL, CT_COL]

    # which columns are required to evaluate the prediction
    PREDICTION_COLS = [CP_CT_COL_COR]

    # which reference genomes are supported
    SUPP_REF_GENOMES = ['hg19', 'hg20']

    @staticmethod
    def predict_functionality(row):
        """
        Return tuple where first element indicates whether the tool predicted functionality
        and the second element indicates whether the tool produced a valid result
        :param row: row representing a mutation from pandas dataframe
        :return: tuple with the No. tools predicting functionality and No. tools producing a result
        """

        if pd.isnull(row[Cravat.CP_CT_COL_COR]):
            return 0.0, 0.0
        else:
            if row[Cravat.CP_CT_COL_COR] <= Cravat.CHASMPLUS_Q_TH:
                return 1.0, 1.0
            else:
                return 0.0, 1.0

    @staticmethod
    def get_annotation(var_df, ds_name, input_fp=None, output_prefix='', output_dir=None, output_fp=None,
                       filter_condition=None, reference_genome='hg19'):
        """
        Add or also run CRAVAT and CHASMplus for the given dataframe
        :param var_df: pandas dataframe
        :param ds_name: name of dataset, possibly related to given filter;
                        used for default naming of input and output files
        :param input_fp: path to input file
        :param output_prefix: output naming prefix
        :param output_dir: path to output directory
        :param output_fp: if cravat_output_fp is given then CRAVAT will not be run and it is assumed to exist already
        :param filter_condition: select only a subset of given dataframe, for example a cancer type or a subject
        :param reference_genome: human reference genome, e.g. hg19 or hg20
        :return: extended pandas dataframe with CHASMPlus pancancer and cancer type specific scores and p-values
        """

        # Cravat can not convert genomic position for mitochondrial DNA
        if filter_condition is None:
            select_indices = np.where(var_df[FUNC_COL] & (var_df.Chromosome != 'MT'))[0]
        else:
            select_indices = np.where(filter_condition & var_df[FUNC_COL] & (var_df.Chromosome != 'MT'))[0]
            
        if len(select_indices) == 0:
            logger.warning('No variants selected for CRAVAT analysis of dataset {}!'.format(ds_name))

        if output_dir is None:
            output_dir = os.path.join('.')
        if input_fp is None:
            input_fp = os.path.join(output_dir, '{}{}'.format(ds_name, Cravat.INPUT_SUFFIX))

        if not os.path.isfile(os.path.abspath(input_fp)):
            Cravat.generate_input_file(
                input_fp, var_df.iloc[select_indices], chromosome_col='Chromosome',
                position_col='StartPosition', reference_col='ReferenceAllele',
                alternate_col='AlternateAllele', subject_col='Subject')
        else:
            logger.debug('CRAVAT input file for dataset {} already exists.'.format(ds_name))

        # get cancer type
        if len(var_df.iloc[select_indices][CT_COL].unique()) != 1:
            # var_df[var_df.Subject == sub_name][CT_COL].unique()
            raise RuntimeError('No cancer type could be inferred for dataset {}: {}'.format(
                ds_name, var_df.iloc[select_indices][CT_COL].unique()))

        cts = var_df.iloc[select_indices][CT_COL].unique()

        if len(cts) == 1:
            ct = cts[0]
        else:
            logger.error('No unique cancer type could be identified: f{cts}')
            ct = None

        if output_fp is None:
            output_fp = os.path.join(output_dir, f'{ds_name}{Cravat.OUTPUT_SUFFIX}')

            if not os.path.isfile(output_fp):
                # run CRAVAT
                Cravat.run(os.path.abspath(input_fp), output_dir, ct, prefix=output_prefix, sub_name=ds_name,
                           reference_genome=reference_genome)
            else:
                logger.debug('CRAVAT output file already exists: '.format(output_fp))

        if os.path.exists(output_fp) and os.path.isfile(output_fp):
            if output_fp.endswith('.xlsx'):

                cravat_df = Cravat.read_results(os.path.abspath(output_fp), cancer_type=ct)
                logger.debug('Read cravat results with {} entries: {}'.format(len(cravat_df), output_fp))

                add_column(var_df, select_indices, NT_VAR_COL, Cravat.CP_COL, cravat_df,
                           sub_key=Cravat.CP_COL, as_float=True)
                add_column(var_df, select_indices, NT_VAR_COL, Cravat.CP_SCORE_COL, cravat_df,
                           sub_key=Cravat.CP_SCORE_COL, as_float=True)

                if ct=='PANCAN':
                    add_column(var_df, select_indices, NT_VAR_COL, Cravat.CP_CT_COL, cravat_df,
                               sub_key=Cravat.CP_COL, as_float=True)
                    add_column(var_df, select_indices, NT_VAR_COL, Cravat.CP_SCORE_CT_COL, cravat_df,
                               sub_key=Cravat.CP_SCORE_COL, as_float=True)
                elif ct is not None and not (isinstance(ct, float) and np.isnan(ct)):
                    add_column(var_df, select_indices, NT_VAR_COL, Cravat.CP_CT_COL, cravat_df,
                               sub_key=Cravat.CP_CT_COL.replace('CT', ct), as_float=True)
                    add_column(var_df, select_indices, NT_VAR_COL, Cravat.CP_SCORE_CT_COL, cravat_df,
                               sub_key=Cravat.CP_SCORE_CT_COL.replace('CT', ct), as_float=True)
                else:
                    logger.warning('No cancer type given for dataset {}: {}'.format(ds_name, ct))

                if len(select_indices) < 10:
                    min_fraction = 0.0
                else:
                    min_fraction = 1.0

                if 'MutationEffect' in var_df.columns:
                    df = var_df.iloc[select_indices, :]
                    min_fraction *= 0.9 * len(df[df.MutationEffect == 'Substitution']) / len(select_indices)
                else:
                    min_fraction *= 0.3

                # check if at least a minimum number of predictions were found
                assert (var_df.iloc[select_indices, :][Cravat.CP_COL].count() >
                        min_fraction * len(select_indices)), \
                    'Only {} Cravat/CHASMplus predictions found for dataset {} but {} indices.'.format(
                        var_df.iloc[select_indices, :][Cravat.CP_COL].count(), ds_name, len(select_indices))

            elif output_fp.endswith('.tsv'):
                cravat_vars = Cravat.read_results(os.path.abspath(output_fp))

                # add CHASM results to the dataframe
                add_column(var_df, select_indices, NT_VAR_COL, 'Chasm', cravat_vars, sub_key=0, as_float=True)

                assert var_df[filter_condition]['Chasm'].count() > 0.3 * len(select_indices), \
                    'Only {} Cravat/Chasm predictions found for dataset {} but {} indices.'.format(
                        var_df[filter_condition]['Chasm'].count(), ds_name, len(select_indices))
        else:
            logger.error('Missing cravat/chasm output for dataset {}! No file {}'.format(ds_name, output_fp))
            raise RuntimeError(
                'Missing cravat/chasm output for dataset {}! No file {}'.format(ds_name, output_fp))

        return var_df

    @staticmethod
    def generate_input_file(input_fp, data_df, chromosome_col='Chromosome', position_col='StartPosition',
                            reference_col='ReferenceAllele', alternate_col='AlternateAllele', var_key_col=NT_VAR_COL,
                            subject_col=None):

        """
        Generate input file for given nonsynonymous variants to run Chasm based on pandas dataframe
        :param input_fp: path to input file
        :param data_df: pandas dataframe with Chromosome, Position, Ref, and Alt columns
        :param chromosome_col: output column name
        :param position_col: output start position column name
        :param alternate_col: column name where alternate allele is stored
        :param reference_col: column name where reference allele is stored
        :param var_key_col: name of column with variant identifier
        :param subject_col: optional name of column with subject names
        """

        if ((chromosome_col is None or position_col is None or reference_col or alternate_col is None)
                and var_key_col is not None):
            # extract information from variant key
            extract = True
            var_key_col = data_df.columns.get_loc(var_key_col)
            chr_col = pos_col = ref_col = alt_col = None
            logger.debug('Extract chromosome, position, and alternate allele information from variant key.')

        else:
            chr_col = data_df.columns.get_loc(chromosome_col)
            pos_col = data_df.columns.get_loc(position_col)
            ref_col = data_df.columns.get_loc(reference_col)
            alt_col = data_df.columns.get_loc(alternate_col)
            var_key_col = None
            extract = False

        if subject_col is not None:
            sub_col = data_df.columns.get_loc(subject_col)
        else:
            sub_col = None

        with open(input_fp, 'w') as file:
            logger.debug('Write substitution variants to CRAVAT/CHASM input file: {}'.format(input_fp))
            tsv_writer = csv.writer(file, delimiter='\t')
            for i, (index, row) in enumerate(data_df.iterrows()):

                if extract:
                    var_key = row[var_key_col]
                    chrom, pos, ref, alt = var_key.split('__')
                else:
                    chrom = row[chr_col]
                    pos = row[pos_col]
                    ref = row[ref_col]
                    alt = row[alt_col]

                tsv_writer.writerow([i, 'chr{}'.format(chrom), '{}'.format(pos), '+',
                                     '-' if ref == '' else ref,
                                     '-' if alt == '' else alt,
                                     '' if sub_col is None else row[sub_col]])

            logger.info('Generated CRAVAT input file from dataframe with {} entries: {}'.format(
                len(data_df), input_fp))

    @staticmethod
    def run(input_fp, output_dir, cancer_type, prefix='', sub_name=None, reference_genome='hg19'):
        """
        Run OpenCRAVAT and CHASMplus
        :param input_fp: path to input file
        :param output_dir: output directory
        :param sub_name: subject name
        :param cancer_type: cancer type
        :param prefix: output naming prefix
        :param reference_genome: human reference genome, e.g. hg19 or hg20
        """

        identifier = '{}{}{}'.format(prefix, cancer_type if sub_name is None else sub_name, Cravat.OUTPUT_SUFFIX)
        cmd = Cravat.CRAVAT_COMMAND.format(identifier, reference_genome, cancer_type, output_dir, input_fp)
        logger.info('Call CRAVAT with identifier {} with cancer type {} with command: {}'.format(
            identifier, cancer_type, cmd))

        with Popen(cmd, stdout=PIPE, shell=True) as p:
            for line in p.stdout:
                logger.debug(line)
            exit_code = p.poll()

        if exit_code is None or exit_code == 0:
            logger.info('Successfully called CRAVAT with output directory: {}'.format(output_dir))
        else:
            logger.error('CRAVAT call was unsuccessful. Is CRAVAT installed?')
            raise RuntimeError('CRAVAT quit with an error: {}'.format(exit_code))

    @staticmethod
    def read_results(cravat_output_fp, cancer_type=None):
        """
        Read Cravat/CHASM results from given filepath
        :param cravat_output_fp: path to cravat results file
        :param cancer_type: cancer type
        :return: pandas dataframe
        """
        if os.path.isfile(cravat_output_fp) and cravat_output_fp.endswith('.xlsx'):

            cravat_df = pd.read_excel(cravat_output_fp, sheet_name='Variant', header=1)
            # create var_key column and use it as index column
            # var key is <chromosome>__<start position>__<reference allele>__<alternate allele>
            cravat_df[NT_VAR_COL] = cravat_df.apply(
                lambda row: '{}__{}__{}__{}'.format(
                    row['Chrom.1'][3:], row['Position.1'], row['Ref Base'], row['Alt Base']), axis=1) # row['Hg19 Chrom'][3:], row['Hg19 Position'], row['Ref Base'], row['Alt Base']), axis=1)
            
            # # check for duplicates
            # if len(cravat_df[cravat_df.duplicated(NT_VAR_COL, keep=False)]) > 0:
            #     df = cravat_df[cravat_df.duplicated(NT_VAR_COL, keep=False)]
            #     # logger.error(df)
            #     logger.debug('Found {} duplicates in CRAVAT output according to the {}!'.format(
            #         len(df), NT_VAR_COL))

            cravat_df = cravat_df.set_index(NT_VAR_COL)

            cravat_df.rename(columns={'P-value': Cravat.CP_COL, 'Score': Cravat.CP_SCORE_COL}, inplace=True)
            cravat_df[Cravat.CP_COL].replace(to_replace='.', value=NAN, inplace=True)
            cravat_df[Cravat.CP_SCORE_COL].replace(to_replace='.', value=NAN, inplace=True)

            if cancer_type is not None and not (isinstance(cancer_type, float) and np.isnan(cancer_type)):
                cravat_df.rename(columns={'P-value.1': Cravat.CP_CT_COL.replace('CT', cancer_type),
                                          'Score.1': Cravat.CP_SCORE_CT_COL.replace('CT', cancer_type)}, inplace=True)

                if Cravat.CP_CT_COL.replace('CT', cancer_type) in cravat_df.columns:
                    cravat_df[Cravat.CP_CT_COL.replace('CT', cancer_type)].replace(to_replace='.', value=NAN,
                                                                                   inplace=True)
                else:
                    logger.error(f'Columns {Cravat.CP_CT_COL} was not found in given Cravat output file. '
                                 + 'Is the cancer type annotation correct?')

                if Cravat.CP_SCORE_CT_COL.replace('CT', cancer_type) in cravat_df.columns:
                    cravat_df[Cravat.CP_SCORE_CT_COL.replace('CT', cancer_type)].replace(to_replace='.', value=NAN,
                                                                                         inplace=True)
                else:
                    logger.error(
                        f'Columns {Cravat.CP_SCORE_CT_COL} was not found in given Cravat output file. '
                        + 'Is the cancer type annotation correct?')

            logger.info('Loaded {} CRAVAT annotated variants from file {}.'.format(len(cravat_df), cravat_output_fp))
            return cravat_df

        else:
            logger.error('Missing Cravat/Chasm output file {}'.format(cravat_output_fp))
            return None
