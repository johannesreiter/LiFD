""" Helper functions to run VEP (McLaren et al, Genome Biology 2016) """
import logging
import os
import numpy as np
import pandas as pd
import csv
from subprocess import Popen, PIPE

from lifd.predictors.predictor import Predictor
from lifd.settings import HOME_DIR, VEP_CACHE, PRD_DIR, VEP_DIR
from lifd.utils import NT_VAR_COL, FUNC_COL, FATHMM_KEY_COL, PP_SCORE_COL, SIFT_SCORE_COL, add_column, NAN

__author__ = 'Johannes Reiter'
__date__ = 'Jan 26, 2019'


# get logger
logger = logging.getLogger('lifd.{}'.format(__name__))


class Vep(Predictor):

    # VEP configuration
    PATH = os.path.join(VEP_DIR)
    INPUT_SUFFIX = '_vep_input.tsv'
    OUTPUT_SUFFIX = '_vep_output.tsv'
    COMMAND = ('./vep --species homo_sapiens --dir {} --assembly {} --cache_version 99 --no_progress '
               + '--offline --no_stats --sift b --symbol --numbers --domains --gene_phenotype --pubmed --protein '
               + '--variant_class --shift_hgvs 1 --check_existing --no_escape --failed 1 --minimal --allele_number '
               + '--pick_order canonical,tsl,biotype,rank,ccds,length --input_file {} --output_file {} '
               + '--polyphen b --buffer_size 50 --force_overwrite --pick')
    IMPACT_COL = 'VEP_impact'

    # default significance threshold, if None, then only the category column is considered for the prediction
    SIGNIFICANCE_TH = 0.05

    # which columns are required to run the predictor
    REQUIRED_COLS = [NT_VAR_COL]

    # which columns are required to evaluate the prediction
    PREDICTION_COLS = [IMPACT_COL]

    # which reference genomes are supported
    SUPP_REF_GENOMES = ['hg19', 'hg20']

    @staticmethod
    def predict_functionality(row):
        """
        Return tuple where first element indicates whether the tool predicted functionality
        and the second element indicates whether the tool produced a valid result
        :param row:
        :return: tuple with the No. tools predicting functionality and No. tools producing a result
        """

        if pd.isnull(row[Vep.IMPACT_COL]):
            return 0.0, 0.0
        else:
            if row[Vep.IMPACT_COL] == 'HIGH':
                return 1.0, 1
            # elif row[Vep.IMPACT_COL] == 'MODERATE':
            #     return 0.5, 1
            else:
                return 0.0, 1

    @staticmethod
    def get_annotation(var_df, ds_name, input_fp=None, output_prefix='', output_dir=None, output_fp=None,
                       filter_condition=None, reference_genome='hg19'):
        """
        Run VEP and FATHMM for given subject's data and add the results to the given dataframe
        :param var_df: pandas dataframe with input data
        :param ds_name: name of dataset, possibly related to given filter;
                        used for default naming of input and output files
        :param input_fp: path to VEP input file
        :param output_prefix: output naming prefix
        :param output_dir: path to output directory
        :param output_fp: path to VEP output file
        :param filter_condition: select only a subset of given dataframe
        :param reference_genome: human reference genome, e.g. hg19 or hg20
        :return: extended pandas dataframe with columns 'PolyPhen', 'Sift', 'impact', 'consequence', 'FATHMM_ID',
        """
        if filter_condition is None:
            select_indices = np.where(var_df[FUNC_COL])[0]
        else:
            select_indices = np.where(filter_condition & (var_df[FUNC_COL]))[0]

        if output_dir is None:
            output_dir = os.path.join('.')
        if input_fp is None:
            input_fp = os.path.join(output_dir, '{}{}'.format(ds_name, Vep.INPUT_SUFFIX))
        if output_fp is None:
            output_fp = os.path.join(output_dir, '{}{}'.format(ds_name, Vep.OUTPUT_SUFFIX))

        if not os.path.isfile(input_fp):
            Vep.generate_input_file(input_fp, var_df.iloc[select_indices],
                                    chromosome_col='Chromosome', gene_col='GeneSymbol',
                                    start_col='StartPosition', end_col='EndPosition',
                                    reference_col='ReferenceAllele', alternate_col='AlternateAllele')

        if not os.path.isfile(output_fp) and os.path.isfile(input_fp):
            Vep.run(os.path.abspath(input_fp), os.path.abspath(output_fp), reference_genome=reference_genome)

        # Extract relevant information from VEP files and generate FATHMM input files if necessary
        if os.path.isfile(output_fp):
            vep_df = Vep.read_results(output_fp)

            if vep_df is not None:
                # polyphen
                add_column(var_df, select_indices, NT_VAR_COL, PP_SCORE_COL, vep_df, sub_key=PP_SCORE_COL,
                           as_value=True, as_float=True)
                add_column(var_df, select_indices, NT_VAR_COL, SIFT_SCORE_COL, vep_df, sub_key=SIFT_SCORE_COL,
                           as_value=True, as_float=True)
                add_column(var_df, select_indices, NT_VAR_COL, Vep.IMPACT_COL, vep_df, sub_key=Vep.IMPACT_COL,
                           as_value=True)
                if FATHMM_KEY_COL not in var_df.columns:
                    add_column(var_df, select_indices, NT_VAR_COL, FATHMM_KEY_COL, vep_df, sub_key=FATHMM_KEY_COL,
                               as_value=True)
                    var_df.replace({FATHMM_KEY_COL: 'nan'}, np.nan, inplace=True)

                # assert var_df[filter_condition][Vep.IMPACT_COL].count() > 0.3 * len(select_indices) - 2, \
                #     'Only {} VEP impact predictions found for case {} but {} indices.'.format(
                #         var_df[filter_condition][Vep.IMPACT_COL].count(), ds_name, len(select_indices))
            else:
                logger.error('No VEP results for case {}.'.format(ds_name))
                raise RuntimeError('No VEP results for case {}. '.format(ds_name))

        return var_df

    @staticmethod
    def generate_input_file(input_fp, data_df, chromosome_col='Chromosome', start_col='StartPosition',
                            end_col='EndPosition', reference_col='ReferenceAllele', alternate_col='AlternateAllele',
                            var_key_col=NT_VAR_COL, gene_col=None):
        """
        Generate input file for given nonsynonymous variants to run VEP based on pandas dataframe
        :param input_fp: path to file
        :param data_df: pandas dataframe with Chromosome, Position, Ref, and Alt columns
        :param chromosome_col: name of column with Chromosome information
        :param start_col: genomic start position of variant
        :param end_col: genomic end position of variant
        :param alternate_col: column name where alternate allele is stored
        :param reference_col: column name where reference allele is stored
        :param var_key_col:
        :param gene_col: optional name of column with gene name
        """

        if end_col is None:
            end_col = None
        else:
            end_col = data_df.columns.get_loc(end_col)

        if ((chromosome_col is None or start_col is None or reference_col is None or alternate_col is None)
                and var_key_col is not None):
            # extract information from variant key
            extract = True
            var_key_col = data_df.columns.get_loc(var_key_col)
            chr_col = pos_col = ref_col = alt_col = None
        else:
            chr_col = data_df.columns.get_loc(chromosome_col)
            pos_col = data_df.columns.get_loc(start_col)
            ref_col = data_df.columns.get_loc(reference_col)
            alt_col = data_df.columns.get_loc(alternate_col)
            var_key_col = None
            extract = False

        if gene_col is not None and gene_col in data_df.columns:
            gene_col = data_df.columns.get_loc(gene_col)
        else:
            gene_col = None

        with open(input_fp, 'w') as vep_file:
            tsv_vep = csv.writer(vep_file, delimiter='\t')
            for index, row in data_df.iterrows():
                if extract:
                    var_key = row[var_key_col]
                    chrom, pos, ref, alt = var_key.split('__')
                else:
                    chrom = row[chr_col]
                    pos = row[pos_col]
                    ref = row[ref_col]
                    alt = row[alt_col]

                if end_col is None:
                    end_pos = pos + ((len(ref) - 1) if alt == '-' or alt == '' else (len(alt) - 1))
                else:
                    end_pos = row[end_col]

                # write VEP input file
                tsv_vep.writerow([chrom, pos, end_pos, '{}/{}'.format(ref, alt),
                                  '' if gene_col is None else '# {}'.format(row[gene_col])])

            logger.info('Generated VEP input file from dataframe with {} entries: {}'.format(
                data_df.shape[0], input_fp))

    @staticmethod
    def run(input_fp, output_fp, reference_genome='hg19'):
        """
        Run tool VEP
        :param input_fp:
        :param output_fp:
        :param reference_genome: human reference genome, e.g. hg19 or hg20
        """

        if reference_genome == 'hg19':
            ref_genome = 'GRCh37'
        elif reference_genome == 'hg20':
            ref_genome = 'GRCh38'
        else:
            raise RuntimeError('VEP does not support reference genome {}.'.format(reference_genome))

        cmd = Vep.COMMAND.format(VEP_CACHE, ref_genome, os.path.abspath(input_fp), os.path.abspath(output_fp))
        logger.info('Call VEP with command {}'.format(cmd))

        with Popen(cmd, stdout=PIPE, cwd=Vep.PATH, shell=True) as p:
            for line in p.stdout:
                logger.info(line)
            exit_code = p.poll()

        if exit_code == 0 or exit_code is None:
            logger.info('Successfully called VEP: {}'.format(output_fp))
        else:
            logger.error('VEP call was unsuccessful. Is VEP installed?')
            raise RuntimeError('VEP quit with an error: {}'.format(exit_code))

    @staticmethod
    def read_results(output_fp):
        """
        Read VEP file and generate FATHMM input file
        :param output_fp:
        :return: read VEP variant dictionary
        """
        if os.path.isfile(output_fp):

            vep_df = pd.read_csv(output_fp, delimiter='\t', header=57)

            # create var_key column and use it as index column
            # var key is <chromosome>__<start position>__<reference allele>__<alternate allele>
            nt_var_keys = list()
            impacts = list()
            polyphen_preds = list()
            polyphen_scores = list()
            sift_preds = list()
            sift_scores = list()
            fathmm_ids = list()

            for index, row in vep_df.iterrows():
                chromosome, pos = row['Location'].split(':')
                ref = row['#Uploaded_variation'][row['#Uploaded_variation'].rfind('_')+1:].split('/')[0]
                if hasattr(row, 'Extra'):
                    if pos.find('-') != -1:
                        # identified insertion => update position calculation
                        if 'VARIANT_CLASS=insertion' in row['Extra']:
                            pos = str(pos[pos.find('-') + 1:])
                        else:
                            pos = pos[:pos.find('-')]
                else:
                    # extract start position of variant
                    pos = pos[:pos.find('-')]

                nt_var_key = '{}__{}__{}__{}'.format(chromosome, pos, ref, row['Allele'])
                nt_var_keys.append(nt_var_key)

                if hasattr(row, 'Extra'):
                    for info in row.Extra.split(';'):
                        if info.startswith('IMPACT='):
                            impacts.append(info[7:])
                            break
                    else:
                        impacts.append(NAN)

                    for info in row.Extra.split(';'):
                        if info.startswith('PolyPhen='):
                            if info.find('(') != -1 and info.find(')') != -1:
                                polyphen_preds.append(info[9:info.find('(')])
                                polyphen_scores.append(float(info[info.find('(') + 1:info.find(')')]))
                            else:
                                polyphen_preds.append(info[9:])
                                polyphen_scores.append(NAN)
                            break
                    else:
                        polyphen_preds.append(NAN)
                        polyphen_scores.append(NAN)

                    for info in row.Extra.split(';'):
                        if info.startswith('SIFT='):
                            if info.find('(') != -1 and info.find(')') != -1:
                                sift_preds.append(info[5:info.find('(')])
                                sift_scores.append(float(info[info.find('(') + 1:info.find(')')]))
                            else:
                                sift_preds.append(info[5:])
                                sift_scores.append(NAN)
                            break
                    else:
                        sift_preds.append(NAN)
                        sift_scores.append(NAN)

                    for info in row.Extra.split(';'):
                        # get protein ID and write it to FATHMM input file
                        if info.startswith('ENSP=') and 'missense_variant' in row.Consequence:
                            protein_id = info[5:]
                            amino_acids = row.Amino_acids.split('/')
                            substitution = '{}{}{}'.format(amino_acids[0], row['Protein_position'], amino_acids[1])
                            fathmm_ids.append(protein_id + '__' + substitution)
                            break
                    else:
                        fathmm_ids.append(NAN)

                else:
                    if row.PolyPhen != '-':
                        polyphen_preds.append(row.PolyPhen[:row.PolyPhen.find('(')])
                        polyphen_scores.append(float(row.PolyPhen[row.PolyPhen.find('(') + 1:row.PolyPhen.find(')')]))
                    else:
                        polyphen_preds.append(NAN)
                        polyphen_scores.append(NAN)

                    if row.SIFT != '-':
                        sift_preds.append(row.SIFT[:row.SIFT.find('(')])
                        sift_scores.append(float(row.SIFT[row.SIFT.find('(') + 1:row.SIFT.find(')')]))
                    else:
                        sift_preds.append(NAN)
                        sift_scores.append(NAN)

                    if row.ENSP != '-' and 'missense_variant' in row.Consequence:
                        protein_id = row.ENSP[5:]
                        amino_acids = row.Amino_acids.split('/')
                        substitution = '{}{}{}'.format(amino_acids[0], row['Protein_position'], amino_acids[1])
                        fathmm_ids.append(protein_id + '__' + substitution)
                    else:
                        fathmm_ids.append(NAN)

            vep_df[NT_VAR_COL] = pd.Series(nt_var_keys, index=vep_df.index)
            vep_df[FATHMM_KEY_COL] = pd.Series(fathmm_ids, index=vep_df.index)
            vep_df[PP_SCORE_COL] = pd.Series(polyphen_scores, index=vep_df.index)
            vep_df[SIFT_SCORE_COL] = pd.Series(sift_scores, index=vep_df.index)
            vep_df[Vep.IMPACT_COL] = pd.Series(impacts, index=vep_df.index)

            # check for duplicates
            if len(vep_df[vep_df.duplicated(NT_VAR_COL, keep=False)]) > 0:
                df = vep_df[vep_df.duplicated(NT_VAR_COL, keep=False)]
                # logger.error(df)
                logger.warning('Found {} duplicates in VEP output according to the {}!'.format(
                    len(df), NT_VAR_COL))

            vep_df.set_index(NT_VAR_COL, inplace=True)
            vep_df.drop_duplicates(keep='first', inplace=True)
            logger.info('Loaded {} VEP annotated variants from file {}.'.format(len(vep_df), output_fp))
            return vep_df

        else:
            logger.error('Missing VEP output file {}'.format(output_fp))
            raise RuntimeError('Missing VEP output file {}'.format(output_fp))
