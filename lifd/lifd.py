""" Methods to predict the functionality of individual mutations """
import logging
import os
import numpy as np
import pandas as pd
import datetime

from lifd.predictors.predictor import Predictor
from lifd.predictors.cravat import Cravat
from lifd.utils import apply_multitest_correction, NT_VAR_COL, PT_VAR_COL, FUNC_COL, MUT_EFFECT_COL, CT_COL
from lifd.settings import DR_LIST_FP
from lifd.version import __version__

# get logger
logger = logging.getLogger(__name__)

try:    # check if varcode and pyensembl is available (necessary for Windows)
    from varcode import Variant  # https://github.com/hammerlab/varcode
    # from pyensembl import ensembl_grch37
    from lifd.mutation_effects import get_mutation_details, annotate_effects, MutationEffects, VAR_EFFECT_COLS
    logger.info('VarCode is used for mutation effect prediction.')

except ImportError:
    # mutation effect prediction can not be performed since VarCode is not avilable
    annotate_effects = None
    MutationEffects = None
    logger.error('VarCode could not be imported! Mutation effect prediction can not be performed automatically.')

__author__ = 'Johannes REITER'
__date__ = 'Jan 19, 2019'


LIFD_COL = 'LiFD'
LIFD_SUP_COL = 'LiFD_support'
LIFD1_SUP_COL = 'LiFD1_support'
LIFD2_SUP_COL = 'LiFD2_support'
LIFD2_SUP_TUPLE_COL = 'LiFD2_results'
# DR_GENE_COL = 'TCGA_driver_gene'

# require that at least one database or the majority of predictors/tools predict functionality
MIN_SUPPORT_TH = 0.5001


class LiFD:

    def __init__(self, databases, predictors, min_support_th=MIN_SUPPORT_TH,
                 driver_gene_fp=DR_LIST_FP, driver_gene_col='DriGeneClf', reference_genome='hg19'):
        """
        Initiate LiFD framework with databases and predictors for analysis
        :param databases:
        :param predictors:
        :param min_support_th:
        :param driver_gene_fp:
        :param driver_gene_col:
        :param reference_genome: human reference genome, e.g. hg19 or hg20
        """

        # provided mutation coordinates are for given reference genome
        self.ref_genome = reference_genome
        # set reference genome for mutation effect analysis
        if annotate_effects is not None:
            MutationEffects.set_ReferenceGenome(reference_genome=self.ref_genome)

        # databases for LiFD1
        if databases is not None and len(databases):
            self._databases = databases
            self.database_support_contribution = 1
            logger.info('Using databases: {}'.format(', '.join(db.name for db in self.databases)))

        else:
            self._databases = list()
            logger.warning('No databases provided for LiFD phase-1 predictions.')

        # predictors for LiFD2
        if predictors is not None and len(predictors):
            self._predictors = predictors
            self.predictor_support_contribution = 1.0 / len(self.predictors)
            for prd in predictors:
                if prd.SUPP_REF_GENOMES is not None and reference_genome not in prd.SUPP_REF_GENOMES:
                    raise RuntimeError('{} does not support reference genome {}.'.format(
                        prd.__name__, reference_genome))
            logger.info('Using predictors: {}'.format(', '.join(prd.__name__ for prd in self.predictors)))

        else:
            self._predictors = list()
            logger.warning('No predictors provided for LiFD phase-2 predictions.')

        # Minimal required support such that a variant is annotated as likely functional
        self.min_support_th = min_support_th

        # define boolean column whether a nonsynonymous mutation occurred in a putative driver gene
        self.driver_gene_clf_col = driver_gene_col

        if driver_gene_fp is None:
            self.driver_gene_fp = None
            self.driver_genes = None
            logger.info('No file with driver genes was given. Using column with driver gene classification: {}'
                        .format(driver_gene_col))
        else:
            # read given driver gene list, default TCGA consensus driver list
            self.driver_gene_fp = driver_gene_fp
            dr_df = pd.read_csv(driver_gene_fp)
            self.driver_genes = set(dr_df.Gene_Symbol.unique())
            logger.info('Using {} driver genes for annotation from {}.'.format(len(self.driver_genes), DR_LIST_FP))

        logger.info('Initialized LiFD {} with {} databases and {} predictors.'.format(
            __version__, len(self.databases), len(self.predictors)))

    @property
    def databases(self):
        return self._databases

    @databases.setter
    def databases(self, dbs):
        self._databases = dbs

    @property
    def predictors(self):
        return self._predictors

    @predictors.setter
    def predictors(self, pdrs):

        if len(pdrs) == 0:
            logger.warning('No predictors provided for {} instance.'.format(self.__name__))

        for p in pdrs:
            if not isinstance(p, Predictor):
                raise AttributeError('Passed object {} has to be an instance of {}!'.format(p, self.__name__))

        self._predictors = pdrs

    def check_required_cols(self, var_df):

        if 'Chromosome' in var_df.columns.values:
            var_df.Chromosome = var_df.Chromosome.astype(str)

        db_req_cols = [db.REQUIRED_COLS for db in self.databases if db.pred_col not in var_df.columns.values]
        prd_req_cols = [prd.REQUIRED_COLS for prd in self.predictors
                        if any(pred_col not in var_df.columns.values for pred_col in prd.PREDICTION_COLS)]

        req_cols_list = db_req_cols + prd_req_cols
        req_cols = set().union(*req_cols_list)
        req_cols.add(FUNC_COL)

        # perform basic mutation effect analysis using VARCODE if it was not performed previously
        for req_col in req_cols:
            if req_col not in var_df.columns.values:
                if req_col in VAR_EFFECT_COLS:
                    if annotate_effects is None:
                        logger.error(
                            '{} column is missing and cannot be added because VARCODE could not be invoked!'.format(
                                FUNC_COL))
                        raise RuntimeError('{} column is missing!'.format(req_col))
                    else:
                        logger.debug('Missing {} column. '.format(req_col))
                        var_df = annotate_effects(var_df)

                else:
                    required_for_dbs = list()
                    for db in self.databases:
                        if req_col in db.REQUIRED_COLS and db.pred_col not in var_df.columns.values:
                            required_for_dbs.append(db.name)
                    if len(required_for_dbs) > 0:
                        raise RuntimeError('Column {} is missing and required to evaluate database {}!'.format(
                            req_col, ', '.join(required_for_dbs)))

                    required_for_prds = list()
                    for prd in self.predictors:
                        if req_col in prd.REQUIRED_COLS and \
                                any(pred_col not in var_df.columns.values for pred_col in prd.PREDICTION_COLS):
                            required_for_prds.append(prd.__name__)
                    if len(required_for_prds) > 0:
                        raise RuntimeError('Column {} is missing and required to invoke predictor {}!'.format(
                            req_col, ', '.join(required_for_prds)))

        return var_df

    def run_lifd(self, var_df=None, output_dir=None, export_fn=None):
        """
        Run two-phase algorithm LiFD for the provided variants and extend the given dataframe with the resulting columns
        :param var_df: pandas dataframe with all required columns for the used databases and predictors
        :param output_dir: directory where output of the algorithm will be stored
        :param export_fn: save the extended dataframe to the given excel filename
        :return: extended dataframe
        """

        var_df = self.check_required_cols(var_df)

        if self.driver_genes is None:
            if self.driver_gene_clf_col not in var_df.columns.values:
                logger.error('{} column is missing '.format(self.driver_gene_clf_col)
                             + 'and cannot be generated because no driver gene list was defined in the settings file!')
                raise RuntimeError('{} boolean column is missing!'.format(self.driver_gene_clf_col))

        else:
            # define variants in driver genes that could possibly have an effect
            var_df.loc[(var_df.GeneSymbol.isin(self.driver_genes)) & var_df[FUNC_COL],
                       self.driver_gene_clf_col] = True
            var_df.loc[~((var_df.GeneSymbol.isin(self.driver_genes)) & var_df[FUNC_COL]),
                       self.driver_gene_clf_col] = False
            logger.info('Added potentially functional driver gene mutation annotation to {} variants.'.format(
                len(var_df[var_df[self.driver_gene_clf_col]])))

        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
            logger.info('Created LIFD output directory: {}'.format(os.path.abspath(output_dir)))
        else:
            logger.info('LIFD output directory: {}'.format(os.path.abspath(output_dir)))

        # original columns in the dataframe
        original_cols = var_df.columns

        # ##################### FIRST PHASE ####################
        if len(self.databases) > 0:
            var_df = self.run_lifd1(var_df)
        else:
            var_df[LIFD1_SUP_COL] = pd.Series(list(np.zeros((len(var_df.index)))), index=var_df.index)
        # ######################################################

        # ##################### SECOND PHASE ###################
        if len(self.predictors) > 0:
            var_df = self.run_lifd2(var_df, output_dir=output_dir)
        else:
            var_df[LIFD2_SUP_TUPLE_COL] = pd.Series(list(np.zeros((len(var_df.index), 2))), index=var_df.index)
            var_df[LIFD2_SUP_COL] = pd.Series(list(np.zeros((len(var_df.index)))), index=var_df.index)
        # ######################################################

        # calculate total LiFD support
        var_df[LIFD_SUP_COL] = var_df.apply(
            lambda row: row[LIFD1_SUP_COL] + row[LIFD2_SUP_COL],  # phase 1 + phase 2
            axis=1)

        # is variant in a driver gene and possibly functional and has sufficient LiFD support?
        var_df[LIFD_COL] = var_df.apply(
            lambda row: row[self.driver_gene_clf_col] and row[FUNC_COL] and (row[LIFD_SUP_COL] >= MIN_SUPPORT_TH),
            axis=1)

        logger.info('LiFD predicts that {} of {} ({:.1%})'.format(
            len(var_df[var_df[LIFD_COL]]), len(var_df[var_df[self.driver_gene_clf_col] & var_df[FUNC_COL]]),
            len(var_df[var_df[LIFD_COL]]) / len(var_df[var_df[self.driver_gene_clf_col] & var_df[FUNC_COL]]))
                    + ' nonsynonymous or splice-site variants in driver genes are functional.')
        logger.info('Positive predictions of LiFD1 {}/{} ({:.1%})'.format(
            len(var_df[(var_df[LIFD1_SUP_COL] >= MIN_SUPPORT_TH) & var_df[self.driver_gene_clf_col]
                       & var_df[FUNC_COL]]),
            len(var_df[var_df[self.driver_gene_clf_col] & var_df[FUNC_COL]]),
            len(var_df[(var_df[LIFD1_SUP_COL] >= MIN_SUPPORT_TH) & var_df[self.driver_gene_clf_col] & var_df[FUNC_COL]])
            / len(var_df[var_df[self.driver_gene_clf_col] & var_df[FUNC_COL]])))
        logger.info('Positive predictions of LiFD2 {}/{} ({:.1%})'.format(
            len(var_df[(var_df[LIFD2_SUP_COL] >= MIN_SUPPORT_TH) & var_df[self.driver_gene_clf_col]
                       & var_df[FUNC_COL]]),
            len(var_df[var_df[self.driver_gene_clf_col] & var_df[FUNC_COL]]),
            len(var_df[(var_df[LIFD2_SUP_COL] >= MIN_SUPPORT_TH) & var_df[self.driver_gene_clf_col] & var_df[FUNC_COL]])
            / len(var_df[var_df[self.driver_gene_clf_col] & var_df[FUNC_COL]])))

        if export_fn is not None:
            if export_fn.endswith('.xlsx'):
                output_fp = os.path.join(output_dir, export_fn)
                writer = pd.ExcelWriter(output_fp, engine='xlsxwriter')
                # add a bold format to use to highlight cells
                wb = writer.book
                bold = wb.add_format({'bold': True})

                # reorder output columns such that the most informative columns come first
                output_cols = list()
                important_cols = ['Subject', CT_COL, 'GeneSymbol', LIFD_COL, LIFD_SUP_COL]
                for c in important_cols:
                    if c in var_df.columns:
                        output_cols.append(c)

                n_freeze_cols = len(output_cols)

                genomic_cols = ['Chromosome', 'StartPosition', 'EndPosition', 'ReferenceAllele',
                                'AlternateAllele', MUT_EFFECT_COL, FUNC_COL]
                if all(c in var_df.columns.values for c in genomic_cols):
                    output_cols = output_cols + genomic_cols

                for c in var_df.columns.values:
                    if c not in output_cols:
                        output_cols.append(c)

                caption = [
                    f'Found {len(var_df[var_df[LIFD_COL]])} LiFDs in {len(var_df)} provided variants.',
                    f'{self.driver_gene_clf_col}: is mutation in the given driver list {self.driver_gene_fp}',
                    f'{LIFD_COL}: is mutation predicted to be a likely functional driver gene mutation',
                    f'{LIFD_SUP_COL}: >=1 if mutation is a known oncogenic variant or otherwise the '
                    + 'fraction of functionality supporting tools',
                    f'{LIFD1_SUP_COL}: number of databases listing variant as oncogenic',
                    f'{LIFD2_SUP_TUPLE_COL}: tuple: first entry denotes number of tools predicting driver mutation'
                    + 'second entry denotes the number of tools that produced a valid prediction',
                    ]
                var_df[output_cols].to_excel(
                    writer, index=False, freeze_panes=(len(caption) + 3, n_freeze_cols), startrow=len(caption) + 2,
                    startcol=0)

                worksheet = writer.sheets['Sheet1']
                # increase column width for better readability
                worksheet.set_column(1, 12, 13)
                worksheet.set_column(13, len(output_cols)-1, 11)
                now = datetime.datetime.now()
                header = f'LiFD {__version__} output generated ' \
                         + f'at {now:%H:%M:%S} on {now:%Y-%m-%d}.'
                worksheet.write(0, 0, header, bold)
                for line_no, line_text in enumerate(caption, 1):
                    worksheet.write(line_no, 0, line_text)

                writer.save()
                writer.close()
                logger.info('Exported {} annotated variants ({} LiFD) to {}.'.format(
                    len(var_df), len(var_df[var_df[LIFD_COL]]), output_fp))
            else:
                logger.error('Could not export LiFD results because given filename did not end with .xlsx: '
                             + str(export_fn))

        return var_df

    def run_lifd1(self, var_df=None, nt_var_key=None, pt_var_key=None):
        """
        Run first phase of LiFD by going through all provided databases and checking whether the variant is present
        If a dataframe with variants is provided, presence columns will be added for each database
        :param var_df: pandas dataframe with variants and the required columns to evaluate the databases
        :param nt_var_key: nucleotide variant key which enables the presence check across databases
        :param pt_var_key: protein variant key which enables the presence check across databases
        :return: either extended dataframe with LiFD1 results or LiFD1 support value
        """

        # run LiFD across whole dataframe
        if var_df is not None and isinstance(var_df, pd.DataFrame):

            for db in self.databases:
                # is database annotation already provided?
                if db.pred_col not in var_df.columns.values:
                    # check if all necessary columns are provided
                    if not all(col in var_df.columns for col in db.REQUIRED_COLS):
                        logger.error('{} database requires these columns: '.format(db.name)
                                     + ', '.join(c for c in db.REQUIRED_COLS))
                        logger.error(
                            'Only these columns were provided: ' + ', '.join(str(col) for col in var_df.columns))
                        raise RuntimeError('Not all required columns are given for database {}: '.format(db.name)
                                           + ', '.join(str(col) for col in db.REQUIRED_COLS))

                    var_df[db.pred_col] = var_df.apply(
                        lambda row: db.in_database(row[NT_VAR_COL], row[PT_VAR_COL]), axis=1)

            # calculate support from databases
            var_df[LIFD1_SUP_COL] = var_df.apply(
                lambda row: sum(db.weight for db in self.databases if row[db.pred_col] >= db.threshold), axis=1)

            return var_df

        elif nt_var_key is not None and pt_var_key is not None:
            # run LiFD1 for a given variant
            raise NotImplementedError
            # lifd1_support = 0.0
            # for db in self.databases:
            #     if db.in_database(nt_var_key, pt_var_key):
            #         lifd1_support += self.database_support_contribution
            #
            # return lifd1_support

        else:
            logger.warning('Either pandas dataframe or variant keys need to be provided to run first phase of LiFD.')

    def run_lifd2(self, var_df=None, nt_var_key=None, pt_var_key=None, output_dir=None):
        """
        Run second phase of LiFD by going through all predictors and evaluating their predictions
        If a dataframe with variants is provided, predictor columns will be added for each predictor
        :param var_df: pandas dataframe with variants and the required columns to evaluate the databases
        :param nt_var_key: nucleotide variant key which enables the presence check across databases
        :param pt_var_key: protein variant key which enables the presence check across databases
        :param output_dir:
        :return: either extended dataframe with LiFD2 results or LiFD2 support value
        """

        # run LiFD2 across whole dataframe
        if var_df is not None and isinstance(var_df, pd.DataFrame):

            for prd in self.predictors:
                # is predictor annotation already provided?
                if any(pred_col not in var_df.columns.values for pred_col in prd.PREDICTION_COLS):
                    # run predictor
                    # check if all necessary columns are provided
                    if not all(col in var_df.columns for col in prd.REQUIRED_COLS):
                        logger.error('{} predictor requires these columns: '.format(prd.__name__)
                                     + ', '.join(c for c in prd.REQUIRED_COLS))
                        logger.error(
                            'Only these columns were provided: ' + ', '.join(str(col) for col in var_df.columns))
                        raise RuntimeError('Not all required columns are given for predictor {}: '.format(prd.__name__)
                                           + ', '.join(str(col) for col in prd.REQUIRED_COLS))

                    for ct in var_df.CancerType.unique():
                        var_df = prd.get_annotation(var_df, ds_name=ct, output_dir=output_dir,
                                                    filter_condition=(var_df.CancerType == ct),
                                                    reference_genome=self.ref_genome)

                    # apply multitest correction for CHASMplus output
                    if prd == Cravat:
                        var_df = apply_multitest_correction(var_df, Cravat.CP_COL)
                        var_df = apply_multitest_correction(var_df, Cravat.CP_CT_COL)

            # calculate support from predictors
            # returns a tuple with the No. tools predicting functionality and No. tools producing a result
            lifd2_supports = list()
            for _, row in var_df.iterrows():
                lifd2_supports.append(np.sum([prd.predict_functionality(row) for prd in self.predictors], axis=0))
            var_df[LIFD2_SUP_TUPLE_COL] = pd.Series(lifd2_supports, index=var_df.index)

            var_df[LIFD2_SUP_COL] = var_df.apply(
                lambda row: 0 if row[LIFD2_SUP_TUPLE_COL][1] == 0 else (
                        row[LIFD2_SUP_TUPLE_COL][0] / row[LIFD2_SUP_TUPLE_COL][1]), axis=1)

            # apply does not work for tuples
            # var_df[LIFD2_SUP_COL] = var_df.apply(
            #     lambda row: np.sum([prd.predict_functionality(row) for prd in self.predictors], axis=0),
            #     axis=1)

            return var_df

        elif nt_var_key is not None and pt_var_key is not None:
            raise NotImplementedError

        else:
            logger.warning('Either pandas dataframe or variant keys need to be provided to run second phase of LiFD.')
