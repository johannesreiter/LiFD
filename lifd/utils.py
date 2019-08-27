""" Helper functions around LiFD """
import logging
import os
import math
import statsmodels.stats.multitest as ssm
import numpy as np
import pandas as pd

from lifd.settings import OUTPUT_DIR

# get logger
logger = logging.getLogger(__name__)

NAN = float('nan')

# column names for identifiers of nucleotide and protein variants
NT_VAR_COL = 'NtVarKey'
PT_VAR_COL = 'PtVarKey'
FATHMM_KEY_COL = 'FatHMM_Key'

PP_PRED_COL = 'PolyPhen_pred'
PP_SCORE_COL = 'PolyPhen'

SIFT_SCORE_COL = 'Sift'

FUNC_COL = 'MaybeFunctional'
MUT_EFFECT_COL = 'MutationEffect'

CT_COL = 'CancerType'


def init(output_dir=OUTPUT_DIR):
    """
    Create output directory for temporary files if not yet present
    :param output_dir: specific output directory has been provided as argument
    """

    output_directory = os.path.join(output_dir, '')

    # create directory for data files
    try:
        if not os.path.exists(os.path.dirname(output_directory)):
            os.makedirs(os.path.dirname(output_directory))
        logger.info('Output directory: {}'.format(os.path.abspath(output_directory)))
    except OSError:
        logger.info('Could not create output folder {} '.format(output_dir))

    return output_directory


def add_column(df, select_indices, df_key_col, new_col_name, dataset, as_value=True, sub_key=None, as_float=False):
    """
    Add column to a given dataframe from a given dataset (dictionary) or dataframe
    :param df: pandas dataframe
    :param select_indices:
    :param df_key_col: unique key of a variant either per subject or the whole dataset
    :param new_col_name: name of new column in original dataframe
    :param dataset: matrix with data or pandas dataframe
    :param as_value: if as_value than the value of the dataset will be added otherwise just True or False
    :param sub_key: sub key in provided dataset, e.g. list index
    :param as_float:
    """
    values = list()
    for index, row in df.iloc[select_indices, :].iterrows():

        key = row[df_key_col]

        if as_value:
            # check format of the provided dataframe
            if type(dataset) == pd.core.frame.DataFrame:

                # add float value
                if (key in dataset.index and sub_key is None
                        and dataset.loc[key] != '' and dataset.loc[key] != 'N/A'
                        and (not isinstance(dataset.loc[key], float) or not np.isnan(dataset.loc[key]))):

                    values.append(float(dataset.loc[key]) if as_float else dataset.loc[key])

                # add float value from sub_key
                elif ((key in dataset.index and sub_key in dataset.loc[key])
                        and ((isinstance(dataset.loc[key][sub_key], float) and not np.isnan(dataset.loc[key][sub_key]))
                             or (not isinstance(dataset.loc[key][sub_key], float) or dataset.loc[key][sub_key] != ''
                                 or dataset.loc[key][sub_key] != 'N/A'))):

                    values.append(float(dataset.loc[key][sub_key]) if as_float else dataset.loc[key][sub_key])

                # no value found
                else:
                    values.append(NAN)
            else:
                # add float value
                if (key in dataset and sub_key is None
                        and dataset[key] != '' and dataset[key] != 'N/A'):
                    values.append(float(dataset[key]) if as_float else dataset[key])
                # add float value from sub_key
                elif (key in dataset and (isinstance(sub_key, int) or (sub_key in dataset[key]))
                        and dataset[key][sub_key] != '' and dataset[key][sub_key] != 'N/A'):
                    values.append(float(dataset[key][sub_key]) if as_float else dataset[key][sub_key])
                # no value found
                else:
                    values.append(NAN)

        else:                                           # only puts True/False
            # check format of the provided dataframe
            if type(dataset) == pd.core.frame.DataFrame:
                if key in dataset.index and ((sub_key is None)
                                             or (sub_key is not None and sub_key in dataset.loc[key])):
                    values.append(True)
                else:
                    values.append(False)
            else:
                if key in dataset:
                    values.append(True)
                else:
                    values.append(False)

    # df[new_col_name] = pd.Series(values, index=df.index)
    if new_col_name not in df.columns:
        df[new_col_name] = NAN

    col_idx = df.columns.get_loc(new_col_name)
    df.iloc[select_indices, col_idx] = values
    logger.info('Added column {} with {} values for {} entries.'.format(
        new_col_name, sum(1 for val in values if (not isinstance(val, float)) or not np.isnan(val)),
        len(select_indices)))


def apply_multitest_correction(var_df, pval_col, fdr=0.1):
    """
    Use Benjamini-Hochberg method to perform multi test correction
    Adds new column to the dataframe with the extension '_corr
    :param var_df: pandas dataframe
    :param pval_col: name of column with the original P-values
    :param fdr: false discovery rate
    :return: extended pandas dataframe with corrected P-value column (pval_col + '_corr')
    """

    select_indices = np.where(var_df[FUNC_COL])[0]

    corr_col = pval_col + '_corr'
    pvals = np.array(var_df[var_df[pval_col].notnull()][pval_col].tolist())
    pvals[pvals == 0.0] = 1e-12
    rej, pval_corr = ssm.multipletests(pvals, alpha=fdr, method='fdr_bh',
                                       is_sorted=False, returnsorted=False)[:2]
    var_keys = var_df[var_df[pval_col].notnull()][NT_VAR_COL].tolist()
    pvals_corr = {key: val for key, val in zip(var_keys, pval_corr)}
    add_column(var_df, select_indices, NT_VAR_COL, corr_col, pvals_corr, sub_key=None, as_float=True)

    # calculate -log of corrected p-values
    var_df['mlog_' + corr_col] = var_df.apply(
        lambda row: math.log(row[corr_col]) * -1 if row[corr_col] != 0.0 else 6, axis=1)

    assert var_df[corr_col].count() > 0.3 * len(select_indices), \
        'Only {} corrected prediction found but {} indices.'.format(
            var_df[corr_col].count(), len(select_indices))

    return var_df
