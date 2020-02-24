""" class to read and process hotspot mutations from Chang et al, Cancer Discovery 2018 """
import logging
import numpy as np
import pandas as pd

from lifd.databases.database import Database
from lifd.utils import PT_VAR_COL

__author__ = 'Johannes Reiter'
__date__ = 'Jan 18, 2019'


# get logger
logger = logging.getLogger('lifd.{}'.format(__name__))

# column name in table file where variant amino acids are stored
VAA_COL = 'Variant Amino Acid'


class HotspotsDB(Database):

    NAME = 'Hotspots'

    # which columns are required to evaluate the database
    REQUIRED_COLS = [PT_VAR_COL]

    def __init__(self, db_source=None, threshold=0.1, weight=1.0, lazy_loading=True):
        """
        Constructor
        :param db_source: path to database file
        :param threshold: threshold such that database would predict functionality
        :param weight: support weight when database meets functionality prediction threshold
        :param lazy_loading: only load database if necessary for running LiFD
        """
        super().__init__(HotspotsDB.NAME, db_source, threshold, weight)

        if not lazy_loading:
            self.read_db()
        else:
            # Load database later if necessary for running LiFD
            self.db_df = None

    def read_db(self):

        if self.db_source is None:
            logger.error('Unable to load {} database without a provided source file!'.format(HotspotsDB.NAME))
            raise RuntimeError('No source file for {} database provided!'.format(HotspotsDB.NAME))
        else:
            logger.info('Reading database {} from file {}.'.format(self.name, self.db_source))

        self.db_df = pd.read_excel(self.db_source)

        self.db_df['Reference_Amino_Acid_#Samples'] = self.db_df.apply(
            lambda row: row.Reference_Amino_Acid.split(':')[1], axis=1)
        self.db_df['Reference_Amino_Acid'] = self.db_df.apply(
            lambda row: row.Reference_Amino_Acid.split(':')[0], axis=1)

        self.db_df['Variant_Amino_Acid_#Samples'] = self.db_df.apply(
            lambda row: row.Variant_Amino_Acid.split(':')[1], axis=1)
        self.db_df['Variant_Amino_Acid'] = self.db_df.apply(
            lambda row: row.Variant_Amino_Acid.split(':')[0], axis=1)

        self.db_df[PT_VAR_COL] = self.db_df.apply(
            lambda row: '{}__{}{}{}'.format(row['Hugo_Symbol'], row.Reference_Amino_Acid, row.Amino_Acid_Position,
                                            row.Variant_Amino_Acid), axis=1)
        self.db_df.set_index(PT_VAR_COL, inplace=True)

        n_rows = self.db_df.shape[0]
        self.db_df = self.db_df[self.db_df.qvalue_pancan <= self.threshold]
        logger.debug('Removed {} cancer hotspots because their q-value was above the threshold of {}.'.format(
            n_rows - self.db_df.shape[0], self.threshold))

        # self.db_df['HotspotKey'] = self.db_df.apply(
        #     lambda row: '{}__{}'.format(row['Hugo Symbol'], row.Codon), axis=1)
        # self.db_df.set_index('HotspotKey', inplace=True)

        logger.info('Loaded {} cancer hotspot variants from file {}.'.format(len(self.db_df), self.db_source))

    def in_database(self, nt_var_key=None, pt_var_key=None):
        if pt_var_key is None or pd.isnull(pt_var_key):
            # 'Hotspot variant check requires a PtVarKey (e.g. KRAS__G12V)')
            return np.nan

        if self.db_df is None:        # database loaded when needed
            self.read_db()

        if pt_var_key in self.db_df.index:
            return True
        elif '{}*'.format(pt_var_key[:-1]) in self.db_df.index:
            logger.info(f'Any codon change of {pt_var_key} was annotated as hotspot.')
            return True
        else:
            return False
