""" class to read and process CGI oncogenic mutations from Tamborero et al, Genome Medicine 2018 """
import logging
import pandas as pd
import json
import os
import pysam

from lifd.databases.database import Database
from lifd.utils import PT_VAR_COL
from lifd.settings import DB_DIR

__author__ = 'Johannes Reiter'
__date__ = 'Jan 19, 2019'


# get logger
logger = logging.getLogger('lifd.{}'.format(__name__))

# column name in TSV file where variant amino acids are stored
VAA_COL = 'Variant Amino Acid'


class CgiDB(Database):

    NAME = 'CGI_Catalog'

    # which columns are required to evaluate the database
    REQUIRED_COLS = [PT_VAR_COL]

    def __init__(self, db_source, threshold=True, weight=1.0, lazy_loading=True):
        """
        Constructor
        :param db_source: path to database file
        :param threshold: threshold such that database would predict functionality
        :param weight: support weight when database meets functionality prediction threshold
        :param lazy_loading: only load database if necessary for running LiFD
        """
        super().__init__(CgiDB.NAME, db_source, threshold, weight)
        
        self.processed_db_fp = None
        
        if not lazy_loading:
            self.read_db()
        else:
            # Load database later if necessary for running LiFD
            self.db_df = None

    def read_db(self):

        if self.db_source is None:
            logger.error('Unable to load {} database without a provided source file!'.format(CgiDB.NAME))
            raise RuntimeError('No source file for {} database provided!'.format(CgiDB.NAME))
        else:
            logger.info('Reading database {} from file {}.'.format(self.name, self.db_source))

        # self.db_df = pd.read_csv(self.db_source, delimiter='\t', encoding='Latin-1')
        # self.db_df['PtVarKey'] = self.db_df.apply(
        #     lambda row: '{}__{}'.format(row.gene, row.protein[2:]), axis=1)
        # self.db_df.set_index('PtVarKey', inplace=True)
        # self.db_df = self.db_df.loc[~self.db_df.index.duplicated(keep='first')]

        # code
        if isinstance(self.db_source, str):
            ext = os.path.splitext(self.db_source)[1]
            self.processed_db_fp = self.db_source.replace(ext, '_processed.json')

        # implemented by Kam Louie; missing additional comments
        if self.processed_db_fp is not None and os.path.exists(self.processed_db_fp) \
                and os.path.isfile(self.processed_db_fp):

            with open(self.processed_db_fp, 'r') as f:
                self.db_df = json.load(f)

            logger.info('Read previously processed CGI database with {} distinct mutations.'.format(len(self.db_df)))

        else:
            tmp = pd.read_csv(self.db_source, delimiter='\t', encoding='Latin-1')
            self.db_df = list()
            for entries in tmp['gdna']:
                entries = entries.split('__')
                for entry in entries:
                    # print(entry)
                    chrom, gen = entry.split(':')
                    chrom = chrom[3:]
                    index = None
                    for i, char in enumerate(gen[2:]):
                        if char.isalpha():
                            index = i
                            break
                    gloc = None
                    gcoords = gen[2:i+2].split('_')
                    if 'dup' in gen and len(gcoords) == 2:
                        gloc = gcoords[1]
                    else:
                        gloc = gcoords[0]
                    mut = gen[i+2:]
                    ref_a = None
                    alt_a = None
                    if '>' in mut:
                        sep_idx = mut.find('>')
                        ref_a = mut[:sep_idx]
                        alt_a = mut[sep_idx+1:]
                    if 'dup' in mut:
                        ref_a = '-'
                        alt_a = mut[3:]
                    if 'del' in mut and 'ins' in mut:
                        del_idx = mut.find('del')
                        ins_idx = mut.find('ins')
                        ref_a = mut[del_idx+3:ins_idx]
                        alt_a = mut[ins_idx+3:]
                        if ref_a is '':
                            start = int(gcoords[0])
                            end = None
                            if len(gcoords) == 1:
                                end = int(gcoords[0])
                            else:
                                end = int(gcoords[1])
                            genome = pysam.Fastafile(os.path.join(DB_DIR, 'hg19.fa'))
                            ref_a = genome.fetch('chr' + chrom, start-1, end)
                    elif 'del' in mut:
                        ref_a = mut[3:]
                        alt_a = '-'
                    elif 'ins' in mut:
                        ref_a = '-'
                        alt_a = mut[3:]

                    nt_var_key = '{}__{}__{}__{}'.format(chrom, gloc, ref_a, alt_a)
                    self.db_df.append(nt_var_key)

            if self.processed_db_fp is not None:
                with open(self.processed_db_fp, 'w') as f:
                    json.dump(self.db_df, f)

        logger.info('Loaded {} CGI oncogenic variants from file {}.'.format(len(self.db_df), self.db_source))

    def in_database(self, nt_var_key=None, pt_var_key=None):
        if nt_var_key is None:
            raise AttributeError('CGI oncogenic variant check requires a NtVarKey.')

        if self.db_df is None:        # database loaded when needed
            self.read_db()

        # var_contents = nt_var_key.split('__') # [chr number, start genome position, ref allele, alt allele]
        # conv_nt_var_key = 'chr' + var_contents[0] + ':g.' + START + _ + END + MUTATION
        # genome = pysam.Fastafile(os.path.join(DB_DIR, 'hg19.fa'))
        # sequence = genome.fetch(chr, start, end)
        if nt_var_key in self.db_df:
            return True
        else:
            return False
