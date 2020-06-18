""" class to read and process variants from genome-wide screen in COSMIC Tate et al., Nucleic Acids Res 2019 """
import logging
import os
import numpy as np
import pandas as pd
from collections import Counter
import pysam

from lifd.utils import NT_VAR_COL, PT_VAR_COL
from lifd.databases.database import Database
from lifd.settings import REF_GENOME_FA_FP

__author__ = 'Johannes Reiter'
__date__ = 'Jan 19, 2019'


# get logger
logger = logging.getLogger('lifd.{}'.format(__name__))

# get nucleotide complement for reverse strand
NT_COMPLEMENT = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}


class CosmicDB(Database):

    NAME = 'COSMIC'

    # which columns are required to evaluate the database
    REQUIRED_COLS = [NT_VAR_COL]

    # default q-value significance threshold
    OCCURRENCES_TH = 4

    def __init__(self, db_source=None, threshold=OCCURRENCES_TH, weight=1.0, lazy_loading=True):
        """
        Constructor
        :param db_source: path to database file
        :param threshold: threshold such that database would predict functionality
        :param weight: support weight when database meets functionality prediction threshold
        :param lazy_loading: only load database if necessary for running LiFD
        """

        if not os.path.isfile(db_source):
            logger.error(f'COSMIC file is not available at: {db_source}')
            logger.error(f'Adapt path to COSMIC file in lifd/settings.py')
            raise RuntimeError(f'COSMIC file is not available at: {db_source}')

        super().__init__(CosmicDB.NAME, db_source, threshold, weight)

        self.processed_db_fp = None

        if lazy_loading:
            # Load database later if necessary for running LiFD
            self.db_df = None
        else:
            self.read_db()

    def read_db(self):

        if self.db_source is None:
            logger.error('Unable to load {} database without a provided source file!'.format(CosmicDB.NAME))
            raise RuntimeError('No source file for {} database provided!'.format(CosmicDB.NAME))
        else:
            logger.info('Reading database {} from file {}.'.format(self.name, self.db_source))

        # check if COSMIC database was processed previously
        # if yes, simply read the processed CSV file
        if isinstance(self.db_source, str):
            ext = os.path.splitext(self.db_source)[1]
            self.processed_db_fp = self.db_source.replace(ext, '_processed.csv')

        if self.processed_db_fp is not None and os.path.exists(self.processed_db_fp) \
                and os.path.isfile(self.processed_db_fp):

            self.db_df = pd.read_csv(self.processed_db_fp, index_col=NT_VAR_COL)

            logger.info('Read previously processed genome/exome-wide sequencing data of '
                        'COSMIC samples with in total {} distinct mutations.'.format(sum(self.db_df.Occurrences)))

        # extract necessary information and save it to CSV file
        else:
            cosmic_df = pd.read_csv(
                self.db_source, delimiter='\t',
                usecols=['ID_sample', 'Mutation ID', 'Mutation CDS', 'Mutation genome position', 'Mutation strand',
                         'Gene name', 'Mutation AA', 'Tumour origin', 'GRCh'])  # only <= v89 is supported
            # HGVS could be useful to get the correct variant change despite the inconsistent mutation strand
            # usecols=['ID_sample', 'MUTATION_ID', 'Mutation CDS', 'Mutation genome position', 'Mutation strand',
            #          'Gene name', 'Mutation AA', 'Tumour origin', 'GRCh'])  # v90 is not supported

            cosmic_df.set_index(['ID_sample', 'Mutation ID'], inplace=True)

            logger.info('Finished reading COSMIC input file {}. Now processing...'.format(self.db_source))

            # nucleotides
            # nts = 'ACGT'
            cosmic_nt_vars = Counter()
            # cosmic_pt_vars = Counter()
            # pt_nt_versions = defaultdict(set)
            # nt_pt_versions = dict()
            gene_names = dict()
            cosmic_samples = set()
            # examples:
            # Mutation CDS                     c.5119_5120GA>TC
            # Mutation AA                             p.E1707>?
            # Mutation genome position    2:168103021-168103022
            liftover = None
            for (sample_id, mutation_id), row in cosmic_df.iterrows():
                # no genomic information available => skip this variant
                mut_gen_pos = row['Mutation genome position']
                if isinstance(mut_gen_pos, float) and np.isnan(mut_gen_pos):
                    continue

                # create Nucleotide Variant Key
                chrom, pos = mut_gen_pos.split(':')
                start_pos = None
                end_pos = None
                if '-' in pos:
                    start_pos = int(pos.split('-')[0])   # start position
                    end_pos = int(pos.split('-')[1])     # end position
                else:
                    start_pos = int(pos)
                    end_pos = int(pos)

                mut_cds = row['Mutation CDS']
                if mut_cds.find('>') != -1:
                    # find start of reference allele
                    refs, alt = mut_cds.split('>')
                    if refs.find('_') != -1:
                        refs = refs[refs.find('_')+1:]
                        ref = ''.join(i for i in refs if not i.isdigit() and i != '+' and i != '-')
                    else:
                        ref = refs[-1]

                    if ref.isdigit():
                        ref = 'del' + ref
                    if alt.isdigit():
                        alt = 'ins' + alt
                    # change_idxs = [ref.find(nt) for nt in nts if ref.find(nt) != -1]
                    # if len(change_idxs) > 0:
                    #     change_idx = min(change_idxs)
                    #     assert change_idx > 0, 'Could not find {} in {}'.format(nts, ref)
                    #     # ref = ref[change_idx:]
                    # else:
                    #     # logger.debug('Change not provided => skip: {}'.format(info))
                    #     alt = None

                # was variant a deletion-insertion?
                elif mut_cds.lower().find('delins') != -1:
                    try:
                        change_idx = mut_cds.lower().find('delins') + 6
                        genome = pysam.Fastafile(REF_GENOME_FA_FP)
                        ref = genome.fetch('chr' + ('Y' if chrom == '23' else chrom), start_pos-1, end_pos).upper()
                        if ref.isdigit() or ref == '?':
                            ref = 'del' + ref
                        alt = mut_cds[change_idx:]
                        if alt.isdigit() or alt == '?':
                            alt = 'del' + alt
                    except KeyError as err:
                        logger.debug(row)
                        logger.warning('PySam: {}'.format(err))
                        continue

                # was variant a deletion?
                elif mut_cds.lower().find('del') != -1:
                    change_idx = mut_cds.lower().find('del') + 3
                    ref = mut_cds[change_idx:]
                    if ref.isdigit() or ref == '?':
                        ref = 'del' + ref
                    alt = '-'

                # was variant an insertion?
                elif mut_cds.lower().find('ins') != -1:
                    change_idx = mut_cds.lower().find('ins') + 3
                    alt = mut_cds[change_idx:]
                    if alt.isdigit() or alt == '?':
                        alt = 'ins' + alt
                    ref = '-'
                else:
                    # logger.debug('Change not provided in COSMIC => skip: {}'.format(row))
                    continue

                if row['Mutation strand'] == '-':
                    if len(ref) == 1 and ref in NT_COMPLEMENT.keys():
                        ref = NT_COMPLEMENT[ref]
                    elif all(nt in NT_COMPLEMENT.keys() for nt in ref):
                        complement_ref = ref[::-1]
                        ref = ''
                        for nt in complement_ref:
                            ref += NT_COMPLEMENT[nt]

                    if len(alt) == 1 and alt in NT_COMPLEMENT.keys():
                        alt = NT_COMPLEMENT[alt]
                    elif all(nt in NT_COMPLEMENT.keys() for nt in alt):
                        complement_alt = alt[::-1]
                        alt = ''
                        for nt in complement_alt:
                            alt += NT_COMPLEMENT[nt]

                if '_ENST' in row['Gene name']:
                    gene_name = row['Gene name'][:row['Gene name'].find('_ENST')]
                else:
                    gene_name = row['Gene name']

                # found information for unique genomic identifier
                nt_key = f'{chrom}__{start_pos}__{ref}__{alt}'
                cosmic_nt_vars[nt_key] += 1
                cosmic_samples.add(sample_id)
                gene_names[nt_key] = gene_name

                if not (all(nt in NT_COMPLEMENT.keys() for nt in ref) or ref == '-' or 'del' in ref.lower()
                        or 'ins' in alt.lower() or ref.isdigit() or alt.isdigit()):

                    logger.warning('Error in formatting of reference {} or alternate allele {}'.format(ref, alt))
                    logger.debug('Row: {}'.format(row))
                    #raise RuntimeError('COSMIC formatting error: {}'.format(row))

                # for debugging of correct formatting of variants only
                # elif not ('del' in ref.lower() or 'ins' in alt.lower()):
                #     try:
                #         v = Variant(contig=chrom, start=pos, ref=ref, alt=alt, ensembl=ensembl_grch37)
                #         top_priority_effect = v.effects().top_priority_effect()
                #     except Exception as err:
                #         print(ref)
                #         print(alt)
                #         print(row)
                #         print(traceback.format_exc())
                #         # raise err

                # else:
                #     nt_key = None
                # create Protein Variant Key
                # if row['Mutation AA'].startswith('p.'):
                #     pt_key = '{}__{}'.format(gene_name, row['Mutation AA'][2:])
                #     cosmic_pt_vars[pt_key] += 1
                #     if nt_key is not None and nt_key is not None:
                #         pt_nt_versions[pt_key].add(nt_key)
                #         if nt_key in nt_pt_versions and pt_key != nt_pt_versions[nt_key]:
                #             logger.info('Variant {} in {} has been mapped to {} but also exists as {}!'.format(
                #                 nt_key, gene_name, nt_pt_versions[nt_key], pt_key))
                #         else:
                #             nt_pt_versions[nt_key] = pt_key

            self.db_df = pd.DataFrame(list(cosmic_nt_vars.items()), columns=[NT_VAR_COL, 'Occurrences'])
            self.db_df.set_index(NT_VAR_COL, inplace=True)
            # self.db_df[PT_VAR_COL] = pd.Series(nt_pt_versions)
            self.db_df['GeneSymbol'] = pd.Series(gene_names)

            self.db_df.to_csv(self.processed_db_fp)
            
            logger.info('Processed genome/exome-wide sequencing data of '
                        '{} COSMIC samples with in total {} distinct mutations.'.format(
                         len(cosmic_samples), sum(cosmic_nt_vars.values())))

    def in_database(self, nt_var_key=None, pt_var_key=None):
        if nt_var_key is None and pt_var_key is None:
            raise AttributeError(f'COSMIC variant check requires either the {NT_VAR_COL} (e.g., 12__25398285__C)'
                                 + f' or the {PT_VAR_COL} (e.g. KRAS__G12V)!')

        if self.db_df is None:        # database loaded when needed
            self.read_db()

        # nucleotide variant gets priority over protein variant
        if nt_var_key is not None:
            if nt_var_key in self.db_df.index:
                return self.db_df.loc[nt_var_key]['Occurrences']
            else:
                return 0
        else:
            if PT_VAR_COL not in self.db_df.columns:
                raise RuntimeError(f'Column {PT_VAR_COL} is not in COSMIC DB.')
            pt_df = self.db_df[self.db_df[PT_VAR_COL] == pt_var_key]
            if len(pt_df) == 0:
                return 0
            else:
                # multiple different genetic mutations led to the same protein variant
                # return the sum of these mutations
                return sum(pt_df['Occurrences'])
