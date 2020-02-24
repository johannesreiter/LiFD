""" Prediction mutation effects by using VarCode """
import logging
import numpy as np
import pandas as pd
import sys
import traceback

from lifd.utils import NT_VAR_COL, PT_VAR_COL, FUNC_COL, MUT_EFFECT_COL, FATHMM_KEY_COL

# https://github.com/hammerlab/varcode
from varcode import Variant
from varcode.effects import predict_variant_effect_on_transcript
from varcode.effects.effect_classes import *

__author__ = 'Johannes REITER'
__date__ = 'March 2, 2017'


# get logger
logger = logging.getLogger('lifd.{}'.format(__name__))

GENE_SYMBOL_COL = 'GeneSymbol'
TRANSCRIPT_ID_COL = 'Transcript_ID'
PROTEIN_ID_COL = 'Protein_ID'

VAR_EFFECT_COLS = {NT_VAR_COL, PT_VAR_COL, GENE_SYMBOL_COL, TRANSCRIPT_ID_COL, PROTEIN_ID_COL,
                   FATHMM_KEY_COL, MUT_EFFECT_COL, FUNC_COL}


class MutationEffects:

    ensembl = None

    @staticmethod
    def set_ReferenceGenome(reference_genome='hg19'):
        """
        Set reference genome for mutation effect (e.g., missense) analysis
        :param reference_genome:
        :return:
        """
        # release 75 uses human reference genome GRCh37/hg19
        if reference_genome == 'hg19':
            from pyensembl import ensembl_grch37
            MutationEffects.ensembl = ensembl_grch37

        # release 77 uses human reference genome GRCh38/hg38
        elif reference_genome == 'hg20':
            from pyensembl import ensembl_grch38
            MutationEffects.ensembl = ensembl_grch38

        # release 54 uses human reference genome GRCh36/hg18
        elif reference_genome == 'hg18':
            from pyensembl import ensembl_grch36
            MutationEffects.ensembl = ensembl_grch36

        else:
            raise RuntimeError('Unknown reference genome: {}'.format(reference_genome)
                               + 'Only reference genomes hg18, hg19, and hg20 are supported.')

        logger.info('Set reference genome to {}.'.format(reference_genome))


def annotate_effects(df):
    """
    Run Varcode/pyensembl and add MutationEffect and Functional columns to the given dataframe
    Requires columns Chromosome, StartPosition, ReferenceAllele, and AlternateAllele
    :return: dataframe extended by mutation effect and whether mutation is potentially functional
    """

    # check if a reference genome has been set
    if MutationEffects.ensembl is None:
        MutationEffects.set_ReferenceGenome()

    req_cols = ['Chromosome', 'StartPosition', 'ReferenceAllele', 'AlternateAllele']
    for req_col in req_cols:
        if req_col not in df.columns.values:
            logger.error(
                '{} column is missing and variant annotation cannot be performed!'.format(
                    req_col))
            logger.error('Provided columns: '+', '.join(c for c in df.columns.values))
            raise RuntimeError('{} column is missing!'.format(req_col))

    if all(c in df.columns.values for c in VAR_EFFECT_COLS):
        logger.debug('All mutation effect columns present in provided dataframe.')
        return df

    nt_keys = list()
    pt_keys = list()
    gene_names = list()
    transcript_ids = list()
    mut_effects = list()
    potentially_functional = list()
    for index, row in df.iterrows():
        nt_var_key = f'{row.Chromosome}__{row.StartPosition}__{row.ReferenceAllele}__{row.AlternateAllele}'
        v = Variant(contig=row.Chromosome, start=row.StartPosition, ref=row.ReferenceAllele,
                    alt=row.AlternateAllele, ensembl=MutationEffects.ensembl)
        pt_var_key, gene_name, transcript_id, mut_effect = get_mutation_details(v)

        nt_keys.append(nt_var_key)
        pt_keys.append(pt_var_key)
        gene_names.append(gene_name)
        transcript_ids.append(transcript_id)
        mut_effects.append(mut_effect)

        if is_potentially_functional(mut_effect):
            potentially_functional.append(True)
        else:
            potentially_functional.append(False)

    if NT_VAR_COL not in df.columns:
        df[NT_VAR_COL] = pd.Series(nt_keys)

    if PT_VAR_COL not in df.columns:
        df[PT_VAR_COL] = pd.Series(pt_keys)

    if GENE_SYMBOL_COL not in df.columns:
        df[GENE_SYMBOL_COL] = pd.Series(gene_names)

    df[TRANSCRIPT_ID_COL] = pd.Series(transcript_ids)

    def _get_protein_id(ts_id):
        if ts_id is None or pd.isnull(ts_id):
            return np.nan
        else:
            return MutationEffects.ensembl.transcript_by_id(ts_id).protein_id

    df[PROTEIN_ID_COL] = df.apply(lambda row: _get_protein_id(row[TRANSCRIPT_ID_COL]), axis=1)

    def _get_fathmm_key(protein_id, pt_key):
        if pd.isnull(protein_id) or pd.isnull(pt_key):
            return np.nan
        else:
            return '{}__{}'.format(protein_id, pt_key.split('__')[1])

    df[FATHMM_KEY_COL] = df.apply(lambda row: _get_fathmm_key(row[PROTEIN_ID_COL], row[PT_VAR_COL]),
                                  axis=1)

    df[MUT_EFFECT_COL] = pd.Series(mut_effects)
    df[FUNC_COL] = pd.Series(potentially_functional)

    return df


def get_mutation_details(variant):
    """
    Return the name of the most likely effect of this mutation
    :param variant: varcode variant, see https://github.com/hammerlab/varcode
    :return: effect name
    """
    try:
        top_priority_effect = variant.effects().top_priority_effect()

        pt_var_key = f'{top_priority_effect.gene_name}__{top_priority_effect.short_description[2:]}'
        gene_name = top_priority_effect.gene_name
        transcript_id = top_priority_effect.transcript_id
        return pt_var_key, gene_name, transcript_id, type(top_priority_effect).__name__
        # return get_variant_effect_longest_transcript(variant)

    except BaseException as e:
        logger.debug('BaseException Error: {}'.format(str(e)))
        if logger.getEffectiveLevel() == logging.DEBUG:
            traceback.print_tb(sys.exc_info()[2])
        logger.warning('Mutation effect for variant {} could not be inferred.'.format(variant.short_description))
        return np.nan, np.nan, np.nan, np.nan


def get_variant_effect_longest_transcript(variant):
    """
    Predict variant effect based on longest transcript and return string
    :param variant: varcode variant, see https://github.com/hammerlab/varcode
    :return: effect name
    """

    try:
        ltc = _get_longest_transcript(variant.transcripts)
        if ltc is None:
            ve = Intergenic(variant).__class__.__name__
        else:
            ve = predict_variant_effect_on_transcript(variant, ltc).__class__.__name__

        return ve

    except:
        logger.warning('Mutation effect for variant {} could not be inferred.'.format(variant.short_description))
        return 'unknown'


def _get_longest_transcript(transcripts):
    """
    Returns the longest transcript from a list of given transcripts
    """
    longest_tc = None
    max_length = 0
    for tc in transcripts:
        if len(tc) > max_length:
            longest_tc = tc
            max_length = len(tc)

    return longest_tc


def _get_all_subclasses(cls):
    """
    Generator of all a class's subclasses
    """
    try:
        for subclass1 in cls.__subclasses__():
            yield subclass1
            for subclass2 in _get_all_subclasses(subclass1):
                yield subclass2
    except TypeError:
        return


def is_potentially_functional(effect_name):
    """
    Can this variant effect the expression of a gene?
    :param effect_name: varcode variant top priority effect, see https://github.com/hammerlab/varcode
    :return: True for non-synonymous mutations and variants that affect splicing
    """
    if effect_name is None or (isinstance(effect_name, float) and np.isnan(effect_name)):
        return np.nan

    # AlternateStartCodon: Replace annotated start codon with alternative start codon (e.g. "ATG>CAG")
    if effect_name == 'AlternateStartCodon':
        return True

    # ComplexSubstitution: Insertion and deletion of multiple amino acids
    elif effect_name == 'ComplexSubstitution':
        return True

    # Deletion:	Coding mutation which causes deletion of amino acid(s)
    elif effect_name == 'Deletion':
        return True

    # ExonLoss:	Deletion of entire exon, significantly disrupts protein
    elif effect_name == 'ExonLoss':
        return True

    # ExonicSpliceSite:	Mutation at the beginning or end of an exon, may affect splicing
    elif effect_name == 'ExonicSpliceSite':
        return True

    # FivePrimeUTR:	Variant affects 5' untranslated region before start codon
    elif effect_name == 'FivePrimeUTR':
        return False

    # FrameShiftTruncation:	A frameshift which leads immediately to a stop codon (no novel amino acids created)
    elif effect_name == 'FrameShiftTruncation':
        return True

    # FrameShift: Out-of-frame insertion or deletion of nucleotides, causes novel protein sequence
    # and often premature stop codon
    elif effect_name == 'FrameShift':
        return True

    # IncompleteTranscript:	Can't determine effect since transcript annotation is incomplete
    # (often missing either the start or stop codon)
    elif effect_name == 'IncompleteTranscript':
        return False

    # Insertion: Coding mutation which causes insertion of amino acid(s)
    elif effect_name == 'Insertion':
        return True

    # Intergenic: Occurs outside of any annotated gene
    elif effect_name == 'Intergenic':
        return False

    # Intragenic: Within the annotated boundaries of a gene but not in a region that's transcribed into pre-mRNA
    elif effect_name == 'Intragenic':
        return False

    # IntronicSpliceSite: Mutation near the beginning or end of an intron but less likely
    # to affect splicing than donor/acceptor mutations
    elif effect_name == 'IntronicSpliceSite':
        return True

    # Intronic:	Variant occurs between exons and is unlikely to affect splicing
    elif effect_name == 'Intronic':
        return False

    # NoncodingTranscript: Transcript doesn't code for a protein
    elif effect_name == 'NoncodingTranscript':
        return False

    # PrematureStop: Insertion of stop codon, truncates protein
    elif effect_name == 'PrematureStop':
        return True

    # Silent: Mutation in coding sequence which does not change the amino acid sequence of the translated protein
    elif effect_name == 'Silent':
        return False

    # SpliceAcceptor: Mutation in the last two nucleotides of an intron, likely to affect splicing
    elif effect_name == 'SpliceAcceptor':
        return True

    # SpliceDonor: Mutation in the first two nucleotides of an intron, likely to affect splicing
    elif effect_name == 'SpliceDonor':
        return True

    # StartLoss: Mutation causes loss of start codon, likely result is that
    # an alternate start codon will be used down-stream (possibly in a different frame)
    elif effect_name == 'StartLoss':
        return True

    # StopLoss:	Loss of stop codon, causes extension of protein by translation of nucleotides from 3' UTR
    elif effect_name == 'StopLoss':
        return True

    # Substitution:	Coding mutation which causes simple substitution of one amino acid for another
    elif effect_name == 'Substitution':
        return True

    # ThreePrimeUTR: Variant affects 3' untranslated region after stop codon of mRNA
    elif effect_name == 'ThreePrimeUTR':
        return False

    elif effect_name == 'unknown':
        return True

    else:
        logger.warning('Variant effect {} is not yet supported.'.format(effect_name))
        return True
