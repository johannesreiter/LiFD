#!/usr/bin/python
"""Script to run LiFD from the command line"""
import logging
import sys
import os
import argparse
import pandas as pd
from pathlib import Path

from lifd.lifd import LiFD
import lifd.utils as utils
from lifd.version import __version__ as lifd_version
import lifd.settings as settings
from lifd.databases.hotspots_database import HotspotsDB
from lifd.databases.cgi_database import CgiDB
from lifd.databases.oncokb_database import OncoKBDB
from lifd.databases.cosmic_db import CosmicDB
from lifd.predictors.vep import Vep
from lifd.predictors.candra import Candra
from lifd.predictors.cravat import Cravat
from lifd.predictors.cgi import Cgi
from lifd.predictors.fathmm import FatHMM

__author__ = 'Johannes Reiter'
__date__ = 'February 21, 2020'

# get logger
logger = logging.getLogger('lifd.{}'.format(__name__))


def usage():
    """
    Give the user feedback on how to call the tool
    Terminates the tool afterwards
    """
    logger.warning("Usage: lifd -i <INPUT_FILE> [--driver_genes=<driver_gene_file>] [--verbose] "
                   "[--hotspots] [--oncokb] [--cosmic] [--oncogenic] [--vep] [--cravat] [--fathmm] [--candra] "
                   "[--cgi] [--cgi_username=<YOUR CGI USERNAME>] [--cgi_token=<YOUR CGI TOKEN>] \n")
    logger.warning("Example: lifd -i lifd_examples/example_variants.xlsx --hotspots --oncokb --cosmic --oncogenic "
                   "--vep --cravat --cgi --fathmm --candra ")
    sys.exit(2)


def main():

    parser = argparse.ArgumentParser(description='Predict likely functional driver (LiFD) mutations.')

    parser.add_argument('-i', '--input', help='path to input file (types: CSV | TSV | XLSX)', type=str, required=True)
    parser.add_argument('--driver_genes', help='path to csv-file with driver genes', type=str,
                        default=settings.DR_LIST_FP)
    parser.add_argument('--verbose', action='store_true', help='Run LiFD in DEBUG logging level.')

    parser.add_argument('--hotspots', action='store_true', help='Annotate with Cancer Hotspots.')
    parser.add_argument('--hotspots_fp', help='path to cancer hotspots file', type=str,
                        default=settings.HOTSPOTS_FP)

    parser.add_argument('--oncokb', action='store_true', help='Annotate with OncoKB.')
    parser.add_argument('--oncokb_fp', help='path to oncokb file of all annotated variants', type=str,
                        default=settings.ONCOKB_ALLVARS_FP)

    parser.add_argument('--cosmic', action='store_true', help='Annotate with COSMIC genome screens data.')
    parser.add_argument('--cosmic_fp', help='path to COSMIC genome screens file', type=str,
                        default=settings.COSMIC_VARS_FP)

    parser.add_argument('--oncogenic', action='store_true', help='Annotate with CGI oncogenic variants.')
    parser.add_argument('--oncogenic_fp', help='path to CGI oncogenic variants file', type=str,
                        default=settings.ONCOGENIC_VARS_FP)

    parser.add_argument('--vep', action='store_true', help='Run VEP on data.')
    parser.add_argument('--cravat', action='store_true', help='Run Cravat/ChasmPLUS on data.')
    parser.add_argument('--fathmm', action='store_true', help='Run FatHMM on data.')
    parser.add_argument('--candra', action='store_true', help='Run Candra on data.')
    parser.add_argument('--cgi', action='store_true', help='Run CGI on data.')
    parser.add_argument('--cgi_username', help='CGI username', type=str)
    parser.add_argument('--cgi_token', help='CGI username', type=str)

    args = parser.parse_args()

    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)
        logger.info('Run LiFD in verbose mode.')

    input_fp = args.input
    input_file_ending = os.path.splitext(input_fp)[1].lower()
    if input_file_ending == '.csv':
        logger.info('Running LiFD {} for csv-file: {}'.format(lifd_version, input_fp))
        var_df = pd.read_csv(input_fp)
    elif input_file_ending == '.tsv':
        logger.info('Running LiFD {} for tsv-file: {}'.format(lifd_version, input_fp))
        var_df = pd.read_csv(input_fp, sep='\t')
    elif input_file_ending == '.xlsx' or input_file_ending == '.xls':
        logger.info('Running LiFD {} for xlsx-file: {}'.format(lifd_version, input_fp))
        var_df = pd.read_excel(input_fp)
    else:
        logger.error('Provided file {} with ending {} is not supported.'.format(input_fp, input_file_ending))
        logger.info('Exiting')
        sys.exit(2)

    # Check required columns
    # CancerType according to TCGA abbreviations:
    # https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/tcga-study-abbreviations
    logger.info('Found {} variants in provided file.'.format(len(var_df)))
    output_dir = utils.init()

    # ############ set up databases for first phase of LiFD ##############
    dbs = list()
    if args.hotspots or args.hotspots_fp != settings.HOTSPOTS_FP:
        hs_db = HotspotsDB(args.hotspots_fp)
        dbs.append(hs_db)

    if args.oncokb:
        oncokb_db = OncoKBDB(args.oncokb_fp)
        dbs.append(oncokb_db)

    if args.cosmic:
        cosmic_db = CosmicDB(args.cosmic_fp)
        dbs.append(cosmic_db)

    if args.oncogenic:
        oncogenic_db = CgiDB(args.oncogenic_fp)
        dbs.append(oncogenic_db)

    # ############# set up predictors for second phase of LiFD #################
    prds = list()
    if args.vep:
        prds.append(Vep)

    if args.cravat:   # ChasmPLUS
        prds.append(Cravat)

    if args.candra:
        prds.append(Candra)

    if args.fathmm:
        prds.append(FatHMM)

    if args.cgi:
        if args.cgi_username and args.cgi_token:
            Cgi.set_login(args.cgi_username, args.cgi_token)
            prds.append(Cgi)
        else:
            try:
                from lifd.cgi_settings import CGI_USER_ID, CGI_TOKEN

                Cgi.set_login(CGI_USER_ID, CGI_TOKEN)
                prds.append(Cgi)

            except (ModuleNotFoundError, ImportError) as e:
                logger.error('No CGI username and token provided and hence CGI can not be called: {}'.format(e))
                logger.error('CGI requires a username and a token which can be provided as arguments '
                             '--cgi_username=<YOUR_CGI_USERNAME> and --cgi_token=<YOUR_CGI_TOKEN> or '
                             'a cgi_settings.py file within the lifd directory with this content is provided: \n'
                             'CGI_USER_ID = <YOUR_CGI_USERNAME>\n'
                             'CGI_TOKEN = <YOUR_CGI_TOKEN>')

    # ############# run LiFD ####################
    lifd = LiFD(databases=dbs, predictors=prds,
                driver_gene_fp=args.driver_genes, driver_gene_col='TCGADrClf', )

    input_p = Path(input_fp)
    output_fn = '{}_LiFDed.xlsx'.format(input_p.stem)
    # output_fp = os.path.join(input_p.parent, output_fn)
    var_df = lifd.run_lifd(var_df, output_dir=output_dir, export_fn=output_fn)

    logger.info('LiFD {} finished annotating {} variants available.'.format(lifd_version, len(var_df)))
