# default settings for running LiFD
import os

# default home directory that contains the static data
HOME_DIR = os.path.expanduser('~')

# default database directory
DB_DIR = os.path.join(HOME_DIR, 'databases')

# default predictor directory
PRD_DIR = os.path.join(HOME_DIR, 'predictors')

VEP_DIR = os.path.join(HOME_DIR, 'software', 'ensembl-vep')

# list of putative driver genes (default: TCGA consensus driver gene list)
DR_LIST_FP = os.path.join(DB_DIR, 'BaileyDing2018_driverconsensus.csv')

# see Tate et al., Nucleic Acids Res 2019
COSMIC_VARS_FP = os.path.join(DB_DIR, 'CosmicGenomeScreensMutantExport_grch37_v89.tsv')

# catalog of validated oncogenic mutations Tamborero et al, Genome Medicine 2018
# https://www.cancergenomeinterpreter.org/mutations
# based on DoCM⁠ (PMID:27684579), ClinVar (PMID:26582918)⁠, OncoKB⁠ (PMID:28890946), and IARC (PMID:17311302) databases
# # 'catalog_of_validated_oncogenic_mutations.tsv'
ONCOGENIC_VARS_FP = os.path.join(DB_DIR, 'TamboreroLopezBigas_cgi_oncogenic_gm2018.tsv')

# Cancer Hotspots from Chang et al, Cancer Discovery 2018
# downloaded from: https://www.cancerhotspots.org/#/download
HOTSPOTS_FP = os.path.join(DB_DIR, 'ChangTaylor_hotspots_cd2018_v2.xls')

# Oncogenic mutation from Chakravarty et al., JCO PO 2017
ONCOKB_INFO_URL = 'https://www.oncokb.org/api/v1/info'
ONCOKB_ALLVARS_URL = 'https://oncokb.org:443/api/v1/utils/allAnnotatedVariants'
# ONCOKB_ALLVARS_FP = None
ONCOKB_ALLVARS_FP = os.path.join(DB_DIR, 'oncoKB_allAnnotatedVariants_v2.1.tsv')

# DEPRECATED: no longer provided as a downloadable file
# see https://github.com/oncokb/oncokb-public/tree/master/data for newest data releases
# ONCOKB_ALLVARS_FP = os.path.join(DB_DIR, 'oncoKB_allAnnotatedVariants_v1.21.txt')

# VEP cache file location
VEP_CACHE = os.path.join(HOME_DIR, '.vep')

# path to fasta file of reference genome
REF_GENOME_FA_FP = os.path.join(DB_DIR, 'hg19.fa')

# default output directory for output files
OUTPUT_DIR = 'LiFD_TMP'