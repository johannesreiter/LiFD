# default settings for running LiFD
import os

# default home directory that contains LiFD_dev
HOME_DIR = os.path.join('/home/groups/jgreiter', 'canary_crest_2019', 'LiFD') # os.path.join('~', 'LiFD')

# default database directory
DB_DIR = os.path.join(HOME_DIR, 'external_dbs')

# default model predictor directory
MDL_DIR = os.path.join(HOME_DIR, 'external_mdls')

# list of putative driver genes (default: TCGA consensus driver gene list)
DR_LIST_FP = os.path.join(DB_DIR, 'BaileyDing2018_driverconsensus.csv')

# see Tate et al., Nucleic Acids Res 2019
COSMIC_VARS_FP = os.path.join(DB_DIR, 'CosmicGenomeScreensMutantExport.tsv')

# catalog of validated oncogenic mutations Tamborero et al, Genome Medicine 2018
# https://www.cancergenomeinterpreter.org/mutations
# based on DoCM⁠ (PMID:27684579), ClinVar (PMID:26582918)⁠, OncoKB⁠ (PMID:28890946), and IARC (PMID:17311302) databases
ONCOGENIC_VARS_FP = os.path.join(DB_DIR, 'catalog_of_validated_oncogenic_mutations.tsv')

# Cancer Hotspots from Chang et al, Cancer Discovery 2018
# downloaded from: https://www.cancerhotspots.org/#/download
HOTSPOTS_FP = os.path.join(DB_DIR, 'ChangTaylor_hotspots_cd2018_v2.xls')

# Oncogenic mutation fromChakravarty et al., JCO PO 2017
# see https://github.com/oncokb/oncokb-public/tree/master/data for newest data releases
ONCOKB_ALLVARS_FP = os.path.join(DB_DIR, 'allAnnotatedVariants.txt')

# VEP cache file location
VEP_CACHE = os.path.join(HOME_DIR, '.vep')

# default output directory for output files
OUTPUT_DIR = 'TMP'