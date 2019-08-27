
from lifd.utils import FUNC_COL, NT_VAR_COL, FATHMM_KEY_COL

__author__ = 'Johannes REITER'
__date__ = 'Jan 19, 2019'


# mirroring mini variants dataframe with a few driver and passenger mutations
CANCER_TYPE = 'PAAD'
VARIANTS = [{'Subject': 'Pat1', 'Chromosome': '12', 'StartPosition': 25398284, 'EndPosition': 25398284,
            'ReferenceAllele': 'C', 'AlternateAllele': 'A', FUNC_COL: True, 'CancerType': CANCER_TYPE,
             NT_VAR_COL: '12__25398284__C__A', FATHMM_KEY_COL: 'ENSP00000256078__G12V', },  # KRAS
            {'Subject': 'Pat2', 'Chromosome': '18', 'StartPosition': 48591889, 'EndPosition': 48591889,
            'ReferenceAllele': 'A', 'AlternateAllele': 'G', FUNC_COL: True, 'CancerType': CANCER_TYPE,
             NT_VAR_COL: '18__48591889__A__G', FATHMM_KEY_COL: 'ENSP00000341551__D351G', },  # SMAD4
            {'Subject': 'Pat2', 'Chromosome': '12', 'StartPosition': 25398284, 'EndPosition': 25398284,
            'ReferenceAllele': 'C', 'AlternateAllele': 'A', FUNC_COL: True, 'CancerType': CANCER_TYPE,
             NT_VAR_COL: '12__25398284__C__A', FATHMM_KEY_COL: 'ENSP00000256078__G12V', },  # KRAS
            {'Subject': 'Pat@', 'Chromosome': '1', 'StartPosition': 107599898, 'EndPosition': 107599898,
            'ReferenceAllele': 'C', 'AlternateAllele': 'A', FUNC_COL: False, 'CancerType': CANCER_TYPE,
             NT_VAR_COL: '1__107599898__C__A', FATHMM_KEY_COL: 'ENSP00000359095__A128A'},  # PRMT6 silent
            ]

N_FUNC_VARS = sum([1 for v in VARIANTS if v[FUNC_COL]])

N_FUNC_UNIQUE_VARS = len({v[NT_VAR_COL] for v in VARIANTS if v[FUNC_COL]})

VARS_MIN_INPUT = [
    {'Subject': 'Pat1', 'Chromosome': '12', 'StartPosition': 25398284, 'EndPosition': 25398284,
     'ReferenceAllele': 'C', 'AlternateAllele': 'A', 'CancerType': CANCER_TYPE, },  # KRAS, missense
    {'Subject': 'Pat1', 'Chromosome': '19', 'StartPosition': 48994757, 'EndPosition': 48994756,
     'ReferenceAllele': '-', 'AlternateAllele': 'G', 'CancerType': CANCER_TYPE, },  # LMTK3, frameshift var, stop gained
    {'Subject': 'Pat1', 'Chromosome': '17', 'StartPosition': 7577568, 'EndPosition': 7577568,
     'ReferenceAllele': 'C', 'AlternateAllele': 'T', 'CancerType': CANCER_TYPE, },  # TP53
    {'Subject': 'Pat1', 'Chromosome': '1', 'StartPosition': 27058029, 'EndPosition': 27058029,
     'ReferenceAllele': 'T', 'AlternateAllele': 'G', 'CancerType': CANCER_TYPE, },  # ARID1A, stop_gained
    {'Subject': 'Pat2', 'Chromosome': '18', 'StartPosition': 48591889, 'EndPosition': 48591889,
     'ReferenceAllele': 'A', 'AlternateAllele': 'G', 'CancerType': CANCER_TYPE, },  # SMAD4
    {'Subject': 'Pat2', 'Chromosome': '12', 'StartPosition': 25398284, 'EndPosition': 25398284,
     'ReferenceAllele': 'C', 'AlternateAllele': 'A', 'CancerType': CANCER_TYPE, },  # KRAS
    {'Subject': 'Pat3', 'Chromosome': '12', 'StartPosition': 25380275, 'EndPosition': 25380275,
     'ReferenceAllele': 'T', 'AlternateAllele': 'A', 'CancerType': CANCER_TYPE, },  # KRAS
    {'Subject': 'Pat@', 'Chromosome': '12', 'StartPosition': 120594740, 'EndPosition': 120594740,
     'ReferenceAllele': 'G', 'AlternateAllele': 'A', 'CancerType': CANCER_TYPE},  # GCN1L1 silent, not func
    ]
