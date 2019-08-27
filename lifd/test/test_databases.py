
import unittest
import os
import shutil
import pandas as pd
from io import StringIO

from lifd.databases.database import Database
from lifd.databases.hotspots_database import HotspotsDB
from lifd.databases.cgi_database import CgiDB
from lifd.databases.oncokb_database import OncoKBDB
from lifd.databases.cosmic_db import CosmicDB

__author__ = 'Johannes REITER'
__date__ = 'Jan 19, 2019'


# mirroring mini hotspot file instead of mocking
CGI_MEM_STRING = (
    'gene\tgdna\tprotein\ttranscript\tinfo\tcontext\tcancer_acronym\tsource\treference\n'
    'KRAS\tchr12:g.25398284C>A\tp.G12V\tENST00000256078\t'
    'CSQN=Missense;codon_pos=25398283-25398284-25398285;ref_codon_seq=GGT;aliases=ENSP00000256078;source=Ensembl\t' 
    'somatic\tOV__NSCLC__PA__THCA__AML__COREAD\tClinVar__DoCM__Martelotto__OncoKB\t' 
    'PMID:19075190__PMID: 20570890__PMID:18316791__PMID:21398618\n'
    'APC\tchr5:g.112174368A>G\tp.N1026S\tENST00000257430\t'
    'CSQN=Missense;codon_pos=112174367-112174368-112174369;ref_codon_seq=AAT;aliases=ENSP00000257430;source=Ensembl\t' 
    'somatic\tCANCER\tOncoKB\tPMID:18166348'
)

# mirroring mini hotspot file instead of mocking
HOTSPOTS_MEM_TSV_v1 = StringIO(
    'Hugo Symbol\tCodon\tAlt Common Codon Usage *\tVariant Amino Acid\tQ-value\tTumor Count\t'
    'Tumor Type Count\tValidation Level [a]\tTumor Type Composition\n'
    + 'BRAF\tV600\tNA\tE:520|K:33|R:4|V:1\t0\t558\t10\tLevel-3\t'
      'skcm:252|thca:237|coadread:44|luad:9|gbm:6|mmyl:5|hgg:2|lgg:1|kirp:1|hnsc:1\n'
    + 'KRAS\tG12\tNA\tD:255|V:227|C:115|R:71|A:43|S:23|G:1|F:1\t0\t736\t21\tLevel-3\t'
      'paad:290|luad:188|coadread:151|ucec:36|stad:14|mmyl:9|ucs:7|cesc:7|ov:6|brca:6|blca:5|'
      'laml:3|esca:3|skcm:2|kirp:2|gbm:2|thca:1|prad:1|npc:1|mds:1|lihc:1\n'
    + 'PIK3CA\tH1047\tNA\tR:236|L:37|Y:6|Q:4\t0\t283\t22\tLevel-3\t'
      'brca:171|ucec:20|coadread:20|hnsc:19|stad:15|ucs:5|skcm:3|luad:3|lihc:3|cesc:3|prad:2|paad:2|'
      'mbl:2|lusc:2|lgg:2|kirc:2|hgg:2|gbm:2|blca:2|ov:1|esca:1|acc:1\n'
    + 'PTEN\tR130\tNA\tQ:40|G:28|*:16|L:5|P:4\t5.33E-59\t93\t13\tLevel-3\t'
      'ucec:56|gbm:7|brca:7|coadread:7|ucs:3|cesc:3|skcm:2|lusm:2|hnsc:2|stad:1|prad:1|lusc:1|kirc:1\n'
)


HOTSPOTS_MEM_STRING_CSV_v2 = (
    'Hugo_Symbol,Amino_Acid_Position,log10_pvalue,Mutation_Count,Reference_Amino_Acid,Total_Mutations_in_Gene,Median_Allele_Freq_Rank,Allele_Freq_Rank,Variant_Amino_Acid,Codon_Change,Genomic_Position,Detailed_Cancer_Types,Organ_Types,Tri-nucleotides,Mutability,mu_protein,Total_Samples,Analysis_Type,qvalue,tm,qvalue_pancan,Is_repeat,seq,length,align100,pad12entropy,pad24entropy,pad36entropy,TP,reason,n_MSK,n_Retro,judgement,inNBT,inOncokb,ref,qvaluect,ct,Samples,,\n'
    + 'BRAF,600,-2670.50439392493,897,V:897,1450,0.277777777777778,,E:833,gTg/gAg:831|Gtg/Atg:29|GTg/AAg:24|gtG/gtA:4|GTg/AGg:4|gTg/gGg:3|gTG/gAA:2,7:140453136_862|7:140453137_29|7:140453135_6,skcm:787:355|thpa:486:293|coad:712:54|coadread:683:47|luad:2057:33|mup:42:16|thap:33:14|thpd:58:9|gbm:688:8|mm:275:6|paad:932:6|macr:38:5|cup:135:5|panet:86:3|hggnos:87:3|hgnec:11:3|ecd:6:3|lggnos:544:3|chol:152:3|read:149:3|blca:852:2|nbl:326:2|lgsoc:17:2|acrm:23:2|apxa:2:2|lch:4:2|gist:122:1|ulms:52:1|gnc:1:1|hnmucm:11:1|apad:17:1|ape:4:1|mdlc:54:1|bcl:9:1|ihch:104:1|ipmn:7:1|esmm:3:1|prcc:335:1|hnsc:643:1,skin:974:357|thyroid:618:316|bowel:1782:113|lung:2761:33|unk:357:21|cnsbrain:2270:20|blood:890:11|pancreas:1059:10|biliarytract:358:4|bladder:958:2|ovaryfallopiantube:699:2|headandneck:988:2|softtissue:739:1|uterus:618:1|breast:2561:1|lymph:366:1|esophagusstomach:1407:1|kidney:1304:1,GTG|ACT|CCC|TCC|CCA|TCA|GCC|CCT|CCG|TCT|GTC|GTA|CTT|ACC,0.0267097457794181,0.0315472574215,24592,"skin,thyroid,pancan,bowel,lung,cnsbrain,pancreas,blood,biliarytract,ovaryfallopiantube,bladder",0,BRAF 600,0,FALSE,,NA,1,1.36223734667793,1.31982163922173,1.33965330506364,TRUE,,287,610,RETAIN,TRUE,TRUE,V,0,skin,thyroid:316|skin:298|bowel:113|lung:33|cnsbrain:20|unk:19|blood:11|pancreas:9|biliarytract:4|bladder:2|ovaryfallopiantube:2|breast:1|headandneck:1|kidney:1|lymph:1|softtissue:1|uterus:1,,\n'
    + 'KRAS,12,-6794.2683848371,2175,G:2175,2885,0.454545454545455,,V:657,gGt/gAt:757|gGt/gTt:657|Ggt/Tgt:380|Ggt/Cgt:166|gGt/gCt:134|Ggt/Agt:67|GGt/TTt:6|gcTGgt/gcCTgt:3|GGt/TCt:1|ggT/ggG:1|GGt/CTt:1|GGt/ATt:1|gcTGgt/gcATgt:1,12:25398284_1557|12:25398285_617|12:25398283_1,paad:932:728|luad:2057:535|coad:712:212|coadread:683:172|uec:339:60|maap:47:35|read:149:34|stad:748:28|cup:135:28|blca:852:27|ucs:117:19|macr:38:16|sbc:33:16|sem:59:11|unk:146:11|esca:556:10|mm:275:9|ehch:27:8|apad:17:8|lusc:346:8|nsgct:152:8|chol:152:7|paasc:8:7|cesc:205:7|brca:1345:6|mov:9:6|lune:37:6|ecad:15:6|nsclcpd:23:6|nsclc:16:6|soc:468:6|ihch:104:5|ampca:14:5|pampca:9:5|prcc:335:4|hgnec:11:4|prad:1366:4|lgsoc:17:4|idc:870:4|skcm:787:4|utuc:76:3|gbm:688:3|blad:7:3|ucp:4:3|ipmn:7:3|aml:198:3|dlbcl:246:2|uccc:14:2|ccov:24:2|ilc:184:2|thpa:486:2|brcanos:84:2|uad:4:2|hgsoc:132:2|luas:10:2|usc:46:2|mcn:2:2|spcc:3:2|sarc:280:2|sarcl:13:2|rms:50:2|npc:66:1|gcemu:5:1|umec:18:1|slct:2:1|mfh:53:1|osmca:1:1|sbmov:2:1|oovg:1:1|mup:42:1|lclc:1:1|ginet:8:1|lxsc:14:1|luca:9:1|sbov:2:1|paac:14:1|fibs:3:1|ansc:30:1|gej:26:1|lupc:1:1|vmm:12:1|scb:2:1|mdlc:54:1|icemu:1:1|hggnos:87:1|thpd:58:1|gist:122:1|ocsc:38:1|gccap:7:1|es:229:1|panet:86:1|pdc:8:1|urcc:48:1|adnos:4:1|aodg:43:1|ssrcc:12:1|imt:1:1|phch:1:1|tet:28:1|hcc:620:1|lggnos:544:1|thym:125:1|mds:28:1,pancreas:1059:745|lung:2761:570|bowel:1782:500|uterus:618:84|unk:357:42|esophagusstomach:1407:40|bladder:958:36|ovaryfallopiantube:699:25|biliarytract:358:21|testis:217:19|breast:2561:15|cervix:239:15|blood:890:13|ampullaofvater:23:10|softtissue:739:8|cnsbrain:2270:6|kidney:1304:5|prostate:1379:4|skin:974:4|headandneck:988:3|thyroid:618:3|lymph:366:2|thymus:162:2|vulvavagina:16:1|bone:297:1|liver:636:1,ACC|CCA|CCC|TCC|TTT|TCA|GTT|GTG|CTC|CCT|CCG,0.0316827126745203,0.0252344822022659,24592,"bowel,lung,pancreas,pancan,uterus,esophagusstomach,bladder,ovaryfallopiantube,biliarytract,cervix,breast,testis,ampullaofvater,blood,softtissue,cnsbrain,kidney,prostate,thyroid,headandneck,skin,thymus",0,KRAS 12,0,FALSE,,NA,1,1.26043970320943,1.2950060272429,1.29721360464951,TRUE,,1260,915,RETAIN,TRUE,TRUE,G,0,bowel,pancreas:269|bowel:160|lung:114|uterus:33|ovaryfallopiantube:14|unk:12|testis:11|bladder:10|biliarytract:7|esophagusstomach:6|cervix:5|breast:4|blood:2|prostate:2|softtissue:2|thyroid:2|ampullaofvater:1|headandneck:1|kidney:1|lymph:1,,\n'
    + 'PTEN,130,-213.140830041323,168,R:168,1101,0.205769230769231,7,*:38,cGa/cAa:72|Cga/Gga:44|Cga/Tga:38|cGa/cCa:6|cGa/cTa:5|ggACGA/ggCAGT:1|cgA/cgT:1|Cga/Aga:1,10:89692905_83|10:89692904_83|10:89692906_1|10:89692903_1,uec:339:75|gbm:688:12|prad:1366:8|coadread:683:7|brca:1345:7|ucs:117:5|idc:870:5|coad:712:4|luad:2057:3|skcm:787:3|cesc:205:3|sclc:224:2|ansc:30:2|ilc:184:2|lusc:346:2|umec:18:2|stad:748:2|hnsc:643:2|chol:152:1|mcc:41:1|cscc:92:1|snsc:4:1|uelms:4:1|ocs:12:1|chrcc:88:1|thpd:58:1|ceas:5:1|aastr:82:1|ccov:24:1|bpscc:1:1|gist:122:1|luas:10:1|gsarc:10:1|cup:135:1|plemeso:37:1|gbc:74:1|hnscup:7:1|ccrcc:771:1|esca:556:1|sarc:280:1,uterus:618:83|breast:2561:14|cnsbrain:2270:14|bowel:1782:13|lung:2761:8|prostate:1379:8|skin:974:5|headandneck:988:4|cervix:239:4|esophagusstomach:1407:3|biliarytract:358:2|ovaryfallopiantube:699:2|kidney:1304:2|softtissue:739:2|thyroid:618:1|penis:6:1|unk:357:1|pleura:76:1,TCG|ACG|TTC|TCC,0.0796107267175334,0.0262309068636857,24592,"pancan,uterus,breast,bowel,prostate,lung,cnsbrain,skin,headandneck,cervix",1.91658796865309e-210,PTEN 130,1.91658796865309e-210,FALSE,NA,NA,1,1.16795714944678,1.2950060272429,1.29741538035332,TRUE,,68,100,RETAIN,TRUE,FALSE,R,1.35644920841025e-113,uterus,uterus:10|cnsbrain:9|prostate:4|breast:3|bowel:2|headandneck:2|skin:2|biliarytract:1|esophagusstomach:1|lung:1|ovaryfallopiantube:1|pleura:1|thyroid:1,,\n'
)

ONCOKB_MEM_STRING = (
    'Isoform\tRefSeq\tEntrez Gene ID\tHugo Symbol\tAlteration\tProtein Change\tOncogenicity\tMutation Effect\t'
    'PMIDs for Mutation Effect\tAbstracts for Mutation Effect\n'
    'ENST00000256078\tNM_033360.2\t3845\tKRAS\tG12V\tG12V\tOncogenic\tGain-of-function\t'
    '"20516123, 20949621, 25359494, 20570890, 12957286, 20147967"\n'
    'ENST00000269305\tNM_000546.5\t7157\tTP53\tG266R\tG266R\tLikely Oncogenic\tLoss-of-function\t16827139\n'
    'ENST00000111111\tNM_000111.1\t0000\tABCD\tG266R\tG266R\tLikely Neutral\tLoss-of-function\t12345678\n'
    'ENST00000111110\tNM_000111.0\t0001\tABCE\tG266R\tG266R\tInconclusive\tLoss-of-function\t12345678\n'
)

COSMIC_MEM_STRING = (
    'Gene name\tAccession Number\tGene CDS length\tHGNC ID\tSample name\tID_sample\tID_tumour\tPrimary site\t'
    'Site subtype 1\tSite subtype 2\tSite subtype 3\tPrimary histology\tHistology subtype 1\tHistology subtype 2\t' 
    'Histology subtype 3\tGenome-wide screen\tMutation ID\tMutation CDS\tMutation AA\tMutation Description\t'
    'Mutation zygosity\tLOH\tGRCh\tMutation genome position\tMutation strand\tSNP\tFATHMM prediction\tFATHMM score\t' 
    'Mutation somatic status\tPubmed_PMID\tID_STUDY\tSample Type\tTumour origin\tAge\n'
    + 'KRAS\tENST00000311936\t567\t6407\tIPMN26\t1691643\t1599910\tpancreas\tNS\tNS\tNS\tother\tadenoma\t'
      'NS\tNS\ty\tCOSM518\tc.34G>C\tp.G12R\tSubstitution - Missense\t\tu\t37\t12:25398285-25398285\t-\tn\tPATHOGENIC\t'
      '.98468\tConfirmed somatic variant\t22158988\tNS\tNS\t\n'
    # mutation on negative strand
    + 'KRAS\tENST00000311936\t567\t6407\tIPMN26\t16916431\t1599910\tpancreas\tNS\tNS\tNS\tother\tadenoma\t'
      'NS\tNS\ty\tCOSM518\tc.34G>C\tp.G12R\tSubstitution - Missense\t\tu\t37\t12:25398285-25398285\t-\tn\tPATHOGENIC\t'
      '.98468\tConfirmed somatic variant\t22158988\tNS\tNS\t\n'
    + 'PTEN\tENST00000371953\t1212\t9588\tPD3989a\t1280816\t1192107\tbreast\tNS\tNS\tNS\tcarcinoma\tductal_carcinoma\t' 
      'NS\tNS\ty\tCOSM4898\tc.950_953delTACT\tp.T319fs*1\tDeletion - Frameshift\t\tu\t37\t10:89720799-89720802\t+\t' 
      '\t\t\tConfirmed somatic variant\t27135926\tNS\tprimary\n'
    + 'MBD5\tENST00000407073\t4485\t20444\tTCGA-A3-3373-01\t1779792\t1683791\tkidney\tNS\tNS\tNS\tcarcinoma\t'
      'clear_cell_renal_cell_carcinoma\t' 
      'NS\tNS\ty\tCOSM475989\tc.2340T>G\tp.L780L\tSubstitution - coding silent\t\tu\t37\t2:149227852-149227852\t+\t' 
      'n\tPATHOGENIC\t.72461\tConfirmed somatic variant\t416\tfresh/frozen - NOS\tprimary\t54\n'
    # long deletion
    + 'TP53_ENST00000413465\tENST00000413465\t858\t11998\tSJHGG003_A\t2307333\t2172537\tcentral_nervous_system\t'
      'brainstem\tNS\tNS\tglioma\tastrocytoma_Grade_IV\tglioblastoma_multiforme\tNS\ty\tCOSM4968949\tc.139_161del23\t'
      'p.P47fs*2\tDeletion - Frameshift\t\tu\t37\t17:7579526-7579548\t-\t\t\t\tConfirmed somatic variant\t24705251\t\t'
      'surgery - NOS\tNS\t15\n'
    # insertion
    + 'CHD3\tENST00000330494\t6003\t1918\tBD6T\t2459924\t2322761\tbiliary_tract\tbile_duct\tNS\tNS\tcarcinoma\t' 
      'NS\tNS\tNS\ty\tCOSM5499281\tc.1612_1613insCC\tp.R540fs*17\tInsertion - Frameshift\t\tu\t37\t' 
      '17:7798765-7798766\t+\t\t\t\tConfirmed somatic variant\t\t658\tNS\tprimary\t35\n'
    # insertion on negative strand to be complemented
    + 'MEMO1\tENST00000295065\t894\t14014\tT1222\t2658242\t2518401\tlarge_intestine\tNS\tNS\tNS\tcarcinoma\t' 
      'adenocarcinoma\tNS\tNS\ty\tCOSM6711426\tc.289_290insTA\tp.R97fs*66\tInsertion - Frameshift\t\tu\t37\t' 
      '2:32145902-32145903\t-\t\t\t\tConfirmed somatic variant\t27149842\t\tNS\tNS\t71.91\n'
    + 'KCNU1_ENST00000399881\tENST00000399881\t3450\t\tI2L-P19Ta-Tumor-Organoid\t2433499\t2296380\tlarge_intestine\t' 
      'colon\tascending\tNS\tcarcinoma\tadenocarcinoma\tNS\tNS\ty\tCOSM1456808\tc.2402delC\tp.P803fs*11\t' 
      'Deletion - Frameshift\thet\tu\t37\t8:36768518-36768518\t+\tn\t\t\tConfirmed somatic variant\t' 
      '25957691\t\tNS\tNS\t82\n'
    + 'ALDH1A3\tENST00000329841\t1539\t409\t59\t2370251\t2233049\tkidney\tNS\tNS\tNS\tcarcinoma\t' 
      'papillary_renal_cell_carcinoma\tNS\tNS\ty\tCOSM5016501\tc.1253_1254insA\tp.I419fs*16\t' 
      'Insertion - Frameshift\t\tu\t37\t15:101447345-101447346\t+\t\t\t\tConfirmed somatic variant\t' 
      '25401301\t\tsurgery fresh/frozen\tprimary\t60\n'
    # longer inseration on negative strand to be complemented
    + 'CEBPA\tENST00000498907\t1077\t1833\tPTC-6C\t2186174\t2054471\tthyroid\tNS\tNS\tNS\tother\tn' 
      'eoplasm\tNS\tNS\ty\tCOSM5446613\tc.568_569insCGCACC\tp.P196_P197insHP\tInsertion - In frame\t\tu\t37\t' 
      '19:33792752-33792753\t-\t\t\t\tConfirmed somatic variant\t\t589\tNS\tNS\t\n'
)

OUTPUT_DIR = os.path.join('..', 'output_testing')


class DatabaseTest(unittest.TestCase):

    def setUp(self):
        if not os.path.exists(OUTPUT_DIR):
            os.makedirs(OUTPUT_DIR)

        self.oncokb_io_tsv = StringIO(ONCOKB_MEM_STRING)
        self.cosmic_io_tsv = StringIO(COSMIC_MEM_STRING)
        self.cgi_io_tsv = StringIO(CGI_MEM_STRING)
        self.hotspots_io_csv = StringIO(HOTSPOTS_MEM_STRING_CSV_v2)

    def tearDown(self):
        if os.path.exists(OUTPUT_DIR) and os.path.isdir(OUTPUT_DIR):
            shutil.rmtree(OUTPUT_DIR)

    def test_Database(self):

        with self.assertRaises(TypeError):
            Database('', '')

    def test_CGIDatabase(self):
        # mirroring mini hotspot file instead of mocking
        cgi_db = CgiDB(self.cgi_io_tsv, lazy_loading=False)

        # does source filename get assigned correctly
        self.assertEqual(cgi_db.db_source, self.cgi_io_tsv)

        # was the TSV file read correctly
        self.assertIsNotNone(cgi_db.db_df)
        self.assertEqual(len(cgi_db.db_df), 2)

        nt_var_key = '12__25398284__C__A'
        self.assertTrue(cgi_db.in_database(nt_var_key=nt_var_key))

        nt_var_key = '5__112174368__A__G'
        self.assertTrue(cgi_db.in_database(nt_var_key=nt_var_key))

        nt_var_key = '6_12574934__C__T'
        self.assertFalse(cgi_db.in_database(nt_var_key=nt_var_key))
        
        with self.assertRaises(AttributeError):
            cgi_db.in_database(None, None)

    def test_CGIDatabase_lazy_loading(self):

        # mirroring mini hotspot file instead of mocking
        cgi_db = CgiDB(self.cgi_io_tsv, lazy_loading=True)

        # does source filename get assigned correctly
        self.assertEqual(cgi_db.db_source, self.cgi_io_tsv)

        # because of lazy loading DB should not yet been loaded
        self.assertIsNone(cgi_db.db_df)
        nt_var_key = '12__25398284__C__A'
        self.assertTrue(cgi_db.in_database(nt_var_key=nt_var_key))
        # after check of first variant, DB should be loaded
        self.assertIsNotNone(cgi_db.db_df)

    def test_CosmicDatabase(self):
        # mirroring mini hotspot file instead of mocking
        cosmic_db = CosmicDB(self.cosmic_io_tsv, lazy_loading=False)

        # does source filename get assigned correctly
        self.assertEqual(cosmic_db.db_source, self.cosmic_io_tsv)

        # was the TSV file read correctly and did the neutral variants get removed?
        self.assertIsNotNone(cosmic_db.db_df)
        self.assertEqual((9, 2), cosmic_db.db_df.shape)

        # test formatting of key when mutations are on negative strand
        nt_var_key = '12__25398285__C__G'
        self.assertEqual(2, cosmic_db.in_database(nt_var_key, None))

        nt_var_key = '12__25398285__G__@'
        self.assertEqual(0, cosmic_db.in_database(nt_var_key, None))

        # test formatting of longer deletions
        nt_var_key = '17__7579526__del23__-'
        self.assertEqual(1, cosmic_db.in_database(nt_var_key, None))

        # test formatting of small insertion
        nt_var_key = '17__7798765__-__CC'
        self.assertEqual(1, cosmic_db.in_database(nt_var_key, None))

        # test formatting of small insertion on negative strand that need to be complemented
        nt_var_key = '2__32145902__-__TA'
        self.assertEqual(1, cosmic_db.in_database(nt_var_key, None))

        # test deletion
        nt_var_key = '8__36768518__C__-'
        self.assertEqual(1, cosmic_db.in_database(nt_var_key, None))

        # test insertion
        nt_var_key = '15__101447345__-__A'
        self.assertEqual(1, cosmic_db.in_database(nt_var_key, None))

        # test complementing of longer insertion on negative strand
        nt_var_key = '19__33792752__-__GGTGCG'
        self.assertEqual(1, cosmic_db.in_database(nt_var_key, None))

        # no longer supported
        # pt_var_key = 'KRAS__G12R'
        # self.assertEqual(cosmic_db.in_database(None, pt_var_key), 2)
        # pt_var_key = 'MBD5__L780L'
        # self.assertEqual(cosmic_db.in_database(None, pt_var_key), 1)
        # pt_var_key = 'KRAS__G12V'
        # self.assertEqual(cosmic_db.in_database(None, pt_var_key), 0)

        with self.assertRaises(RuntimeError):
            cosmic_db.in_database(None, 'KRAS__G12R')

        with self.assertRaises(AttributeError):
            cosmic_db.in_database(None, None)

    # @unittest.skip('Cosmic is not yet implemented')
    def test_CosmicDatabase_lazy_loading(self):

        # mirroring mini hotspot file instead of mocking
        cosmic_db = CosmicDB(self.cosmic_io_tsv, lazy_loading=True)

        # does source filename get assigned correctly
        self.assertEqual(cosmic_db.db_source, self.cosmic_io_tsv)

        # because of lazy loading DB should not yet been loaded
        self.assertIsNone(cosmic_db.db_df)
        nt_var_key = '12__25398285__C__G'
        self.assertEqual(2, cosmic_db.in_database(nt_var_key, None))
        # after check of first variant, DB should be loaded
        self.assertIsNotNone(cosmic_db.db_df)

    def test_OncoKBDatabase(self):

        # mirroring mini hotspot file instead of mocking
        oncokb_db = OncoKBDB(db_source=self.oncokb_io_tsv, lazy_loading=False)

        # does source filename get assigned correctly
        self.assertEqual(oncokb_db.db_source, self.oncokb_io_tsv)

        # was the TSV file read correctly and did the neutral variants get removed?
        self.assertIsNotNone(oncokb_db.db_df)
        self.assertEqual(oncokb_db.db_df.shape, (2, 10))

        pt_var_key = 'KRAS__G12V'
        self.assertTrue(oncokb_db.in_database(None, pt_var_key))

        pt_var_key = 'TP53__G266R'
        self.assertTrue(oncokb_db.in_database(None, pt_var_key))

        pt_var_key = 'KRAI__G12V'
        self.assertFalse(oncokb_db.in_database(None, pt_var_key))

        pt_var_key = 'AVCD__G266R'
        self.assertFalse(oncokb_db.in_database(None, pt_var_key))

        with self.assertRaises(AttributeError):
            oncokb_db.in_database(None, None)

    def test_OncoKBDatabase_lazy_loading(self):

        # mirroring mini hotspot file instead of mocking
        oncokb_db = OncoKBDB(db_source=self.oncokb_io_tsv, lazy_loading=True)

        # does source filename get assigned correctly
        self.assertEqual(oncokb_db.db_source, self.oncokb_io_tsv)

        # because of lazy loading DB should not yet been loaded
        self.assertIsNone(oncokb_db.db_df)
        pt_var_key = 'KRAS__G12V'
        self.assertTrue(oncokb_db.in_database(None, pt_var_key))
        # after check of first variant, DB should be loaded
        self.assertIsNotNone(oncokb_db.db_df)

    def test_HotspotsDatabase(self):

        test_hs_df = pd.read_csv(self.hotspots_io_csv)

        test_hs_fp = os.path.join(OUTPUT_DIR, 'hotspots_test.xlsx')
        test_hs_df.to_excel(test_hs_fp, index=False)

        hdb = HotspotsDB(test_hs_fp, lazy_loading=False)

        # does source filename get assigned correctly
        self.assertEqual(hdb.db_source, test_hs_fp)

        # was the TSV file read correctly
        self.assertIsNotNone(hdb.db_df)
        self.assertEqual(hdb.db_df.shape, (3, 43))

        pt_var_key = 'KRAS__G12V'
        self.assertTrue(hdb.in_database(None, pt_var_key))

        pt_var_key = 'PTEN__R130*'
        self.assertTrue(hdb.in_database(None, pt_var_key))

        pt_var_key = 'KRAI__G12V'
        self.assertFalse(hdb.in_database(None, pt_var_key))

        # with self.assertRaises(AttributeError):
        #     hdb.in_database(None, None)

        self.assertTrue(pd.isnull(hdb.in_database(None, None)))

    def test_HotspotsDatabase_lazy_loading(self):

        test_hs_df = pd.read_csv(self.hotspots_io_csv)

        test_hs_fp = os.path.join(OUTPUT_DIR, 'hotspots_test.xlsx')
        test_hs_df.to_excel(test_hs_fp, index=False)

        hdb = HotspotsDB(test_hs_fp, lazy_loading=True)

        # does source filename get assigned correctly
        self.assertEqual(hdb.db_source, test_hs_fp)

        # because of lazy loading DB should not yet been loaded
        self.assertIsNone(hdb.db_df)
        pt_var_key = 'KRAS__G12V'
        self.assertTrue(hdb.in_database(None, pt_var_key))
        # after check of first variant, DB should be loaded
        self.assertIsNotNone(hdb.db_df)
