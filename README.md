# LiFD
## Pipeline to identify LiFDs (Likely Functional Drivers)
See corresponding papers for full details Reiter et al, *[An analysis of genetic heterogeneity in untreated cancers](https://doi.org/10.1038/s41568-019-0185-x)* (Nature Reviews Cancer, 2019) and Reiter, Makohon-Moore, Gerold et al, *[Minimal functional driver gene heterogeneity among untreated metastases](https://doi.org/10.1126/science.aat7171)* (Science, 2018) 

========

### What is LiFD?
LiFD is a two-phase algorithm to predict likely functional driver (LiFD) mutations that integrates information from multiple databases and bioinformatic methods. 
During the first phase of LiFD, variants that are present in OncoKB, the catalog of validated oncogenic mutations (CGI, Cancer Genome Interpreter; driver prediction “known”), known cancer hotspots, or present at least 4 times in COSMIC (Catalogue of Somatic Mutations in Cancer) are annotated as functional. 
If a variant is not annotated as functional in the first phase, LiFD uses CHASMplus, FATHMM, CanDrA+, CGI, and VEP to predict the functional consequences of individual mutations in the second phase.  
By default, LiFD requires a q-value of at most 0.1 for CHASMplus. 
For FATHMM, LiFD uses the recommended threshold of at most -0.75. 
For CanDrA, LiFD uses a significance threshold of 0.05 for driver predictions. 
For CGI, LiFD requires “predicted” in its driver classification column and for VEP LiFD requires “HIGH” in its predicted impact column. 
If the majority (>50%) of the methods that produce a valid result predict functionality, LiFD annotates the mutation as likely functional and otherwise as unlikely functional.
Default reference genome is hg19. Only some of the pooled tools support hg20.

### <a name="releases"> Releases
* LiFD 0.1.0 2019-08-27: Initial release.
* LiFD 0.1.1 2020-02-24: Integrated OncoKB upgrade to curl API. Added option to run LiFD from the command line.

### <a name="installation"> Installation and Setup
1. Install Python 3.6 ([https://www.python.org/downloads](https://www.python.org/downloads)). Check installation with ```which python3.6```. Load Python 3.6 with ```ml python/3.6``` if installing LiFD onto a remote server/cluster.
2. Install required packages with ```pip install numpy scipy pandas statsmodels xlrd openpyxl xlsxwriter```; NumPy ([http://www.numpy.org](http://www.numpy.org)), SciPy ([http://www.numpy.org](http://www.numpy.org)), pandas ([http://pandas.pydata.org/](http://pandas.pydata.org/)).
3. Install PyEnsembl (https://github.com/openvax/pyensembl) and Varcode (https://github.com/openvax/varcode), used for mutation effect annotation, with ```pip install pyensembl varcode```. Run ```pyensembl install --release 75 --species human``` to get the latest genome database for hg19/GRCh37.
4. Open a terminal and clone the repository from GitHub with ```git clone https://github.com/johannesreiter/LiFD.git```
and install LiFD to your python environment by running ```pip3.6 install -e <LiFD_directory>```
5. Test installation by opening a python shell ```python3.6``` and execute these commands ```import lifd``` and ```lifd.__version__```. As output you should get your current LiFD version. Exit the shell with ```exit()```.

### <a name="run_simple"> Running LiFD when functional annotation is provided
LiFD takes as input either a CSV or an Excel file where the following columns are required for the different methods. 
See directory ```examples``` and files ```SupplementaryTable1.xlsx``` or ```SupplementaryTable3.xlsx```.
If boolean column ```MaybeFunctional``` (denoting whether variant is a nonsynonymous or splice-site mutation) is not provided, then VarCode will be invoked to assess mutation effects.
Either a boolean column ```DriGeneClf``` (denoting whether variant occurred in putative driver gene) needs to be provided or a CSV file with list of putative driver genes needs to be given (e.g., ```examples/BaileyDing2018_driverconsensus.csv```).
Default settings can be configured in ```lifd/settings.py```
  - CGI (known): Requires boolean column ```In_CGI_Catalog```
  - COSMIC: Requires numerical column ```In_COSMIC```
  - Hotspots: Requires boolean column ```In_Hotspots```
  - OncoKB: Requires boolean column ```In_OncoKB```
  - CGI (prediction): Requires column ```CGI_driver```
  - CanDrA: Requires ```CanDrA_clf``` and ```CanDrA_sig```
  - CHASMplus: Requires ```CHASMplus_CT_Pvalue_corr```
  - FATHMM: Requires column ```FATHMM_score```
  - VEP: Requires column ```VEP_impact```
  
  
##### Usage: 
```shell
$ lifd -i <INPUT_FILE> [--driver_genes=<driver_gene_file>] [--verbose] [--hotspots] [--oncokb] [--cosmic] [--oncogenic] [--vep] [--cravat] [--fathmm] [--candra] [--cgi] [--cgi_username=<YOUR CGI USERNAME>] [--cgi_token=<YOUR CGI TOKEN>]
```

##### Parameters:
- ```-i <input_file>```: Path to input file (types: CSV | TSV | XLSX)
- ```--driver_genes=<driver_gene_file>``` Path to csv-file with driver genes
- ```--verbose```: Run LiFD in DEBUG logging level.
- ```--hotspots```: Annotate with Cancer Hotspots from Chang et al, Cancer Discovery 2018.
- ```--oncokb```: Annotate with oncogenic mutations (OncoKB) from Chakravarty et al., JCO PO 2017.
- ```--cosmic```: Annotate with COSMIC (Tate et al., Nucleic Acids Res 2019).
- ```--oncogenic```: Annotate with catalog of validated oncogenic mutations from Tamborero et al, Genome Medicine 2018.
- ```--vep```: Annotate with VEP.
- ```--cravat```: Annotate with Cravat/ChasmPLUS.
- ```--fathmm```: Annotate with FatHMM.
- ```--candra```: Annotate with Candra.
- ```--cgi```: Annotate with CGI which also requires the following two arguments with a username and a token
- ```--cgi_username=<YOUR CGI USERNAME>```: your CGI username to connect to the CGI interface
- ```--cgi_token=<YOUR CGI TOKEN>```: your CGI token to connect to the CGI interface


##### Example:
```shell
$ lifd -i lifd_examples/example_variants.xlsx --hotspots --oncokb --cosmic --oncogenic --vep --cravat --cgi --fathmm --candra
```
For an ipython notebook example see ```lifd_examples/intraprimary_heterogeneity_analysis.ipynb``` or ```lifd_examples/intermetastatic_heterogeneity_analysis.ipynb```


### <a name="run_complex"> Running LiFD when functional annotation is missing
To automatically invoke the various methods, various dependencies need to be configured which can be very cumbersome.

LiFD takes as input either a CSV, TSV, or an Excel file with the following required columns: ```Chromosome```, ```StartPosition```, ```EndPosition```, ```ReferenceAllele```, and ```AlternateAllele```.
See ```lifd_examples/example_variants.xlsx``` for format.
Some predictors also require a ```CancerType``` column according to the TCGA abbreviations ([https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/tcga-study-abbreviations](https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/tcga-study-abbreviations)). 
If boolean column ```MaybeFunctional``` (denoting whether variant is a nonsynonymous or splice-site mutation) is not provided, then VarCode will be invoked.
For an ipython notebook example see ```lifd_examples/lifd_examples.ipynb```.
To run LiFD from the command line see above.

Download the following databases into a directory and set variable ```DB_DIR``` with path to the databases in ```lifd/settings.py```:
  - OncoKB v1.21: apply to get access to their database at ([https://www.oncokb.org/](https://www.oncokb.org/)); previously the data was freely available ([https://github.com/oncokb/oncokb-public/tree/master/data](https://github.com/oncokb/oncokb-public/tree/master/data)). Set variable ```ONCOKB_ALLVARS_FP``` with path to downloaded file in ```lifd/settings.py``` accordingly.
  - Cancer Hotspots V2 ([https://www.cancerhotspots.org/#/download](https://www.cancerhotspots.org/#/download)). Set variable ```HOTSPOTS_FP``` with path to downloaded file in ```lifd/settings.py``` accordingly. Note that in the original study, version 1 of the hotspots were used.
  - COSMIC Mutation Data Genome Screens v89 ([https://cancer.sanger.ac.uk/cosmic/download?genome=37](https://cancer.sanger.ac.uk/cosmic/download?genome=37)). Set variable ```COSMIC_VARS_FP``` with path to downloaded file in ```lifd/settings.py``` accordingly. Note that in COSMIC release v90 many things were changed and the strand information is now always positive inconsistent with previous releases. Release v91 will fix this issue again.
  - CGI Catalog of Validated Oncogenic Mutations ([https://www.cancergenomeinterpreter.org/mutations](https://www.cancergenomeinterpreter.org/mutations)). Set variable ```ONCOGENIC_VARS_FP``` with path to downloaded file in ```lifd/settings.py``` accordingly.

Install Pysam (https://github.com/pysam-developers/pysam), used to find reference alleles for indels, with ```pip3.6 install pysam```. Download a fasta file for hg19 (http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz) into the database directory.

Setup the predictive tools (and their cancer-specific models) to create annotations for LiFD automatically:
  - CanDrA.v+ ([https://bioinformatics.mdanderson.org/main/CanDrA](https://bioinformatics.mdanderson.org/main/CanDrA))
  - CRAVAT and CHASMplus ([https://chasmplus.readthedocs.io/en/latest/](https://chasmplus.readthedocs.io/en/latest/))
  - CGI is a web tool and no installation is required. However, a login needs to be created at ([https://www.cancergenomeinterpreter.org](https://www.cancergenomeinterpreter.org)). A file ```cgi_settings.py``` containing the credentials is expected at ```lifd/cgi_settings.py``` in the following format:
    ```
    CGI_USER_ID = '<YOUR USERNAME>'
    CGI_TOKEN = '<YOUR CGI TOKEN>'
    ```
  - FATHMM ([http://fathmm.biocompute.org.uk/downloads.html](http://fathmm.biocompute.org.uk/downloads.html))
    - The files for FATHMM (```fathmm.py``` and ```config.ini```) should be located under ```DB_DIR``` in the folder ```fathmm/cgi-bin```. The following edits should be made to ```fathmm.py``` to make the file compatible with Python 3, if desired:
      - Change all lines "except Exception, e:" to "except Exception as e:"
      - Change line "import ConfigParser" to "import configparser as ConfigParser"
    - FATHMM requires MySQL (https://www.mysql.com/downloads/) or MariaDB. If using MariaDB, which is common for remote servers/clusters, make the additional change below to replace the line that initializes dbCursor:
        ```
        mariadb_connection = mariadb.connect(
          host = str(Config.get("DATABASE", "HOST")), 
          port = int(Config.get("DATABASE", "PORT")), 
          unix_socket = int(Config.get("DATABASE", "UNIX_SOCKET")), 
          user = str(Config.get("DATABASE", "USER")), 
          password = str(Config.get("DATABASE", "PASSWD")), 
          database = str(Config.get("DATABASE", "DB")))
        dbCursor = mariadb_connection.cursor(dictionary = True)
        ```
        Load MariaDB with ```ml system mariadb```.
        If using MySQL, start server with ```mysql.server start```
    - Check the line ```COMMAND = 'python3 fathmm.py -w Cancer {} {}'``` in ```lifd/predictors/fathmm.py``` and correct the line to reflect the version of Python. Note that a separate installation of Python 2 will be necessary to run fathmm.py if the above corrections are not made.
    - As written on the FATHMM installation page, place information about the user, database, etc. in ```config.ini``` file as follows:
    ```
    [DATABASE]
    HOST   = <MySQL Host>
    PORT   = <MySQL Port>
    USER   = <MySQL Username>
    PASSWD = <MySQL Password>
    DB     = fathmm
    ```
      It may be necessary to add ```UNIX_SOCKET = <MySQL/MariaDB Socket>``` to the FATHMM config file, especially for a remote server implementation. On a remote server, run ```mysqld_safe``` on the command line to start a mysqld server, which blocks for MySQL commands. Terminal multiplexers such as tmux are espsecially useful for this task.
    - Download MySQLdb with ```pip install MySQL-python``` for Python 2 and ```pip install mysqlclient``` for Python 3 if using MySQL. For MariaDB, download mysql.connector with ```pip install mysql-connector-python```. Delete all ```import MySQLdb statements``` and insert ```import mysql.connector as mariadb``` in the downloaded ```fathmm.py``` file.
    
  - VEP ([https://uswest.ensembl.org/info/docs/tools/vep/script/vep_download.html](https://uswest.ensembl.org/info/docs/tools/vep/script/vep_download.html))
    - VEP requires Perl (https://www.perl.org/get.html) and MySQL (https://www.mysql.com/downloads/). In a remote server, load Perl and MariaDB with ```ml perl``` and ```ml system mariadb```, respectively. Install the required Perl packages (Archive::Zip, DBD::mysql, DBI) with ```cpanm [package]``` before running the VEP setup file with ```perl INSTALL.pl```. The directory of the cache can be assigned with an optional argument (```--CACHEDIR [dir]```) during setup. Set variable ```VEP_CACHE``` with path to VEP cache (initially set to the default location of the cache).
    - Note that some dependencies of various Perl packages may fail during either installation of the requirements or the setup, such as ```Test::Pod::Coverage``` or ```XML::DOM::XPath```; keep track of these packages and install them separately, possibly using ```--force``` argument to circumnavigate outdated tests. Dependencies of dependencies may also fail, in which case this process should be repeated.
    - Note that many packages may already exist on a remote server. Load some of the packages with ```ml libgd``` and ```ml biology htslib```.



========

### Problems?
If you have any questions, you can contact us ([https://github.com/johannesreiter](https://github.com/johannesreiter)) and we will try to help.


### License
Copyright (C) 2019 Johannes Reiter

LiFD is licensed under the GNU General Public License, Version 3.
This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, 
version 3 of the License.
There is no warranty for this free software.

========

Author: Johannes G. Reiter, Stanford University, [https://reiterlab.stanford.edu]('https://reiterlab.stanford.edu')