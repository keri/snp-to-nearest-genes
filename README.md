# SNPs to Nearest Genes
Returns nearest gene names and numbers associated with one or more SNP inputs

# Purpose
After not seeing a simple solution to what I was looking for, I wrote this a very simple python script using the entrez library. It makes calls to the dbSNP database at [The National Library of Medicine](https://www.ncbi.nlm.nih.gov/), returning two csv files. Feel free to use either the script or notebook with instructions for use below.

## Input:

jupyter notebook: 
```
-- use a list of rs###'s for each snp
-- change the output filename when instantiating 
        GeneOutputFile = './nearestGenes.csv'
        SNPInfoOutputFile = './SNPInfo.csv' 
```        
command line script: 
```
-- use space separated rs###'s for each as an argument
-- the output file path is an argument. However, if you would like to change the file name attached to that path, do so on line 127 and 128 when files are created
```
## Output:

Output filepaths are parameters you provide.
```
default for nearest gene file =  current directory './nearestGenes.csv' 
default for comprehensive SNP file = current directory './SNPInfo'

```
Nearest gene information (as per ncbi)
```
GeneInfo schema:

    GENE_ID : int ID#
    NAME : string common name 
    ACC : string NC_000... 
    snpID : string rs###
```    

Comprehensive information about each SNP

```
SNPInfo schema:

    ACC : string         
    ALLELE: Bool 
    ALLELE_ORIGIN: string : rs###, however most of the time is blank 
    CHR : string chromosome snp found on 
    CHRPOS : string chr:bp
    CHRPOS_PREV_ASSM : string chr:bp from previous assembly
    CHRPOS_SORT : string chromosome position but with 000 in front
    CITED_SORT : string string
    CLINICAL_SIGNIFICANCE : string
    CLINICAL_SORT : string 
    CREATEDATE : data record of snp was created
    DOCSUM : string summary of all columns
    FXN_CLASS : string comma separate (i.e. 'intron_variant,genic_downstream_transcript_variant')
    GENES : list of dictionaries of genes (i.e. [{'NAME': 'LRRC8D', 'GENE_ID': '55144'}])
    GLOBAL_MAFS : list of dictionary (i.e. {'STUDY': '1000Genomes', 'FREQ': 'C=0.049521/.....)
    GLOBAL_POPULATION : string
    GLOBAL_SAMPLESIZE : string
    HANDLE : string reported biobanks with snp? (i.e. 1000GENOMES,EVA_UK10K_TWINSUK,....)
    IDList : list of uid's generated from first get_idList(snp)
    MERGED_SORT : string
    ORIG_BUILD : string
    SNP_CLASS : string (i.e. 'snv')
    SNP_ID : string
    SNP_ID_SORT : string
    SPDI : string (i.e. 'NC_000005.10:95305173:T:C')
    SS : string (i.e. '10263196,82343251,82652005,112205736,165508732')
    SUSPECTED : string
    TAX_ID : string
    TEXT : string (i.e. 'MergedRs=6896334')
    UPDATEDATE : Date
    UPD_BUILD : string
    VALIDATED : string (i.e. 'by-frequency,by-alfa,by-cluster')
    snpID : string rs###
```    
## Installations you will need

```
from Bio import Entrez
from collections import OrderedDict
import pandas as pd

```


# Running snp-to-gene-notebook.ipynb 

You will need:

```
email address : ncbi requires this but you don't need an account. They like to keep track of who is using their databases
file path for snps : or you can create a list and name it snp_list so the functions work
output file path : the default will be your current directory with nearestGenes.csv and SNPInfo.csv as file names.

```
Run the function:

```

SNPInfo, GeneInfo = get_closest_genes(snp_list)

```
   
# Running the command line script snp_to_gene.py

cd into the folder with snp_to_gene.py and run below
note: the file path is used for both files so don't include a filename, only up to folder you would like to store them in. If you would like to change the filename from the default nearestGenes.csv and SNPInfo.csv you can do that on lines 127 and 128: 

```
python snp_to_gene.py --h for help

python snp_to_gene.py --snplist rs11955986 rs17498135 --email your@email.com --filepath any/path/you/choose
