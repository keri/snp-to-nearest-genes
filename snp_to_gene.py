
import argparse
# defined command line options
# this also generates --help and error handling

CLI=argparse.ArgumentParser()
CLI.add_argument(
  "--snplist",  # name on the CLI - drop the `--` for positional/required parameters
  nargs="*",  # 0 or more values expected => creates a list
  type=str,
  default=['rs6896334'],  # default if nothing is provided
)

CLI.add_argument(
  "--email",
  type=str,  # any type/callable can be used here
  default='./',
)

CLI.add_argument(
  "--filepath",
  type=str,  # any type/callable can be used here
  default='./',
)

args = CLI.parse_args()

snp_list = args.snplist
email = args.email
filepath = args.filepath

# load python modules
from Bio import Entrez
from collections import OrderedDict
import pandas as pd

#initialize default parameters
Entrez.email = email # provide your email address
db = 'snp'          # set search to database
paramEutils = { 'usehistory':'Y',  "retmode":'xml',"retmax": 10} #Use Entrez search history to cache results


def get_idList(snp_id):
	'''Returns uid's for each snp to be used in downstream searches'''
	#retmode must be in xml format which creates a dict object to get the idList

	eSearch = Entrez.esearch(db=db, term=snp_id, **paramEutils)
	res = Entrez.read(eSearch)
	idList = res['IdList'] #returns a list of id's to create summaries for
	return(idList)

def parse_return_xml(dict_obj):
    '''Returns a dictionary of keys from entrez summary result using 'snp' database
    Entrez.esummary(db=db, id='rs####', retmax=retmax, retmode=retmode)'''
    
    record = {}
    # info to collect is three layers in
    first_order = dict_obj['DocumentSummarySet']
    second_order = first_order['DocumentSummary'] 
    
    #return the dictionary of summary
    record = second_order[0]
    
    return(record)

def get_closest_genes(snp_list):
    '''Calls ncbi snp database with a list of rs## (snp ids), and returns a dataframe
    with a row of information for each snp id and a dataframe for each protein if that snp has a 
    associated with it.
    schema for the snp df:
    snpid : text
    gene_name : text
    gene_id : gene_id list
    snp_id : snp_id from list
    maf_list : mafs in different populations
    pop_sample_size : population sample size
    handle : comma separated string of (I think) institutions that submitted info about snp
    acc : protein number to be used downstream
    spdi: protein number : base position of variant
    prot_class : class of variant 
    validated : how validated (frequency, cluster, etc)
    docsum : summary of all columns
    idList : uid list
    index : rsid

	schema for protein df:
	GENE_ID : #
	NAME : common gene name
	ACC : gene accession number
	index : rsid
    '''

    #instantiates a new dataframe for each list of snps
    snp_records = pd.DataFrame()
    genes = pd.DataFrame()

    #iterate over each snp to create row for dataframe
    for snp in snp_list:
        #every snp may have several ids associated with it
        idList = get_idList(snp)
        if len(idList) > 0:
            eSummary = Entrez.esummary(db=db, id=idList[0], **paramEutils)
            summary_dict = Entrez.read(eSummary)
            record = parse_return_xml(summary_dict)
            record['IDList'] = idList
            #Create new dataframe for one record with index as snp
            new_record = pd.DataFrame()
            new_record = new_record.append(record, ignore_index=True)
            new_record.index = [snp]
            #append to existing dataframe with all snps
            snp_records = snp_records.append(new_record)
            # separate gene information into separate df to simplify data
            if len(record['GENES']) > 0: #parse genes out into separate dataframe
                new_gene = pd.DataFrame(data=record['GENES'],index=range(len(record['GENES'])))
                new_gene['snpId'] = snp
                genes = genes.append(new_gene)

    snp_records['snpId'] = snp_records.index

    return(snp_records, genes)

if __name__ == '__main__':

	snp_records, gene_records = get_closest_genes(snp_list)
	snp_records.to_csv(filepath+'nearestGenes.csv')
	gene_records.to_csv(filepath+'SNPInfo.csv')

