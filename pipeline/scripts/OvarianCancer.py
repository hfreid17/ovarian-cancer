#################################################################
#################################################################
############### Ovarian Cancer Analysis - Python Support ############
#################################################################
#################################################################
##### Author: Hannah Freid
##### Affiliation: Ma'ayan Laboratory,
##### Icahn School of Medicine at Mount Sinai

#############################################
########## 1. Load libraries
#############################################
##### 1. Python modules #####
from ruffus import *
import sys
import pandas as pd
import numpy as np
import scipy.stats as ss
import json

##### 2. Custom modules #####
# Pipeline running

#############################################
########## 2. General Setup
#############################################
##### 1. Variables #####

# TCGA expression normalized (fpkm) data file
fpkm_file = "rawdata/TCGA-OV-fpkm-uq.txt"

########################################################
########################################################
########## S1. Finding Gene Signatures for Each Sample
########################################################
########################################################

##### TCGA-fpkm expression data that is downloaded is made into a dataframe and this is used to find gene signatures for each sample.
### Input: TCGA-fpkm expression data (.txt file)
### Output: Dictionary of the top 500 and bottom 500 genes for each sample


############################################################
########## 1. Reading in data and Z-scores across rows and columns of dataframe
############################################################

# @follows(mkdir("s1-expression.dir")) ## may add subdirectories later

@transform(fpkm_file,
            suffix(".txt"),
            "_zscore.txt")

def zscore(infile, outfile):


    # Read the infile using pandas
    df_expfpkm = pd.read_table(infile).set_index('gene_symbol')

    # perform the z-score normalization (axis=0 is along rows)
    df_expfpkm_zscore = df_expfpkm.T.apply(ss.zscore, axis=0).T.dropna().apply(ss.zscore, axis=0)

    # Save the Z-score normalized df to outfile
    df_expfpkm_zscore.to_csv(outfile, sep='\t')

    
############################################################
########## 2. Creating dictionary of gene signatures
############################################################

@transform(zscore, 
            suffix('_zscore.txt'), 
            "_signatures.json")


def signatures(infile, outfile):
    
    # create empty dictionary
    dict_signatures = {}

    df_expfpkm_zscore = pd.read_table(infile).set_index('gene_symbol')
    
    
    # Get number of samples
    n = len(df_expfpkm_zscore.columns)
    i = 0
    
    # Loop through samples 
    for sample in df_expfpkm_zscore.columns[:]:

        # Print status
        i += 1
        print('Doing sample {sample} ({i}/{n})...'.format(**locals()))

         # Extract column
        col = df_expfpkm_zscore[sample].sort_values(ascending=False) # sort values in decreasing order
        
        genesets = {
            "top":col.index[:500].tolist(),
            "bottom":col.index[-500:].tolist()
            
            }
            
        # Extract the top 500 genes from the index
        dict_signatures[sample] = genesets

    # open file
    f = open(outfile, 'w')

    # write to file
    f.write(json.dumps(dict_signatures, ensure_ascii=False, indent=4))     
    

    # close file
    f.close()




######################################################################
######################################################################
########## S3. Find L1000fwd Reversing and Mimicking Drugs-result ids
######################################################################
######################################################################

##### Gene signatures for each sample is used on L1000fwd to find result ids for samples for reversing and mimicking drugs for signatures.
### Input: Dictionary of the top 500 and bottom 500 genes for each sample
### Output: Dictionary of the L1000fwd reversing drug result ids for each sample

############################################################
########## 1. Reversing drug result ids
############################################################

@transform(signatures,
            suffix('_signatures.json'),
            "_L1000reversingresultids.json")


def L1000reversingresultids(infile, outfile):


    
    # query the L1000 signatures using gene sets
    import json, requests
    from pprint import pprint

    L1000FWD_URL = 'http://amp.pharm.mssm.edu/L1000FWD/'

    # create empty dictionary for L1000FWD result ids
    l1000fwd_reversing_resultids = {}

    with open(infile) as infile2:
        signature_data = json.load(infile2)

    # Loop through dictionary
    for sample, genesets in signature_data.items():

        payload = {
            # put downregulated genes in the upregulated spot to find drugs that reverse the signature
            # same with upregulated gewnes in the downregulated spot
            'up_genes': genesets['bottom'],
            'down_genes': genesets['top']
        }
    
        response = requests.post(L1000FWD_URL + 'sig_search', json=payload)
        
        if response.status_code == 200:
            l1000fwd_reversing_resultids[sample] = response.json() #saving results to a dictionary
            json.dump(response.json(), open('api3_result.json', 'w'), indent=4)      


     # open file
    f = open(outfile, 'w')

    # write to file
    f.write(json.dumps(l1000fwd_reversing_resultids, ensure_ascii=False, indent=4))

    # close file
    f.close()


############################################################
########## 2. Mimicking drug result ids
############################################################


@transform(signatures,
            suffix('_signatures.json'),
            "_L1000mimickingresultids.json")


def L1000mimickingresultids(infile, outfile):
    
    # query the L1000 signatures using gene sets


    import json, requests
    from pprint import pprint

    L1000FWD_URL = 'http://amp.pharm.mssm.edu/L1000FWD/'

    # create empty dictionary for L1000FWD result ids
    l1000fwd_mimicking_resultids = {}

    with open(infile) as infile2:
        signature_data = json.load(infile2)

    # Loop through dictionary
    for sample, genesets in signature_data.items():

        payload = {
            # put downregulated genes in the upregulated spot to find drugs that reverse the signature
            # same with upregulated gewnes in the downregulated spot
            'up_genes': genesets['top'],
            'down_genes': genesets['bottom']
        }
    
        response = requests.post(L1000FWD_URL + 'sig_search', json=payload)
        
        if response.status_code == 200:
            l1000fwd_mimicking_resultids[sample] = response.json() #saving results to a dictionary
            json.dump(response.json(), open('api3_result.json', 'w'), indent=4)      


     # open file
    f = open(outfile, 'w')

    # write to file
    f.write(json.dumps(l1000fwd_mimicking_resultids, ensure_ascii=False, indent=4))

    # close file
    f.close()



#######################################################
#######################################################
########## S4. L1000fwd Retrieve Top 50 Results
#######################################################
#######################################################

##### Takes L1000fwd result ids for each sample and gives top 50 signature ids. Does this for both the reversing and mimicking drugs.
### Input: Dictionary of L1000fwd result ids for each sample
### Output: Dictionary of the L1000fwd top 50 signature ids for each sample.

########################################################################
########## 1. Top 50 signature ids for reversing drugs for each sample
#######################################################################


@transform(L1000reversingresultids,
            suffix('_L1000reversingresultids.json'),
            "_L1000reversingtop50ids.json")


def L1000reversingtop50(infile, outfile):
    
    import json, requests
    from pprint import pprint

    with open(infile) as infile2:
        l1000reversingresultids2 = json.load(infile2)
    
    
    # create empty dictionary
    l1000_reversing_top50 = {}
    
    
    L1000FWD_URL = 'http://amp.pharm.mssm.edu/L1000FWD/'
    
    
    for k in l1000reversingresultids2:
        result_id = str(l1000reversingresultids2[k])
        result_id_edited = result_id[14:40]
        
        response = requests.get(L1000FWD_URL + 'result/topn/' + result_id_edited)
        
        if response.status_code == 200:
            l1000_reversing_top50[k] = response.json() #adding top 50 results to dictionary
            json.dump(response.json(), open('api4_result.json', 'w'), indent=4)


    # open file
    f = open(outfile, 'w')

    # write to file
    f.write(json.dumps(l1000_reversing_top50, ensure_ascii=False, indent=4))

    # close file
    f.close()




########################################################################
########## 2. Top 50 signature ids for mimicking drugs for each sample
#######################################################################
@transform(L1000reversingresultids,
            suffix('_L1000mimickingresultids.json'),
            "_L1000mimickingtop50ids.json")


def L1000mimickingtop50(infile, outfile):
    
    import json, requests
    from pprint import pprint

    with open(infile) as infile2:
        l1000mimickingresultids2 = json.load(infile2)
    
    
    # create empty dictionary
    l1000_mimicking_top50 = {}
    
    
    L1000FWD_URL = 'http://amp.pharm.mssm.edu/L1000FWD/'
    
    
    for k in l1000mimickingresultids2:
        result_id = str(l1000mimickingresultids2[k])
        result_id_edited = result_id[14:40]
        
        response = requests.get(L1000FWD_URL + 'result/topn/' + result_id_edited)
        
        if response.status_code == 200:
            l1000_mimicking_top50[k] = response.json() #adding top 50 results to dictionary
            json.dump(response.json(), open('api4_result.json', 'w'), indent=4)


    # open file
    f = open(outfile, 'w')

    # write to file
    f.write(json.dumps(l1000_mimcking_top50, ensure_ascii=False, indent=4))

    # close file
    f.close()




#######################################################
#######################################################
########## S5. L1000fwd Get Single Signature by ID
#######################################################
#######################################################

##### Takes L1000fwd top 50 signature ids for each sample and ouputs signature--use this to extract drug.
### Input: Dictionary of the L1000fwd top 50 signatue ids for each sample.
### Output: Dataframe of the 50 top drugs for each sample--index is drug

##############################################################
########## 1. Getting single singature by id--reversing drugs
##############################################################

@transform(L1000reversingtop50,
            suffix('_L1000reversingtop50ids.json'),
            "_L1000reversingdrugsignatureids.json")


def L1000reversingdrugsignatures(infile, outfile):
    
    import json, requests
    from pprint import pprint

    with open(infile) as infile2:
        l1000reversingtop50ids2 = json.load(infile2)

    # create empty dictionary
    l1000_reversing_signatureids = {}
    
    
    L1000FWD_URL = 'http://amp.pharm.mssm.edu/L1000FWD/'


    for sample in l1000reversingtop50ids2:
        sig_id = l1000reversingtop50ids2[sample] # might need to take part of the string/might be a nested dictionary
        response = requests.get(L1000FWD_URL + 'sig/' + sig_id)
        if response.status_code == 200:
            l1000_reversing_signatureids = response.json()
            json.dump(response.json(), open('api2_result.json', 'wb'), indent=4)

    # open file
    f = open(outfile, 'w')

    # write to file
    f.write(json.dumps(l1000_reversing_signatureids, ensure_ascii=False, indent=4))

    # close file
    f.close()



##############################################################
########## 2. Getting perturbation/drug--reversing drugs
##############################################################




##############################################################
########## 3. Getting single singature by id--mimicking drugs
##############################################################

@transform(L1000mimickingtop50,
            suffix('_L1000mimickingtop50ids.json'),
            "_L1000mimickingdrugsignatureids.json")


def L1000mimickingdrugsignatures(infile, outfile):
    
    import json, requests
    from pprint import pprint

    with open(infile) as infile2:
        l1000mimickingtop50ids2 = json.load(infile2)

    # create empty dictionary
    l1000_mimicking_signatureids = {}
    
    
    L1000FWD_URL = 'http://amp.pharm.mssm.edu/L1000FWD/'


    for sample in l1000mimickingtop50ids2:
        sig_id = l1000mimickingtop50ids2[sample] # might need to take part of the string/might be a nested dictionary
        response = requests.get(L1000FWD_URL + 'sig/' + sig_id)
        if response.status_code == 200:
            l1000_mimicking_signatureids = response.json()
            json.dump(response.json(), open('api2_result.json', 'wb'), indent=4)

    # open file
    f = open(outfile, 'w')

    # write to file
    f.write(json.dumps(l1000_mimicking_signatureids, ensure_ascii=False, indent=4))

    # close file
    f.close()




##############################################################
########## 4. Getting perturbation/drug--mimicking drugs
##############################################################