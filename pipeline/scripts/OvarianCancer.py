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
import matplotlib.pyplot as plt
import seaborn as sns

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
### Output: Nested dictionary of the L1000fwd top 50 signature ids for each sample.

########################################################################
########## 1. Top 50 signature ids for reversing drugs for each sample
#######################################################################


@transform(L1000reversingresultids,
            suffix('_L1000reversingresultids.json'),
            "_L1000top50ids.json")


def L1000top50(infile, outfile):
    
    import json, requests
    from pprint import pprint

    with open(infile) as infile2:
        l1000reversingresultids2 = json.load(infile2)
    
    
    # create empty dictionary
    l1000_reversing_top50 = {}
    
    
    L1000FWD_URL = 'http://amp.pharm.mssm.edu/L1000FWD/'
    # L1000FWD_URL = 'http://146.203.54.243:31278/custom/SigineDMOA'
    
    
    for k in l1000reversingresultids2:
        result_id = str(l1000reversingresultids2[k])
        result_id_edited = result_id[15:39] # getting just result id from string
        
        response = requests.get(L1000FWD_URL + 'result/topn/' + result_id_edited)
        
        if response.status_code == 200:
            l1000_reversing_top50[k] = response.json() # adding top 50 results to dictionary
            json.dump(response.json(), open('api4_result.json', 'w'), indent=4)


    # open file
    f = open(outfile, 'w')

    # write to file
    f.write(json.dumps(l1000_reversing_top50, ensure_ascii=False, indent=4))

    # close file
    f.close()


######################################################################
########## 2. Making top 50 results nested dictionary into dataframe
######################################################################


@transform(L1000top50,
            suffix('_L1000top50ids.json'),
            "_L1000_df_signature_ids.txt")


def signatureids_todf(infile, outfile):
    
    import json, requests
    from pprint import pprint

    with open(infile) as infile2:
        top50_dict = json.load(infile2)


    temp_df_signatureids = pd.DataFrame(top50_dict).T

    df_list = []

    for sample, rowdata in temp_df_signatureids.iterrows(): # row id is sample, other data is rowdata
        for direction in ["similar", "opposite"]:
            df = pd.DataFrame(rowdata[direction])
            df["direction"] = direction
            df["sample"] = sample
            df_list.append(df)
    df_signatureids = pd.concat(df_list)

    df_signatureids2 = df_signatureids.set_index("sig_id")
    
    df_signatureids2.to_csv(outfile, sep="\t")






################################################################
################################################################
########## S5. L1000fwd Get Single Signature by Signature ID
################################################################
################################################################

##### Takes L1000fwd top 50 signature ids for each sample and ouputs signature--stored in nested dictionary, which is later made into a dataframe.
### Input: Dataframe of the L1000fwd top 50 signatue ids for each sample.
### Outputs: Nested dictionary of drug signatures for each sample and dataframe of drug signatures.

##############################################################
########## 1. Getting single singature by id
##############################################################

@transform(signatureids_todf,
            suffix('_L1000_df_signature_ids.txt'),
            "_L1000drugsignatures.json")


def L1000drugsignatures(infile, outfile):
    
    import json, requests
    from pprint import pprint

    L1000FWD_URL = 'http://amp.pharm.mssm.edu/L1000FWD/'

    df_sigids = pd.read_table(infile)
    
    list_sigids = list(df_sigids["sig_id"])

    dict_drugsignatures = {}

    for signature_id in list_sigids:
        response = requests.get(L1000FWD_URL + 'sig/' + signature_id)
        if response.status_code == 200:
            dict_drugsignatures[signature_id] = (response.json())
#       json.dump(response.json(), open('api2_result.json', 'w'), indent=4)


    # open file
    f = open(outfile, 'w')

    # write to file
    f.write(json.dumps(dict_drugsignatures, ensure_ascii=False, indent=4))

    # close file
    f.close()


###################################################################
########## 2. Converting signature nested dictionary to dataframe
###################################################################

@transform(L1000drugsignatures,
            suffix('_L1000drugsignatures.json'),
            "_L1000_df_drugsignature.txt")


def drugsignature_todf(infile, outfile):

    with open(infile) as infile2:
        drugsignatures_dict = json.load(infile2)

    df_drugsignatures = pd.DataFrame(drugsignatures_dict).T

    df_drugsignatures.to_csv(outfile, sep="\t")






##########################################
#########################################
########## S6. Finding MOA of Drugs
########################################
########################################

##### Merging the dataframes of signature ids and drug signatures, using signature id (sig_id) as the index.
### Input: Dataframe of drug signatures.
### Outputs: Dictionary and dtaframe of drugs and their MOAs.

##############################################################
########## 1. Searching L1000fwd for MOA--creating dictionary
##############################################################

@transform(drugsignature_todf,
            suffix('_L1000_df_drugsignature.txt'),
            "_L1000_dict_MOA.json")


def L1000MOA(infile, outfile):

    import json, requests
    from pprint import pprint
    
    df_signatures_temp = pd.read_table(infile)

    list_pertid = (df_signatures_temp["pert_id"]).tolist()
    list_pertid_short = list_pertid[0:5]

    dict_MOA = {}

    L1000FWD_URL = 'http://amp.pharm.mssm.edu/dmoa/report/'

    for item in list_pertid_short:
        id_brd = item
        response = requests.get(L1000FWD_URL + id_brd + "/json")

        if response.status_code == 200:
            dict_MOA[item] = response.json()
            json.dump(response.json(), open('api1_result.json', 'w'), indent=4)
    

    # open file
    f = open(outfile, 'w')

    # write to file
    f.write(json.dumps(dict_MOA, ensure_ascii=False, indent=4))

    # close file
    f.close()

###################################################################
########## 2. Converting MOA dictionary to dataframe
###################################################################

# @transform(L1000MOA,
#             suffix('_L1000drugsignatures.json'),
#             "_L1000_df_drugsignature.txt")


# def MOA_todf(infile, outfile):

#      with open("rawdata/TCGA-OV-fpkm-uq_L1000top50ids.json") as infile2:
#         top50_dict2 = json.load(infile2)

#     list_samples = top50_dict2.keys()

#     df_MOA = pd.DataFrame()
#     df_MOA["pert_id"] = list_pertid
#     df_MOA["MOA"] = list_MOA
#     df_MOA["sample"] = list_samples
#     df_MOA.set_index("sample")

#     df_drugsignatures = pd.DataFrame(drugsignatures_dict).T

#     df_drugsignatures.to_csv(outfile, sep="\t")






#####################################################################################################
#####################################################################################################
########## S7. Clustermap of MOAs of Drugs
#####################################################################################################
#####################################################################################################





#####################################################################################################
#####################################################################################################
########## S7. Merging Dataframes of Drug Signatures and Signature Ids.
#####################################################################################################
#####################################################################################################

##### Merging the dataframes of signature ids and drug signatures, using signature id (sig_id) as the index.
### Inputs: Dataframe of signature ids (from top 50 results) and drug signatures.
### Output: Merged dataframe of the two input dataframes with the signature id (sig_id) as the index.

######################################################################
########## 1. Merging the dataframes
######################################################################


@transform(drugsignature_todf,
            suffix('_L1000_df_drugsignature.txt'),
            add_inputs(signatureids_todf),
            "_L1000_drugsig_sigids_merged.txt")


def mergeddfs(infiles, outfile):
    print(infiles)
    # split(infiles)
    signature_file, id_file = infiles



    df_sig_ids = pd.read_table(id_file)
    df_drugsignatures = pd.read_table(signature_file)


    df_merged = df_sig_ids.merge(df_drugsignatures, on = "sig_id")

    df_merged.to_csv(outfile, sep="\t")





#########################################
#########################################
########## S8. Plotting
#########################################
#########################################

# ##### Takes dataframe of signatures and signature ids and makes various plots using seaborn.
# ### Input: Dataframe of signatures and signature ids.
# ### Output: Various plots 

# ##############################################################
# ########## 1. testing out clustermap--reversing drugs
# ###############################################################

@transform(mergeddfs,
            suffix('_L10000_drugsig_sigids_merged.csv'),
            "_L1000_clustermap.png")


def clustermaptest(infile, outfile):

    df2 = pd.read_table(infile)
    # df2["pvals"] = np.log10(df2["pvals"] + 1)

    # fig, ax = plt.subplots(figsize=(20,10))
    clustermap = sns.clustermap(np.log10(df2+1), z_score=0, cmap="RdBu_r")

    # clustermap2 = clustermap.get_figure()
    # clustermap2.savefig(outfile, dpi=400)