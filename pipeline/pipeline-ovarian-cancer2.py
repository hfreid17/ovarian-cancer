#################################################################
#################################################################
############### Ovarian Cancer Analysis ############
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
import glob
import requests


##### 2. Custom modules #####
# Pipeline running
sys.path.append('pipeline/scripts')
import OvarianCancer as P

#############################################
########## 2. General Setup
#############################################
##### 1. Variables #####

# TCGA expression normalized (fpkm) data file
fpkm_file = "rawdata/TCGA-OV-fpkm-uq.txt"


##### 2. R Connection #####
from rpy2.robjects import r, pandas2ri
pandas2ri.activate()
r.source('pipeline/scripts/ovarian-cancer.R')



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




###################################################################################
###################################################################################
########## S2. Find L1000fwd Reversing Drugs-result ids--alone drugs
###################################################################################
###################################################################################

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



#############################################################
#############################################################
########## S3. L1000fwd Retrieve Top 50 Results--alone drugs
############################################################
############################################################

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






#############################################################################
#############################################################################
########## S4. L1000fwd Get Single Signature by Signature ID--alone drugs
#############################################################################
#############################################################################

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






##################################################
##################################################
########## S5. Finding MOA of Drugs--alone drugs
#################################################
#################################################

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




#####################################################################################################
#####################################################################################################
########## S6. Merging Dataframes of Drug Signatures and Signature Ids--alone drugs
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

    ## getting just reversing alone drugs (similar)
    df_merged2 = df_merged[df_merged["direction"]=="similar"].copy(deep=True)


    df_merged2.to_csv(outfile, sep="\t")





#########################################
#########################################
########## S7. Plotting--alone drugs
#########################################
#########################################


### I DID PLOTTING ON A JUPYTER NOTEBOOK



##########################################
##########################################
########## S8. Chea3 Queries
##########################################
##########################################

# ##### Takes z-scored expression data and queries chea3 (transformation into TF).
# ### Input: Zscored expression data.
# ### Output: Separate files of TF from chea3 queries for each sample.

# ###############################
# ########## 1. Query Chea3
# ###############################

@mkdir("rawdata/chea.dir")

@subdivide(zscore,
        formatter(),
        'rawdata/chea.dir/*_chea.txt',
        'rawdata/chea.dir/')

def runChea(infile, outfiles, outfile_root):

    # making dictionary of top 300 most differentially expressed genes for each sample
    df_zscore = pd.read_table(infile).set_index("gene_symbol")

    df_zscore_abs = df_zscore.abs()

    dict_signatures_300 = {}

    for sample in df_zscore.columns:
        col = df_zscore_abs[sample].sort_values(ascending=False) # sort values in decreasing order
        dict_signatures_300[sample] = col.index[:300]
    

    # looping through samples to query chea3
    for tcga_sample, geneset_forchea in dict_signatures_300.items():
        print('Doing {tcga_sample}...'.format(**locals()))
        outfile = '{outfile_root}{tcga_sample}_chea.txt'.format(**locals())
        chea_results_r = r.Rrun_chea(geneset_forchea, tcga_sample, library="ReMap")
        chea_results_df = pandas2ri.ri2py(chea_results_r)
        chea_results_df.to_csv(outfile, sep = "\t")

    




#########################################
#########################################
########## S9. Chea3 Analysis
#########################################
#########################################

# ##### 
# ### Input:
# ### Output: 

# #########################################
# ########## 1. Converting to TF space
# #########################################


# @merge(runChea, 
#         glob.glob("rawdata/chea.dir/*"),
#         "rawdata/tfspace.txt")

@merge(runChea,
        "rawdata/tfspace.txt")


def transformation_toTF(infiles, outfile):


    # make each df, use pivot, add to list of dfs, concatenate dfs in for loop

    list_dfs = []

    for infile in infiles:
        temp_df = pd.read_table(infile)
        temp_df2 = pd.DataFrame()
        temp_df2["sample_id"] = temp_df["set1"]
        temp_df2["TF"] = temp_df["TF"]
        temp_df2["p_val"] = temp_df["FET.p.val"]
        transformed_df = temp_df2.pivot(index = "sample_id", columns="TF", values = "p_val")
        list_dfs.append(transformed_df)

    df_final = pd.DataFrame()


    for df in list_dfs:
        df_final = pd.concat([df_final, df])
    
    df_final.to_csv(outfile, sep="\t")
    



    
######################################################################
######################################################################
########## S10. Finding Combination Drugs--Updating Gene Signatures
#####################################################################
#####################################################################

# ##### Takes list of currently used ovarian cancer drugs (from Mt. Sinai EMR data) and updates gene signature for each sample based on genes affected/regulated by the currently used drugs.
# ### Input: List of currently used drugs, ovarian cancer patient's gene signatures
# ### Output: List of updated gene signatures for each sample (genes affected by currently used drugs are removed)

# #############################################################################
# ########## 1. Creating updated gene signatures for each sample for each drug
# #############################################################################

@mkdir("rawdata/combinationdrugs.dir")

@subdivide(signatures,
        formatter(),
        'rawdata/combinationdrugs.dir/*_updatedsignatures.json',
        'rawdata/combinationdrugs.dir/')

# this was from when I reran this step to get the gene signatures for each drug
# @subdivide(signatures,
#         formatter(),
#         'rawdata/combinationdrugs.dir/*_signature_foreachdrug_forplots3.json',
#         'rawdata/combinationdrugs.dir/')


def pairdrugs_updatedsignatures(infile, outfiles, outfile_root):

    import json, requests
    from pprint import pprint

    list_currentovdrugs = ["paclitaxel", "gemcitabine", "cisplatin", "doxorubicin", "topotecan", "docetaxel", "etoposide", "cyclophosphamide", "ifosfamide", "irinotecan", "pemetrexed", "tamoxifen"]
    
    df_sigmetadata = pd.read_csv("/Users/maayanlab/Documents/Ovarian Cancer Project/ovarian-cancer/CD_signature_metadata.csv").set_index("sig_id")
    
    with open(infile) as infile2:
        genesignatures = json.load(infile2)
    
    df_genesig = pd.DataFrame(genesignatures)

    df_genesig = df_genesig.T

    L1000FWD_URL = 'http://amp.pharm.mssm.edu/L1000FWD/'

    dict_updated_genesig = {}
    
    for drug in list_currentovdrugs:

        # finding pertids for current ov cancer drugs to be used to find updated signature ids in metadata file

        response = requests.get(L1000FWD_URL + 'synonyms/' + drug)
        if response.status_code == 200:
            dict_current_pertids = response.json()
            json.dump(response.json(), open('api1_result.json', 'w'), indent=4)

    

        # finding updated signature ids from pertids for ov cancer drugs in metadata file

        df_current_pertids = pd.DataFrame(dict_current_pertids)

        list_current_pertids = list(df_current_pertids["pert_id"])

        list_current_sigid_updated = []


        for sig, rowdata in df_sigmetadata.iterrows():
            temp_pertid = rowdata["pert_id"]
            if temp_pertid in list_current_pertids:
                list_current_sigid_updated.append(sig)

    
        # getting drug signature from l1000

        dict_current_signatures = {}

        for item in list_current_sigid_updated:
            response = requests.get(L1000FWD_URL + 'sig/' + item)
            if response.status_code == 200:
                dict_current_signatures[item] = (response.json())
                json.dump(response.json(), open('api2_result.json', 'w'), indent=4)

        
        # converting nested dictionary of drug signatures to df

        df_current_signatures = pd.DataFrame(dict_current_signatures).T



        # extracting up and down genes for drug from df and for patients from df to be used for subtracting from each patient's signature

        current_down = []
        current_up = []

        for sigid, rowdata in df_current_signatures.iterrows():
            current_down_temp = rowdata["down_genes"]
            current_up_temp = rowdata["up_genes"]
            current_down = current_down + current_down_temp
            current_up = current_up + current_up_temp

    # ## not from original pipeline--put in after for analysis--if this step ends here and file is written to output, can get the drug induced gene signature for each drug

    #     dict_3 = {}
    #     dict_3["up"] = current_up
    #     dict_3["down"] = current_down
    #     outfile = '{outfile_root}{drug}_signature_foreachdrugs_forplots3.json'.format(**locals())

    #     print(outfile)

    #     # open file
    #     f = open(outfile, 'w')

    #     # write to file
    #     f.write(json.dumps(dict_3, ensure_ascii=False, indent=4))

    #     # close file
    #     f.close() 
        

        
        for sample, rowdata in df_genesig.iterrows():
            sample_bottom_temp = rowdata["bottom"]
            sample_top_temp = rowdata["top"]

            # subtracting drug genes from sample genes
            for gene in sample_bottom_temp:
                if gene in current_up:
                    sample_bottom_temp.remove(gene)
            for gene2 in sample_top_temp:
                if gene2 in current_down:
                    sample_top_temp.remove(gene2)
            
            updated_genesets = {
                "top":sample_top_temp,
                "bottom":sample_bottom_temp
            }
            
            # Extract the top 500 genes from the index
            dict_updated_genesig[sample] = updated_genesets
        
        outfile = '{outfile_root}{drug}_updatedsignatures.json'.format(**locals())

        # open file
        f = open(outfile, 'w')

        # write to file
        f.write(json.dumps(dict_updated_genesig, ensure_ascii=False, indent=4))

        # close file
        f.close()

    



########################################################################
########################################################################
########## S11. Finding Combination Drugs from Updated Gene Signatures
#######################################################################
#######################################################################

# ##### Takes updated gene signatures and predicts drug pairs for each of these drugs for each patient based on the patient's gene signature.
# ### Inputs: Updated list of gene signatures (for each drug)
# ### Output: List of predicted drugs to be paired with each currently used drug for each patient.

# ##################################################
# ########## 1. Query L1000 to find reversing drugs
# #################################################



@transform(pairdrugs_updatedsignatures,
            suffix("_updatedsignatures.json"),
            "_resultids.json")


def pairdrugs_queryL1000(infile, outfile):

    print(outfile)

    # query the L1000 signatures using gene sets

    L1000FWD_URL = 'http://amp.pharm.mssm.edu/L1000FWD/'


    # create empty dictionary for L1000FWD result ids
    l1000_pair_reversing_resultids = {}

    with open(infile) as infile2:
        updatedsig = json.load(infile2)

    # Loop through dictionary
    for sample, genesets in updatedsig.items():

        payload = {
            # put downregulated genes in the upregulated spot to find drugs that reverse the signature
            # same with upregulated gewnes in the downregulated spot
            'up_genes': genesets['bottom'],
            'down_genes': genesets['top']
        }
        
        response = requests.post(L1000FWD_URL + 'sig_search', json=payload)
            
        if response.status_code == 200:
            l1000_pair_reversing_resultids[sample] = response.json() #saving results to a dictionary
            json.dump(response.json(), open('api3_result.json', 'w'), indent=4)      

    # open file
    f = open(outfile, 'w')

    # write to file
    f.write(json.dumps(l1000_pair_reversing_resultids, ensure_ascii=False, indent=4))

    # close file
    f.close()



#####################################################################################
######### 2. Top 50 signature ids for reversing drugs for each drug for each sample
####################################################################################


@transform(pairdrugs_queryL1000,
            suffix("_resultids.json"),
            "_top50results.json")


def pairdrugs_top50(infile, outfile):

    print(outfile)

    with open(infile) as infile2:
        drugresults = json.load(infile2)
    
    
    # create empty dictionary
    l1000_pair_reversing_top50 = {}
    
    
    L1000FWD_URL = 'http://amp.pharm.mssm.edu/L1000FWD/'
    
    
    
    for k in drugresults:
        result_id = str(drugresults[k])
        result_id_edited = result_id[15:39] # getting just result id from string
        
        response = requests.get(L1000FWD_URL + 'result/topn/' + result_id_edited)
        
        if response.status_code == 200:
            l1000_pair_reversing_top50[k] = response.json() # adding top 50 results to dictionary
            json.dump(response.json(), open('api4_result.json', 'w'), indent=4)


    # open file
    f = open(outfile, 'w')

    # write to file
    f.write(json.dumps(l1000_pair_reversing_top50, ensure_ascii=False, indent=4))

    # close file
    f.close()




######################################################################
########## 3. Making top 50 results nested dictionary into dataframe
######################################################################


@transform(pairdrugs_top50,
            suffix('_top50results.json'),
            "_top50results_todf.txt")


def pairdrugs_top50_todf(infile, outfile):
    
    print(outfile)

    with open(infile) as infile2:
        top50_dict = json.load(infile2)


    temp_df_signatureids = pd.DataFrame(top50_dict).T

    df_list = []

    for sample, rowdata in temp_df_signatureids.iterrows(): # row id is sample, other data is rowdata
        for direction in ["similar"]:
            df = pd.DataFrame(rowdata[direction])
            df["sample"] = sample
            df_list.append(df)
    df_signatureids = pd.concat(df_list)

    df_signatureids2 = df_signatureids.set_index("sig_id")
    
    df_signatureids2.to_csv(outfile, sep="\t")




##############################################################
########## 4. Getting single singature by id
##############################################################

@transform(pairdrugs_top50_todf,
            suffix('_top50results_todf.txt'),
            "_predicteddrugsignature.json")


def pairdrugs_predicteddrugsignature(infile, outfile):
    
    print(outfile)

    L1000FWD_URL = 'http://amp.pharm.mssm.edu/L1000FWD/'

    df_top50 = pd.read_table(infile)
    
    list_sigids2 = list(df_top50["sig_id"])

    dict_pair_predicteddrugsignatures = {}

    for signature_id in list_sigids2:
        response = requests.get(L1000FWD_URL + 'sig/' + signature_id)
        if response.status_code == 200:
            dict_pair_predicteddrugsignatures[signature_id] = (response.json())
#       json.dump(response.json(), open('api2_result.json', 'w'), indent=4)


    # open file
    f = open(outfile, 'w')

    # write to file
    f.write(json.dumps(dict_pair_predicteddrugsignatures, ensure_ascii=False, indent=4))

    # close file
    f.close()





########################################################################
########################################################################
########## S12. Combination Drugs Analysis
#######################################################################
#######################################################################

# ##### Takes dictionary of pair drug signatures for each current drug and makes dataframe of drugs suggested to be paired with each drug for each patient.
# ### Inputs: Updated list of gene signatures (for each drug)
# ### Output: List of suggested drugs to be paired with each currently used drug for each patient.

# ##################################################################################################
# ########## 1. Make dictionary into dataframe with fields of interest and pvals from top50 results
# ##################################################################################################

@transform(pairdrugs_predicteddrugsignature,
            regex(r'(.*)_predicteddrugsignature.json'),
            add_inputs(r'\1_top50results_todf.txt'),
            r"\1_df_foranalysis.txt")


def pairdrugs_df_analysis(infiles, outfile):

    signaturefile, top50file = infiles

    print(outfile)

    with open(signaturefile) as infile2:
        pair_drugsignatures_dict = json.load(infile2)

    df_pair_drugsignatures = pd.DataFrame(pair_drugsignatures_dict).T.set_index("sig_id")

    df_pair_signatureids = pd.read_table(top50file).set_index("sig_id")

    df_pair_foranalysis = pd.DataFrame()

    list_sigids = df_pair_drugsignatures.index.tolist()
    list_pertdesc = (df_pair_drugsignatures["pert_desc"]).tolist()
    list_pertdose = (df_pair_drugsignatures["pert_dose"]).tolist()
    list_pertid = (df_pair_drugsignatures["pert_id"]).tolist()
    list_perttime = (df_pair_drugsignatures["pert_time"]).tolist()
    df_pair_foranalysis["sig_id"] = list_sigids
    df_pair_foranalysis["pert_desc"] = list_pertdesc
    df_pair_foranalysis["pert_dose"] = list_pertdose
    df_pair_foranalysis["pert_id"] = list_pertid
    df_pair_foranalysis["pert_time"] = list_perttime
    df_pair_foranalysis.set_index("sig_id")

    df_pair_foranalysis2 = df_pair_foranalysis.merge(df_pair_signatureids, how="inner", on="sig_id")

    df_pair_foranalysis2 = df_pair_foranalysis2.set_index("sig_id")
    
    
    df_pair_foranalysis2.to_csv(outfile, sep="\t")






# ##################################################################################################
# ########## 2. Correcting pertnames using metadata file
# ##################################################################################################   


## I DID NOT RUN THIS STEP WITH THE CORRECT CODE BELOW IN THE PIPELINE--THE CODE I ORIGINALLY RAN FOR THIS STEP IN THE PIPELINE WAS
## WRONG AND I RAN A CORRECTED VERSION (AS SEEN BELOW) ON A JUPYTER NOTEBOOK--THE CORRECTED FILES ARE DRUGNAME_WORSTPROG_COUNTS2.TXT (THESE FILES INCLUDE OTHER STEPS THAT COME AFTER THIS STEP TOO)

@transform(pairdrugs_df_analysis,
            suffix('_df_foranalysis.txt'),
            "_df_correctpertnames_test.txt")


def pairdrugs_correct_pertnames(infile, outfile):

    print(outfile)

    df_pair = pd.read_table(infile)


    ## to be able to use novel drug data and not have same drug listed under multiple brd ids

    df_drugsmetadata = pd.read_csv("/Users/maayanlab/Documents/Ovarian Cancer Project/ovarian-cancer/rawdata/Drugs_metadata.csv").set_index("pert_iname")
    df_merge_grouped = df_merge.groupby("pert_iname").size().rename("counts").to_frame().sort_values(by="counts", axis=0, ascending=False)

    
    # save to file
    df_test3.to_csv(outfile, sep="\t", index = False)






# ##################################################################################################
# ########## 3. Getting FDA approval info and MOA using metadata file made by zichen
# ##################################################################################################   

# @transform(pairdrugs_correct_pertnames,
#             suffix('_df_correctpertnames_counts_sorted.txt'),
#             "_df_final_analysis2.txt")


# def pairdrugs_finaldf_analysis(infile, outfile):

#     print(outfile)

#     df_pair_correctnames2 = pd.read_table(infile).set_index("pert_names")


#     ## to be able to use novel drug data and not have same drug listed under multiple brd ids

#     df_fdamoa_data = pd.read_csv("/Users/maayanlab/Documents/Ovarian Cancer Project/ovarian-cancer/L1000FWD_drugs_7-30-18.csv")



#     df_pair_correctnames2_merge = df_pair_correctnames2.merge(df_fdamoa_data, left_on="pert_names", right_on="pert_iname").drop_duplicates("pert_iname")

#     # finding pert_name from pertids for predicted drugs using metadata file --will give name or pertid if no name available

#     # list_pair_pertnames = list(df_pair_correctnames2.index)

#     # list_pair_MOA = []
#     # list_pair_FDA = []
#     # list_pair_target = []


#     # for name, rowdata in df_fdamoa_data.iterrows():
#     #     temp_MOA = rowdata["MOA"]
#     #     temp_FDA = rowdata["Phase"]
#     #     temp_target = rowdata["Target"]
#     #     if name in list_pair_pertnames:
#     #         list_pair_MOA.append(temp_MOA)
#     #         list_pair_FDA.append(temp_FDA)
#     #         list_pair_target.append(temp_target)
#     #         list_pair_pertnames.remove(name)

#     # df_pair_correctnames2["MOA"] = list_pair_MOA
#     # df_pair_correctnames2["Phase"] = list_pair_FDA
#     # df_pair_correctnames2["target"] = list_pair_target
#     # # print(len(list_pair_MOA))
#     # # print(len(list_pair_FDA))
#     # # print(len(list_pair_target))
#     # # print(len(df_pair_correctnames2.index.values))


#     # save to file
#     df_pair_correctnames2_merge.to_csv(outfile, sep="\t")















##################################################################################
#################################################################################
########## S13. Finding Combination and Alone Drugs for Worst Prognosis Patients--## I RAN THESE STEPS ON A JUPYTER NOTEBOOK AT THE END AND THEN PUT IT IN THE PIPELINE BELOW--BUT I HAD NOT RUN THIS PART OF THE PIPELINE 
#################################################################################
#################################################################################

# ##### Takes dictionary of pair drug signatures for each current drug and ???
# ### Inputs: Updated list of gene signatures (for each drug)
# ### Output: List of predicted drugs to be paired with each currently used drug for each patient.

# ##################################################################################################
# ########## 1. Finding Worst Prognosis Patients
# ##################################################################################################

@transform(pairdrugs_correct_pertnames,
            '/Users/maayanlab/Downloads/clinical.cart.2018-06-28/clinical.csv.txt',
            add_inputs("/Users/maayanlab/Documents/Ovarian Cancer Project/TCGA-OV-fpkm-uq.txt")
            "/rawdata/clinical_prognosis_worst.txt")

def findingworstprogpatients(infiles, outfile):


    clinical, exp = infiles

    ## loading clinical TCGA data

    df_clinical = pd.read_csv(clinical).set_index('submitter_id')
    clinical_ids = df_clinical.index.values.tolist()


    # finding overlapping samples in clinical and expression data

    # loading expression data to get list of overlapping samples
    df_expression = pd.read_table(exp).set_index('gene_symbol')


    # making list of patient (submitter) ids to use for the clinical analysis

    sub_ids_exp = list(df_expression.columns.values)
    sub_ids_corrected = []
    for item in sub_ids_exp:
        new_item = item[0:12]
        sub_ids_corrected.append(new_item)

    overlapping_ids = []
    for item2 in sub_ids_corrected:
        if item2 in clinical_ids:
            overlapping_ids.append(item2)


    # of these overlapping patients--find group of lowest survival/worst prognosis (dying within two years)

    df_temp = pd.DataFrame()
    df_temp["submitter_id"] = overlapping_ids
    df_clinical_cut = df_clinical.merge(df_temp, left_on="submitter_id", right_on="submitter_id").set_index("submitter_id")
    df_clinical_prognosis = pd.DataFrame()
    list_days_to_death = df_clinical_cut["days_to_death"].tolist()
    df_clinical_prognosis["submitter_id"] = df_clinical_cut.index
    df_clinical_prognosis["days_to_death"] = list_days_to_death

    for number, rowdata in df_clinical_prognosis.iterrows():
        if rowdata["days_to_death"] == "--":
            df_clinical_prognosis = df_clinical_prognosis.drop(labels=number, axis=0)


    days_to_death_updated = df_clinical_prognosis["days_to_death"].tolist()
    days_to_death_updated_int = []
    for item in days_to_death_updated:
        item_int = int(item)
        days_to_death_updated_int.append(item_int)
    df_clinical_prognosis["days_to_death2"] = days_to_death_updated_int
    df_clinical_prognosis = df_clinical_prognosis.sort_values(by="days_to_death2", axis=0)
    for number, rowdata in df_clinical_prognosis.iterrows():
        if rowdata["days_to_death2"] > 730:
            df_clinical_prognosis = df_clinical_prognosis.drop(labels=number, axis=0)
    


    # saving to file
    df_clinical_prognosis.to_csv(outfile, sep="\t")




# ##################################################################################################
# ########## 2. Finding Alone Drugs for Worst Prognosis Patients
# ##################################################################################################

@transform(findingworstprogpatients,
            suffix("/rawdata/TCGA-OV-fpkm-uq_L1000_drugsig_sigids_merged.txt"),
            "/rawdata/alone_worst_prog_counts.txt")


def alonedrugs_worstprog(infile, outfile):


    ## getting alone drugs
    df_alone = pd.read_table(infile).set_index("sig_id")

    ## FIXED THIS STEP IN THE PIPELINE--SO IF PIPELINE IS RERUN THIS WOULDN'T BE NEEDED--BUT NEEDED IT WHEN I RAN THIS PART ON JUPYTER
    # ## 1--getting just reversing alone drugs (similar)
    # df_alone2 = df_alone[df_alone["direction"]=="similar"].copy(deep=True)

    ## 2--getting just worst prognosis patient alone drugs
    list_alone_samples = df_alone["sample"].tolist()
    list_alone_samples_edited = []
    for item in list_alone_samples:
        temp_item = item[0:12]
        list_alone_samples_edited.append(temp_item)
    df_alone["samples2"] = list_alone_samples_edited

    df_alone_clinical = df_alone2.merge(df_clinical_prognosis, left_on="samples2", right_on="submitter_id")


    ## 3--# to be able to use novel drug data and not have same drug listed under multiple brd ids


    # finding pert_name from pertids for predicted drugs using metadata file --will give name or pertid if no name available

    df_alone3 = pd.DataFrame()

    df_alone3["pert_name"] = list_pair_pertnames2
    
    df_alone4 = df_alone3.groupby('pert_name').size().rename('counts').to_frame()
    df_alone4 = df_alone4.sort_values(by="counts", axis=0, ascending=False)

        
    df_alone4.to_csv(outfile, sep="\t")






# ##################################################################################################
# ########## 3. Finding Combination Drugs for Worst Prognosis Patients
# ##################################################################################################

@transform(alonedrugs_worstprog,
            suffix("/combinationdrugs.dir/_df_foranalysis.txt"),
            "/rawdata/_alone_worst_prog_counts.txt")


def combinationdrugs_worstprog(infile, outfile):


    ## finding combination drugs for patients with the worst prognosis



    # FROM WHEN I RAN THIS ON JUPYTER
    # ## finding combination drugs for patients with the best prognosis
    # infiles = glob.glob('rawdata/combinationdrugs.dir/*_df_foranalysis.txt')

    df_drugsmetadata = pd.read_csv("/Users/maayanlab/Documents/Ovarian Cancer Project/ovarian-cancer/rawdata/Drugs_metadata.csv").set_index("pert_iname")

    for infile in infiles:
        df_results = pd.read_table(infile)

        list_results_samples = list(df_results["sample"])
        list_results_samples_edited = []
        for item in list_results_samples:
            temp_item = item[0:12]
            list_results_samples_edited.append(temp_item)
        df_results["samples2"] = list_results_samples_edited

        df_clinical_prognosis.read_table("/Users/maayanlab/Documents/Ovarian Cancer Project/ovarian-cancer/rawdata/clinical_prognosis_worst.txt")

        df_results_merged = df_results.merge(df_clinical_prognosis, left_on="samples2", right_on="submitter_id")
        


        df_merge2 = df_results_merged.merge(df_drugsmetadata.reset_index(), on="pert_id")
        
        df_merge3 = df_merge2.groupby("pert_iname").size().rename("counts").to_frame().sort_values(by="counts", axis=0, ascending=False)


        
        # drug_name = os.path.basename(infile).split('_')[0]
        
        # outfile = './rawdata/combinationdrugs.dir/{drug_name}_worstprog_counts2.txt'.format(**locals())

        
        df_merge3.to_csv(outfile, sep="\t")





##################################################
##################################################
########## Run pipeline
##################################################
##################################################
pipeline_run([sys.argv[-1]], multiprocess=1, verbose=1)
print('Done!')