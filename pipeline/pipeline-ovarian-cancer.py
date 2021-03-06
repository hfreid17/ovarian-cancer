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

@originate(None)
def test(infile):
    r.r_function(['this', 'is', 'alist'])





#########################################
#########################################
########## S8. Chea3 Queries
#########################################
#########################################

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
    




# @merge(runChea,
#             regex('chea2.dir/*_chea.txt'), # is this the right way to do this?
#             "TFcounts.json")

# def countTF(infiles, outfile):

#     list_TFs = []

#     # creating master list of TFs
#     for file in infiles:
#         temp_df = pd.read_table(file)
#         list_TFs.append(temp_df["TF"])
    
#     # counting occurrances of each TF
#     dict_TFcounts = {}
#     for item in list_TFs:
#         count = 0
#         for item2 in list_TFs:
#             if item == item2:
#                 count += 1
#             dict_TFcounts[item] = count

    
#     # open file
#     f = open(outfile, 'w')

#     # write to file
#     f.write(json.dumps(dict_TFcounts, ensure_ascii=False, indent=4))

#     # close file
#     f.close()


    

# # #########################################
# # ########## 2. TF survival analysis
# # #########################################

# @transform(runChea,
#             regex('chea2.dir/*_chea.txt'), # is this the right way to do this?
#             "??")

# def countTF(infiles, outfile):

#     df_tfs_sample = pd.DataFrame()
    
#     for file in infiles:
#         temp_df = pd.read_table(file)
#         list_TFs_temp = temp_df["TF"]
#         list_samplename = temp_df["set1"]
#         samplename = list_samplename[1]
#         samplename = samplename[0:12] # this may need to be 11?
#         df_tfs_sample[samplename] = list_TFs_temp

#     df_tfs_sample = df_tfs_sample.T

#     df_survival = pd.read_table(PATH??)

#     # add column label for submitter id for df_tfs_sample--call it "patient_barcode"

#     df_survival_TF = df_tfs_sample.merge(df_survival, how="inner", on="patient_barcode" )






    
######################################################################
######################################################################
########## S10. Finding Combination Drugs--Updating Gene Signatures
#####################################################################
#####################################################################

# ##### Takes list of currently used ovarian cancer drugs (from Mt. Sinai EMR data) and updates gene signature for each sample based on genes affected/regulated by the currently used drugs.
# ### Input: List of currently used drugs, ovarian cancer patient's gene signatures
# ### Output: List of updated gene signatures for each sample (genes affected by currently used drugs are removed)

# #########################################
# ########## 1. 
# #########################################

@mkdir("rawdata/combinationdrugs.dir")

# @subdivide(signatures,
#         formatter(),
#         'rawdata/combinationdrugs.dir/*_updatedsignatures.json',
#         'rawdata/combinationdrugs.dir/')

@subdivide(signatures,
        formatter(),
        'rawdata/combinationdrugs.dir/*_signature_foreachdrug_forplots3.json',
        'rawdata/combinationdrugs.dir/')


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

    ## not from original pipeline--put in after for analysis

        dict_3 = {}
        dict_3["up"] = current_up
        dict_3["down"] = current_down
        outfile = '{outfile_root}{drug}_signature_foreachdrugs_forplots3.json'.format(**locals())

        print(outfile)

        # open file
        f = open(outfile, 'w')

        # write to file
        f.write(json.dumps(dict_3, ensure_ascii=False, indent=4))

        # close file
        f.close() 
        

        
#         for sample, rowdata in df_genesig.iterrows():
#             sample_bottom_temp = rowdata["bottom"]
#             sample_top_temp = rowdata["top"]

#             # subtracting drug genes from sample genes
#             for gene in sample_bottom_temp:
#                 if gene in current_up:
#                     sample_bottom_temp.remove(gene)
#             for gene2 in sample_top_temp:
#                 if gene2 in current_down:
#                     sample_top_temp.remove(gene2)
            
#             updated_genesets = {
#                 "top":sample_top_temp,
#                 "bottom":sample_bottom_temp
#             }
            
#             # Extract the top 500 genes from the index
#             dict_updated_genesig[sample] = updated_genesets
        
#         outfile = '{outfile_root}{drug}_updatedsignatures.json'.format(**locals())

#         # open file
#         f = open(outfile, 'w')

#         # write to file
#         f.write(json.dumps(dict_updated_genesig, ensure_ascii=False, indent=4))

#         # close file
#         f.close()

    



# ########################################################################
# ########################################################################
# ########## S10. Finding Combination Drugs from Updated Gene Signatures
# #######################################################################
# #######################################################################

# # ##### Takes updated gene signatures and predicts drug pairs for each of these drugs for each patient based on the patient's gene signature.
# # ### Inputs: Updated list of gene signatures (for each drug)
# # ### Output: List of predicted drugs to be paired with each currently used drug for each patient.

# # ##################################################
# # ########## 1. Query L1000 to find reversing drugs
# # #################################################



# @transform(pairdrugs_updatedsignatures,
#             suffix("_updatedsignatures.json"),
#             "_resultids.json")


# def pairdrugs_queryL1000(infile, outfile):

#     print(outfile)

#     # query the L1000 signatures using gene sets

#     L1000FWD_URL = 'http://amp.pharm.mssm.edu/L1000FWD/'


#     # create empty dictionary for L1000FWD result ids
#     l1000_pair_reversing_resultids = {}

#     with open(infile) as infile2:
#         updatedsig = json.load(infile2)

#     # Loop through dictionary
#     for sample, genesets in updatedsig.items():

#         payload = {
#             # put downregulated genes in the upregulated spot to find drugs that reverse the signature
#             # same with upregulated gewnes in the downregulated spot
#             'up_genes': genesets['bottom'],
#             'down_genes': genesets['top']
#         }
        
#         response = requests.post(L1000FWD_URL + 'sig_search', json=payload)
            
#         if response.status_code == 200:
#             l1000_pair_reversing_resultids[sample] = response.json() #saving results to a dictionary
#             json.dump(response.json(), open('api3_result.json', 'w'), indent=4)      

#     # open file
#     f = open(outfile, 'w')

#     # write to file
#     f.write(json.dumps(l1000_pair_reversing_resultids, ensure_ascii=False, indent=4))

#     # close file
#     f.close()



########################################################################
########## 2. Top 50 signature ids for reversing drugs for each sample
#######################################################################


# @transform(pairdrugs_queryL1000,
#             suffix("_resultids.json"),
#             "_top50results.json")


# def pairdrugs_top50(infile, outfile):

#     print(outfile)

#     with open(infile) as infile2:
#         drugresults = json.load(infile2)
    
    
#     # create empty dictionary
#     l1000_pair_reversing_top50 = {}
    
    
#     L1000FWD_URL = 'http://amp.pharm.mssm.edu/L1000FWD/'
    
    
    
#     for k in drugresults:
#         result_id = str(drugresults[k])
#         result_id_edited = result_id[15:39] # getting just result id from string
        
#         response = requests.get(L1000FWD_URL + 'result/topn/' + result_id_edited)
        
#         if response.status_code == 200:
#             l1000_pair_reversing_top50[k] = response.json() # adding top 50 results to dictionary
#             json.dump(response.json(), open('api4_result.json', 'w'), indent=4)


#     # open file
#     f = open(outfile, 'w')

#     # write to file
#     f.write(json.dumps(l1000_pair_reversing_top50, ensure_ascii=False, indent=4))

#     # close file
#     f.close()




# ######################################################################
# ########## 3. Making top 50 results nested dictionary into dataframe
# ######################################################################


# @transform(pairdrugs_top50,
#             suffix('_top50results.json'),
#             "_top50results_todf.txt")


# def pairdrugs_top50_todf(infile, outfile):
    
#     print(outfile)

#     with open(infile) as infile2:
#         top50_dict = json.load(infile2)


#     temp_df_signatureids = pd.DataFrame(top50_dict).T

#     df_list = []

#     for sample, rowdata in temp_df_signatureids.iterrows(): # row id is sample, other data is rowdata
#         for direction in ["similar"]:
#             df = pd.DataFrame(rowdata[direction])
#             df["sample"] = sample
#             df_list.append(df)
#     df_signatureids = pd.concat(df_list)

#     df_signatureids2 = df_signatureids.set_index("sig_id")
    
#     df_signatureids2.to_csv(outfile, sep="\t")




# ##############################################################
# ########## 4. Getting single singature by id
# ##############################################################

# @transform(pairdrugs_top50_todf,
#             suffix('_top50results_todf.txt'),
#             "_predicteddrugsignature.json")


# def pairdrugs_predicteddrugsignature(infile, outfile):
    
#     print(outfile)

#     L1000FWD_URL = 'http://amp.pharm.mssm.edu/L1000FWD/'

#     df_top50 = pd.read_table(infile)
    
#     list_sigids2 = list(df_top50["sig_id"])

#     dict_pair_predicteddrugsignatures = {}

#     for signature_id in list_sigids2:
#         response = requests.get(L1000FWD_URL + 'sig/' + signature_id)
#         if response.status_code == 200:
#             dict_pair_predicteddrugsignatures[signature_id] = (response.json())
# #       json.dump(response.json(), open('api2_result.json', 'w'), indent=4)


#     # open file
#     f = open(outfile, 'w')

#     # write to file
#     f.write(json.dumps(dict_pair_predicteddrugsignatures, ensure_ascii=False, indent=4))

#     # close file
#     f.close()





# ########################################################################
# ########################################################################
# ########## S11. Combination Drugs Analysis
# #######################################################################
# #######################################################################

# # ##### Takes dictionary of pair drug signatures for each current drug and ???
# # ### Inputs: Updated list of gene signatures (for each drug)
# # ### Output: List of predicted drugs to be paired with each currently used drug for each patient.

# # ##################################################################################################
# # ########## 1. Make dictionary into dataframe with fields of interest and pvals from top50 results
# # ##################################################################################################

# @transform(pairdrugs_predicteddrugsignature,
#             regex(r'(.*)_predicteddrugsignature.json'),
#             add_inputs(r'\1_top50results_todf.txt'),
#             r"\1_df_foranalysis.txt")


# def pairdrugs_df_analysis(infiles, outfile):

#     signaturefile, top50file = infiles

#     print(outfile)

#     with open(signaturefile) as infile2:
#         pair_drugsignatures_dict = json.load(infile2)

#     df_pair_drugsignatures = pd.DataFrame(pair_drugsignatures_dict).T.set_index("sig_id")

#     df_pair_signatureids = pd.read_table(top50file).set_index("sig_id")

#     df_pair_foranalysis = pd.DataFrame()

#     list_sigids = df_pair_drugsignatures.index.tolist()
#     list_pertdesc = (df_pair_drugsignatures["pert_desc"]).tolist()
#     list_pertdose = (df_pair_drugsignatures["pert_dose"]).tolist()
#     list_pertid = (df_pair_drugsignatures["pert_id"]).tolist()
#     list_perttime = (df_pair_drugsignatures["pert_time"]).tolist()
#     df_pair_foranalysis["sig_id"] = list_sigids
#     df_pair_foranalysis["pert_desc"] = list_pertdesc
#     df_pair_foranalysis["pert_dose"] = list_pertdose
#     df_pair_foranalysis["pert_id"] = list_pertid
#     df_pair_foranalysis["pert_time"] = list_perttime
#     df_pair_foranalysis.set_index("sig_id")

#     df_pair_foranalysis2 = df_pair_foranalysis.merge(df_pair_signatureids, how="inner", on="sig_id")

#     df_pair_foranalysis2 = df_pair_foranalysis2.set_index("sig_id")
    
    
#     df_pair_foranalysis2.to_csv(outfile, sep="\t")






# # ##################################################################################################
# # ########## 2. Correcting pertnames using metadata file
# # ##################################################################################################   


# @transform(pairdrugs_df_analysis,
#             suffix('_df_foranalysis.txt'),
#             "_df_correctpertnames_test.txt")


# def pairdrugs_correct_pertnames(infile, outfile):

#     print(outfile)

#     df_pair = pd.read_table(infile).set_index("pert_id")


#     ## to be able to use novel drug data and not have same drug listed under multiple brd ids

#     df_drugsmetadata = pd.read_csv("/Users/maayanlab/Documents/Ovarian Cancer Project/ovarian-cancer/rawdata/Drugs_metadata.csv").set_index("pert_iname")


#     # finding pert_name from pertids for predicted drugs using metadata file --will give name or pertid if no name available


#     list_pair_pertids = list(df_pair.index)

#     list_pair_pertnames = []


#     for name, rowdata in df_drugsmetadata.iterrows():
#         temp_pertid = rowdata["pert_id"]
#         if temp_pertid in list_pair_pertids:
#             list_pair_pertnames.append(name)

#     df_test3 = pd.DataFrame()

#     df_test3["pert_name"] = list_pair_pertnames


#     # save to file
#     df_test3.to_csv(outfile, sep="\t", index = False)

#     # dict_pair_pertnames_counts = {}

#     # for item in list_pair_pertnames:
#     #     count = 0
#     #     for item2 in list_pair_pertnames:
#     #         if item == item2:
#     #             count += 1
#     #         dict_pair_pertnames_counts[item] = count



    

#     # # creating new dataframe of sorted counts of correct pertnames--can be used to make plots
#     # pair_correct_counts_keys = list(dict_pair_pertnames_counts.keys())
#     # pair_correct_counts_values = list(dict_pair_pertnames_counts.values())
#     # df_pair_correctpertnames = pd.DataFrame()
#     # df_pair_correctpertnames["pert_names"] = pair_correct_counts_keys
#     # df_pair_correctpertnames["counts"] = pair_correct_counts_values
#     # df_pair_correctpertnames = df_pair_correctpertnames.sort_values(by="counts", axis=0, ascending=False)



#     # # save to file
#     # df_pair_correctpertnames.to_csv(outfile, sep="\t", index = False)

#     # # open file
#     # f = open(outfile, 'w')

#     # # write to file
#     # f.write(json.dumps(dict_pair_pertnames_counts, ensure_ascii=False, indent=4))

#     # # close file
#     # f.close()






# # ##################################################################################################
# # ########## 3. Getting FDA approval info and MOA using metadata file made by zichen
# # ##################################################################################################   

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





# ##################################################################################################
# ########## 4. Counts for patients with worst prognosis
# ##################################################################################################   

@transform(pairdrugs_df_analysis,
            suffix('_df_foranalysis.txt'),
            "_df_worst_prognosis_counts.txt")


def pairdrugs_worstprognosis(infile, outfile):

    print(outfile)

    ## repeating step from before of correcting pertnames
    df_pair = pd.read_table(infile).set_index("pert_id")

    list_pair_samples = df_pair["samples"].tolist()
    list_pair_samples_edited = []
    for item in list_pair_samples:
        temp_item = item[0:12]
        list_pair_samples_edited.append(temp_item)
    df_pair["samples2"] = list_pair_samples_edited

    df_clinical_prognosis.read_table("/Users/maayanlab/Documents/Ovarian Cancer Project/ovarian-cancer/rawdata/clinical_prognosis_worst.txt")

    df_pair_merged = df_pair.merge(df_clinical_prognosis, left_on="samples2", right_on="submitter_id")

    


#     ## to be able to use novel drug data and not have same drug listed under multiple brd ids

#     df_drugsmetadata = pd.read_csv("/Users/maayanlab/Documents/Ovarian Cancer Project/ovarian-cancer/rawdata/Drugs_metadata.csv").set_index("pert_iname")


#     # finding pert_name from pertids for predicted drugs using metadata file --will give name or pertid if no name available


#     list_pair_pertids = list(df_pair.index)

#     list_pair_pertnames = []


#     for name, rowdata in df_drugsmetadata.iterrows():
#         temp_pertid = rowdata["pert_id"]
#         if temp_pertid in list_pair_pertids:
#             list_pair_pertnames.append(name)

#     df_test3 = pd.DataFrame()

#     df_test3["pert_name"] = list_pair_pertnames


#     # save to file
#     df_test3.to_csv(outfile, sep="\t", index = False)

#     # dict_pair_pertnames_counts = {}

#     # for item in list_pair_pertnames:
#     #     count = 0
#     #     for item2 in list_pair_pertnames:
#     #         if item == item2:
#     #             count += 1
#     #         dict_pair_pertnames_counts[item] = count



    

#     # # creating new dataframe of sorted counts of correct pertnames--can be used to make plots
#     # pair_correct_counts_keys = list(dict_pair_pertnames_counts.keys())
#     # pair_correct_counts_values = list(dict_pair_pertnames_counts.values())
#     # df_pair_correctpertnames = pd.DataFrame()
#     # df_pair_correctpertnames["pert_names"] = pair_correct_counts_keys
#     # df_pair_correctpertnames["counts"] = pair_correct_counts_values
#     # df_pair_correctpertnames = df_pair_correctpertnames.sort_values(by="counts", axis=0, ascending=False)



#     # # save to file
#     # df_pair_correctpertnames.to_csv(outfile, sep="\t", index = False)

#     # # open file
#     # f = open(outfile, 'w')

#     # # write to file
#     # f.write(json.dumps(dict_pair_pertnames_counts, ensure_ascii=False, indent=4))

#     # # close file
#     # f.close()










##################################################
##################################################
########## Run pipeline
##################################################
##################################################
pipeline_run([sys.argv[-1]], multiprocess=1, verbose=1)
print('Done!')