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
import sys
import pandas as pd
import numpy as np
import scipy.stats as ss

##### 2. Custom modules #####
# Pipeline running

#############################################
########## 2. General Setup
#############################################
##### 1. Variables #####

# TCGA expression normalized (fpkm) data file
expfpkmfile = /Users/maayanlab/Documents/Ovarian Cancer Project/ovarian-cancer/rawdata/TCGA-OV-fpkm-uq.txt

#######################################################
#######################################################
########## S1. Download Data and Make Dataframe
#######################################################
#######################################################

##### TCGA-fpkm expression data that is downloaded is made into a dataframe.
### Input: TCGA-fpkm expression data (.txt file)
### Output: Dataframe of the TCGA expression data--the index is the gene symbol and the column headers are sample ids.

#############################################
########## 1. Creating dataframe
#############################################

@transform(input = expfpkmfile, output = /rawdata/df_expfpkm.txt)
df_expfpkm = pd.read_table(expfpkmfile).set_index('gene_symbol')
df_expfpkm.to_csv("/rawdata/df_expfpkm.txt", sep=" ")




#######################################################
#######################################################
########## S2. Find Gene Signatures for Each Sample
#######################################################
#######################################################

##### TCGA-fpkm expression dataframe is used to find the gene signature for each sample.
### Input: Dataframe of the TCGA expression data--the index is the gene symbol and the column headers are sample ids.
### Output: Dictionary of the top 500 and bottom 500 genes for each sample

############################################################
########## 1.Z-scores across rows and columns of dataframe
############################################################

@transform(input = /rawdata/df_expfpkm.txt, output = /rawdata_expfpkm_zscore.txt)
# axis=0 is along rows
df_expfpkm_zscore = df_expfpkm.T.apply(ss.zscore, axis=0).T.dropna().apply(ss.zscore, axis=0)

############################################################
########## 2. Creating dictionary of gene signatures
############################################################

# define empty dictionary
results = {}


# Loop through samples
for sample in df_expfpkm_zscore.columns[5:]:
    
    # Extract column
    col = df_expfpkm_zscore[sample].sort_values(ascending=False) # sort values in decreasing order
    
    genesets = {
        "top":col.index[:500].tolist(),
        "bottom":col.index[-500:].tolist()

    }
    # Extract the top 500 genes from the index
    results[sample] = genesets
















