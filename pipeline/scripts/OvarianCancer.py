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

#######################################################
#######################################################
########## S1. Download Data and Make Dataframe
#######################################################
#######################################################
##### TCGA-fpkm expression data that is downloaded is made into a dataframe.
### Input: TCGA-fpkm expression data (.txt file)
### Output: Dataframe of the TCGA expression data--the index is the gene symbol and the column headers are sample ids.
#############################################
########## 1. Creating Dataframe
#############################################

expfpkmfile = /Users/maayanlab/Documents/Ovarian Cancer Project/ovarian-cancer/rawdata/TCGA-OV-fpkm-uq.txt
df_expfpkm = pd.read_table(expfpkmfile).set_index('gene_symbol')


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








