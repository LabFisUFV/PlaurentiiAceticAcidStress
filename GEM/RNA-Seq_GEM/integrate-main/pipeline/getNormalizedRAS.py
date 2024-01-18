#!/usr/bin/env python
# coding: utf-8

import os
import pandas as pd
import genericLib as gL

# setting working dirs
workingDirs = gL.setWorkingDirs()
OUTDIR = workingDirs[2]

# setting input data
inputFileName = 'papla-GEM_RAS'
outputFileName = 'papla-GEM_wNormalizedRAS'

# Load RAS dataset
RAS = pd.read_csv(os.path.join(OUTDIR, inputFileName + '.csv'), sep="\t", index_col = 'Rxn')
## Remaining NaN in the matrix will be substituted with 1 to indicate that a missing RAS value is originally present for these reactions due to emtpy GPR rules
RAS = RAS.fillna(1)

# Normalization for each reaction on the sample having the highest RAS. Here samples from the transcriptomics dataset in rawData directory are used.
# Change sample names and number of replicates (nReplicas variable) of each sample accordingly if the input transcriptomic dataset is different.
nReplicas = 3.0
RAS['mean_Par_wo']=(RAS['Par_wo_1']+RAS['Par_wo_2']+RAS['Par_wo_3'])/nReplicas
RAS['mean_Par_Ac']=(RAS['Par_Ac_1']+RAS['Par_Ac_2']+RAS['Par_Ac_3'])/nReplicas
RAS['mean_ATS_wo']=(RAS['ATS_wo_1']+RAS['ATS_wo_2']+RAS['ATS_wo_3'])/nReplicas
RAS['mean_ATS_Ac']=(RAS['ATS_Ac_1']+RAS['ATS_Ac_2']+RAS['ATS_Ac_3'])/nReplicas

RASNormalized = RAS[['mean_Par_wo','mean_Par_Ac','mean_ATS_wo','mean_ATS_Ac']]
maxValues = RASNormalized.max(axis=1)

RASNormalized['norm_Par_wo']=RASNormalized['mean_Par_wo']/maxValues
RASNormalized['norm_Par_Ac']=RASNormalized['mean_Par_Ac']/maxValues
RASNormalized['norm_ATS_wo']=RASNormalized['mean_ATS_wo']/maxValues
RASNormalized['norm_ATS_Ac']=RASNormalized['mean_ATS_Ac']/maxValues


RASNormalized = RASNormalized.fillna(1)

## Mask rows where all values were originally equal to 0, by putting 0: if all values in rows are equal to 0, then the maximum value in this row is 0 and 0 / 0 is
# undetermineted which means NaN that before was substituted with 1. However, it's important to mantain the 0 value to indicate that the RAS is not missing but it's 0 and
# consequently the reactions is off.
RASNormalized.loc[(RAS.eq(0).all(1))] = 0

# save output matrix
RASNormalized.to_csv(os.path.join(OUTDIR,outputFileName + '.csv'), sep = '\t', index = True)
