# Individual differences of neurophysiological brain activity are heritable  

## Introduction

Brain-fingerprint is a novel technique used to explore the inter-individual differrences in the brain. We are capable of accurately differentiating individuals in a cohort based only on the resting-stage spectral brain activy using magnetoencephalography (MEG).
However, the influence of genes on those differences remain unknown.
We thus explored how heritable were the fingerprints, and what genes relate to it.

## Description of the project

Here are the different steps we did :

- We apply the brain-fingerprint method to diffenrentiate both individuals and twin pairs : PSD_accuracy.py
- To compare with chance-level, we followed the same pipeline, but with the empty room recording : PSD_accuracy_empty_room.py
- We futhermore explored the same method, but restricting ourselves to shorter segments of recording : PSD_accuracy_30_sec.py
- We observed the correlations between individuals, and the distributions of correlation within twin pairs : PSD_correlations.py
- To evaluate how useful are the different features from the spectral brain for fingerprinting, we computed the intraclass correlation (ICC) : PSD_ICC_Fingerprint.py
- We computed then the heritability of the same features, using the Falconer's formula : PSD_Heritability.py
- We also computed the the heritability of the anatomical features, to see the influence of anatomy on the spectrel-brain heritability : PSD_Heritability_Anatomy.py

All of those steps are merged together in main.py

Finally, statistical analyses and plots were performed : plots_Schaefer.R


## Results

See Preprint : LINK


## How to use this repository

In addition to the provided data on this github repository, you will require :
- to request the access to the resticted data for the dataset for the MEG subjects (https://db.humanconnectome.org/), and add the non-resticted data to the folder "Data" as "All_Data.csv", as well as the restricted dataset containing the information for the twin pairs, to name "All_Data_RESTRICTED.csv".
- to download the "Schaefer_30_second" folder if you want to run the analysis using the short segments and to upload it in "Data
- to download the "Schaefer_artifacts_correction" folder if you want to run the analysis using the PSDs after regressing out the artifacts



All the results to conduct the ananlyses can be computed  using **main.py** :

Example :
'''
python3 main.py --acc True -- only_gt True 
'''

Here are the different parameters :
- --acc : either if we compute the accuracy or not
- --acc_30_sec : either if we compute the accuracy or not for the 30-sec PSDs
- --acc_empty_room : either if we compute the accuracy of the empty room or not
- --corr : either if we compute the correlations or not
- --heritability : either if we compute the heritability or not
- --heritability_anatomy : either if we compute the hertiability of anatomy or not
- --fingerprint : either if we compute the ICC for fingerprinting or not
- --data : path of the folder containing the data
- --data_30_sec :  path of the folder containing the 30-seconds PSDs
- --results : path for the results
- --only_gt : choice of if we use only twin pairs based on genetic test ("True"), or if we also use the self reported twins ("False")
- --correlation_type : correlation used for heritability (pearson or icc)
- --n_resample : number of bootstraps
- --icc_without_twins : either ICC for fingerprinting is computed with or without twins