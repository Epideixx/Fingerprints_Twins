# Fingerprints Twins

Study of the heritability of neurophysiological fingerprints by using MEG records from Twins.
The used Dataset comes from the HCP study, containing 89 MEG recordings on MZ Twins, DZ Twins and others.

This project is part of a study on heritability at the NeuroSPEED, McGill.

(Description and result to come ...)

## Steps :

1. Get the datasets as csv (row : subject_ID, column : ROI_frequency)
If the signals are given for each subject and each recording, the file conversion_Matlab_to_csv.py will convert to the desired format

2. PSD_accuracy.py will compute the accuracy of the identification of the subject themself as well as the identification of the twin.
We use a bootstrap method, computing the accuracy for 1000 bootstraps of 90 % of the original dataset

3. PSD_correlation.py will compute the correlation between the ...

4. PSD_ICC_Heritability.py ...