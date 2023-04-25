# ---   PARAMETERS   ---
# Frequency bands
BROADBAND = (0.0, 150.0)
DELTA = (0.0, 4.0)
THETA = (4.0, 8.0)
ALPHA = (8.0, 13.0)
BETA = (13.0, 30.0)
GAMMA = (30.0, 50.0)
HIGH_GAMMA = (50.0, 150.0)

# Parameters 
DATA_PATH = "New_Data"
FOLDER_RESULTS = "Results_Log_DK/ICC_Only"
N_RESAMPLE = 1000

NB_ROI = 200
NB_FREQ = 301

# ---   DEPENDENCIES   ---

import os
import re
import numpy as np
import pandas as pd
from scipy.stats import pearsonr, zscore
from tqdm import tqdm

# To plot
import seaborn as sns
sns.set_theme(style="white")
from matplotlib import pyplot as plt

# ---   IMPORT DATA   ---

record_1 = pd.read_csv(os.path.join(DATA_PATH, "record_1_DK.csv"), index_col="Subject_ID")
record_2 = pd.read_csv(os.path.join(DATA_PATH, "record_2_DK.csv"), index_col="Subject_ID")
record_3 = pd.read_csv(os.path.join(DATA_PATH, "record_3_DK.csv"), index_col="Subject_ID")

# Log the data

record_1 = np.log(record_1)
record_2 = np.log(record_2)
record_3 = np.log(record_3)

# ---   USEFULL ITEMS   ---

### Annotate Data such that twins are linked ###

# Import extra data (confidential)
all_data_restricted_filename = "new_Data/All_Data_RESTRICTED.csv"
all_data_restricted = pd.read_csv(all_data_restricted_filename)

# Remove subjects who don't have a MEG recording
subjects_with_MEG = record_1.index
mask = [True if subj_id in subjects_with_MEG else False for subj_id in all_data_restricted["Subject"]]
all_data_restricted = all_data_restricted[mask]
print("Number of subjects with MEG :", len(all_data_restricted))

# Stock individuals in a dictionnary identified by the Family ID, and giving the pairs of twins as well as the type of twins
twins_dict = {} # Items in twins_dict => familly ID : {type_twins (MZ/DZ/NT), list_ID_subject_same_family}
for i in all_data_restricted.index:
    subj_ID = all_data_restricted.loc[i, "Subject"]
    fam_ID = all_data_restricted.loc[i, "Family_ID"]

    if fam_ID not in twins_dict:
        if len(all_data_restricted.loc[i, "ZygosityGT"]) >= 2:
            twins_dict[fam_ID] = {"type" : all_data_restricted.loc[i, "ZygosityGT"], "subjects" : [subj_ID] }
        elif all_data_restricted.loc[i, "ZygositySR"] in ["MZ", "DZ"]:
            twins_dict[fam_ID] = {"type" : all_data_restricted.loc[i, "ZygositySR"], "subjects" : [subj_ID] }
        elif all_data_restricted.loc[i, "ZygositySR"] == "NotMZ":
            twins_dict[fam_ID] = {"type" : "DZ", "subjects" : [subj_ID] }
        else : 
            twins_dict[fam_ID] = {"type" : "NT", "subjects" : [subj_ID] } #NT = NoTwin
    else : 
        twins_dict[fam_ID]["subjects"].append(subj_ID)

# Go from Subject ID to its row number in record and vice-versa 
# (usefull when using the numpy conversion)
subjects_row_to_id = {i : id for i, id in enumerate(record_1.index)}
subject_id_to_row = {id : i for i, id in enumerate(record_1.index)}

# Create dictionnary to rename individuals to :
# - Twin MZ 1A and Twin MZ 1B
# - Twin DZ 1A and Twin DZ 1B
# - NotTwin 1

count_MZ = 1
count_DZ = 1
count_NT = 1

rename_twins = {} # {subject_ID : new_name}
for twins in twins_dict.values():
    if twins["type"] == "MZ" and len(twins["subjects"]) >= 2:
        rename_twins[twins["subjects"][0]] = "Twin_MZ_" + str(count_MZ) + "A"
        rename_twins[twins["subjects"][1]] = "Twin_MZ_" + str(count_MZ) + "B"
        count_MZ += 1
    elif twins["type"] == "DZ" and len(twins["subjects"]) >= 2:
        rename_twins[twins["subjects"][0]] = "Twin_DZ_" + str(count_DZ) + "A"
        rename_twins[twins["subjects"][1]] = "Twin_DZ_" + str(count_DZ) + "B"
        count_DZ += 1
    else : 
        rename_twins[twins["subjects"][0]] = "NotTwin" + str(count_NT)
        count_NT += 1
        if len(twins["subjects"]) >= 2: 
            rename_twins[twins["subjects"][1]] = "NotTwin" + str(count_NT)
            count_NT += 1


# Lists of IDs depending on the category of the subject
ids_MZ = [k for k, v in rename_twins.items() if "Twin_MZ" in v]
ids_DZ = [k for k, v in rename_twins.items() if "Twin_DZ" in v]
ids_NT = [k for k, v in rename_twins.items() if "NotTwin" in v]

def f7(seq):
    """ Remove duplicates of a list and preserve the order. """
    seen = set()
    seen_add = seen.add
    return [x for x in seq if not (x in seen or seen_add(x))]

list_ROI = f7([c.replace(re.search("_[0-9]+\.[0-9]*", c).group(0), "") for c in record_1.columns])
n_ROI = len(list_ROI)
freqs = sorted(set([float(re.search("[0-9]+\.[0-9]*", c).group(0)) for c in record_1.columns]))
n_freqs = len(freqs)

bands = [DELTA, THETA, ALPHA, BETA, GAMMA, HIGH_GAMMA]
bands_names = ["DELTA", "THETA", "ALPHA", "BETA", "GAMMA", "HIGH GAMMA"]

# Summary 

nb_pair = len([pair["subjects"] for pair in list(twins_dict.values()) if len(pair["subjects"]) >= 2])
nb_pair_MZ = len(ids_MZ)//2
nb_pair_DZ = len(ids_DZ)//2
nb_pair_not_twins = len([pair["subjects"] for pair in list(twins_dict.values()) if len(pair["subjects"]) >= 2 and pair["type"] == "NT"])

print("The total number of participants is : ", len(all_data_restricted))
print("We will work with {} pairs of siblings : {} pairs of MZ twins, {} pairs of DZ twins and {} pairs of normal siblings".format(nb_pair, nb_pair_MZ, nb_pair_DZ, nb_pair_not_twins))

# ---   FUNCTIONS   ---

def compute_icc(df, n = 89, k = 2):
    df_b = n-1
    df_w = n*(k-1)

    x = df.to_numpy()
    x_w_mean = x.mean(axis = 1)
    x_g_mean = x.mean()
    ss_t = ((x - x_g_mean) ** 2).sum()
    ss_w = ((x - np.expand_dims(x_w_mean, axis = 1)) ** 2).sum()
    ss_b = ss_t - ss_w
    ms_b = ss_b / df_b
    ms_w = ss_w / df_w
    res = (ms_b - ms_w) / (ms_b + ((k-1)*ms_w))

    return res

# ---   MAIN ICC   ---

save_file = os.path.join(FOLDER_RESULTS, "ICC_only_MZ.csv")

if not(os.path.exists(save_file)) : 

    zscore_1 = zscore(record_1, axis = 1).rename(columns = {col : col + "_A" for col in record_1.columns})
    zscore_2 = zscore(record_2, axis = 1).rename(columns = {col : col + "_B" for col in record_1.columns})
    zscore_3 = zscore(record_3, axis = 1).rename(columns = {col : col + "_C" for col in record_1.columns})

    mask_MZ = [True if subj_id in ids_MZ else False for subj_id in record_1.index]
    zscore_1_MZ = zscore_1[mask_MZ].reindex(ids_MZ)
    zscore_2_MZ = zscore_2[mask_MZ].reindex(ids_MZ)
    zscore_3_MZ = zscore_3[mask_MZ].reindex(ids_MZ)


    n_subs = np.sum(mask_MZ)
    n_measurements = 3

    icc = np.zeros(n_ROI * n_freqs)
    for i, feature in tqdm(enumerate(record_1.columns), total = len(record_1.columns)):
        df = pd.concat([zscore_1_MZ[feature + "_A"], zscore_2_MZ[feature + "_B"], zscore_3_MZ[feature + "_C"]], axis = 1)
        icc[i] = compute_icc(df, n = n_subs, k = n_measurements)

    icc = np.reshape(icc, (n_ROI, n_freqs))
    icc = pd.DataFrame(icc, columns=freqs, index=list_ROI)
    save_file = os.path.join(FOLDER_RESULTS, "ICC_only_MZ.csv")
    icc.to_csv(save_file, index_label="ROI")

    icc_avg_per_band = icc.copy()
    for ind, band in enumerate(bands):
        cols = [c for c in icc.columns if float(c) >= band[0] and float(c) < band[1]]
        icc_avg_per_band[bands_names[ind]] = np.mean(icc[cols], axis = 1)

    icc_avg_per_band = icc_avg_per_band[bands_names]
    icc_avg_per_band["total_avg"] = np.mean(icc_avg_per_band[bands_names], axis = 1)

    save_file = os.path.join(FOLDER_RESULTS, "ICC_only_MZ_avg.csv")
    icc_avg_per_band.to_csv(save_file, index_label="ROI")



save_file = os.path.join(FOLDER_RESULTS, "ICC_only_DZ.csv")

if not(os.path.exists(save_file)) : 

    zscore_1 = zscore(record_1, axis = 1).rename(columns = {col : col + "_A" for col in record_1.columns})
    zscore_2 = zscore(record_2, axis = 1).rename(columns = {col : col + "_B" for col in record_1.columns})
    zscore_3 = zscore(record_3, axis = 1).rename(columns = {col : col + "_C" for col in record_1.columns})

    mask_DZ = [True if subj_id in ids_DZ else False for subj_id in record_1.index]
    zscore_1_DZ = zscore_1[mask_DZ].reindex(ids_DZ)
    zscore_2_DZ = zscore_2[mask_DZ].reindex(ids_DZ)
    zscore_3_DZ = zscore_3[mask_DZ].reindex(ids_DZ)


    n_subs = np.sum(mask_DZ)
    n_measurements = 3

    icc = np.zeros(n_ROI * n_freqs)
    for i, feature in tqdm(enumerate(record_1.columns), total = len(record_1.columns)):
        df = pd.concat([zscore_1_DZ[feature + "_A"], zscore_2_DZ[feature + "_B"], zscore_3_DZ[feature + "_C"]], axis = 1)
        icc[i] = compute_icc(df, n = n_subs, k = n_measurements)

    icc = np.reshape(icc, (n_ROI, n_freqs))
    icc = pd.DataFrame(icc, columns=freqs, index=list_ROI)
    save_file = os.path.join(FOLDER_RESULTS, "ICC_only_DZ.csv")
    icc.to_csv(save_file, index_label="ROI")

    icc_avg_per_band = icc.copy()
    for ind, band in enumerate(bands):
        cols = [c for c in icc.columns if float(c) >= band[0] and float(c) < band[1]]
        icc_avg_per_band[bands_names[ind]] = np.mean(icc[cols], axis = 1)

    icc_avg_per_band = icc_avg_per_band[bands_names]
    icc_avg_per_band["total_avg"] = np.mean(icc_avg_per_band[bands_names], axis = 1)

    save_file = os.path.join(FOLDER_RESULTS, "ICC_only_DZ_avg.csv")
    icc_avg_per_band.to_csv(save_file, index_label="ROI")