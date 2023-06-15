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
DATA_PATH = "new_Data"
FOLDER_RESULTS = "Results_Log_Schaefer_test"
ONLY_GT = True #Do not change the result here as we do not consider the relation between different participants
WITHOUT_TWINS = True
N_RESAMPLE = 10


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

# --- PATH TO SAVE THE RESULTS ---

if not os.path.exists(FOLDER_RESULTS):
    os.mkdir(FOLDER_RESULTS)

if ONLY_GT:
    FOLDER_RESULTS = os.path.join(FOLDER_RESULTS, "Only_GT")
else :
    FOLDER_RESULTS = os.path.join(FOLDER_RESULTS, "Include_SR")
if not os.path.exists(FOLDER_RESULTS):
    os.mkdir(FOLDER_RESULTS)

FOLDER_RESULTS = os.path.join(FOLDER_RESULTS, "PSD_ICC_Fingerprint")
if not os.path.exists(FOLDER_RESULTS):
    os.mkdir(FOLDER_RESULTS)

# ---   IMPORT DATA   ---

record_1 = pd.read_csv(os.path.join(DATA_PATH, "record_1_Schaefer.csv"), index_col="Subject_ID")
record_2 = pd.read_csv(os.path.join(DATA_PATH, "record_2_Schaefer.csv"), index_col="Subject_ID")
record_3 = pd.read_csv(os.path.join(DATA_PATH, "record_3_Schaefer.csv"), index_col="Subject_ID")

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
if ONLY_GT :
    for i in all_data_restricted.index:
        subj_ID = all_data_restricted.loc[i, "Subject"]
        fam_ID = all_data_restricted.loc[i, "Family_ID"]

        if fam_ID not in twins_dict:
            if len(all_data_restricted.loc[i, "ZygosityGT"]) >= 2:
                twins_dict[fam_ID] = {"type" : all_data_restricted.loc[i, "ZygosityGT"], "subjects" : [subj_ID] }
            else : 
                twins_dict[fam_ID] = {"type" : "NT", "subjects" : [subj_ID] } #NT = NoTwin
        else : 
            # We have only one pair of "normal" siblings and still not sure what they are
            # So, if classified as Not twins, it means the relation is not certain
            twins_dict[fam_ID]["subjects"].append(subj_ID)
else :
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
nb_not_twins = len(ids_NT)
nb_pair_not_twins = len([pair["subjects"] for pair in list(twins_dict.values()) if len(pair["subjects"]) >= 2 and pair["type"] == "NT"])

print("The total number of participants is : ", len(all_data_restricted))
print("We will work with {} pairs of siblings : {} pairs of MZ twins, {} pairs of DZ twins and {} pairs of siblings with uncertain relation".format(nb_pair, nb_pair_MZ, nb_pair_DZ, nb_pair_not_twins))

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

# ---   MAIN   ---

if WITHOUT_TWINS :

    iccs = []

    # Bootstrap
    for k in tqdm(range(N_RESAMPLE), total=N_RESAMPLE): 

        zscore_1 = zscore(record_1, axis = 1).rename(columns = {col : col + "_A" for col in record_1.columns})
        zscore_2 = zscore(record_2, axis = 1).rename(columns = {col : col + "_B" for col in record_1.columns})
        zscore_3 = zscore(record_3, axis = 1).rename(columns = {col : col + "_C" for col in record_1.columns})

        ids_to_keep = [list(twins_dict.values())[i]['subjects'][np.random.randint(len(list(twins_dict.values())[i]['subjects']))] for i in range(len(twins_dict.values()))]
        mask = [True if id in ids_to_keep else False for id in zscore_1.index]

        zscore_1 = zscore_1[mask]
        zscore_2 = zscore_2[mask]
        zscore_3 = zscore_3[mask]

        print(zscore_1.shape[0])

        n_subs = zscore_1.shape[0]
        n_measurements = 3

        icc = np.zeros(n_ROI * n_freqs)
        for i, feature in enumerate(record_1.columns):
            df = pd.concat([zscore_1[feature + "_A"], zscore_2[feature + "_B"], zscore_3[feature + "_C"]], axis = 1)
            icc[i] = compute_icc(df, n = n_subs, k = n_measurements)

        icc = np.reshape(icc, (n_ROI, n_freqs))
        icc = pd.DataFrame(icc, columns=freqs, index=list_ROI)

        iccs.append(icc)

    # Average all of the bootstraps
    iccs = np.array(iccs)
    icc = np.mean(iccs, axis = 0)

    icc = pd.DataFrame(icc, columns=freqs, index=list_ROI)
    save_file = os.path.join(FOLDER_RESULTS, "ICC_without_twins.csv")
    icc.to_csv(save_file, index_label="ROI")


    icc_avg_per_band = icc.copy()
    for ind, band in enumerate(bands):
        cols = [c for c in icc.columns if float(c) >= band[0] and float(c) < band[1]]
        icc_avg_per_band[bands_names[ind]] = np.mean(icc[cols], axis = 1)

    icc_avg_per_band = icc_avg_per_band[bands_names]
    icc_avg_per_band["total_avg"] = np.mean(icc_avg_per_band[bands_names], axis = 1)

    save_file = os.path.join(FOLDER_RESULTS, "ICC_without_twins_avg_per_band.csv")
    icc_avg_per_band.to_csv(save_file, index_label="ROI")

else :
    save_file = os.path.join(FOLDER_RESULTS, "ICC.csv")

    zscore_1 = zscore(record_1, axis = 1).rename(columns = {col : col + "_A" for col in record_1.columns})
    zscore_2 = zscore(record_2, axis = 1).rename(columns = {col : col + "_B" for col in record_1.columns})
    zscore_3 = zscore(record_3, axis = 1).rename(columns = {col : col + "_C" for col in record_1.columns})

    n_subs = len(record_1.index)
    n_measurements = 3

    icc = np.zeros(n_ROI * n_freqs)
    for i, feature in tqdm(enumerate(record_1.columns), total = len(record_1.columns)):
        df = pd.concat([zscore_1[feature + "_A"], zscore_2[feature + "_B"], zscore_3[feature + "_C"]], axis = 1)
        icc[i] = compute_icc(df, n = n_subs, k = n_measurements)

    icc = np.reshape(icc, (n_ROI, n_freqs))
    icc = pd.DataFrame(icc, columns=freqs, index=list_ROI)
    save_file = os.path.join(FOLDER_RESULTS, "ICC.csv")
    icc.to_csv(save_file, index_label="ROI")

    icc_avg_per_band = icc.copy()
    for ind, band in enumerate(bands):
        cols = [c for c in icc.columns if float(c) >= band[0] and float(c) < band[1]]
        icc_avg_per_band[bands_names[ind]] = np.mean(icc[cols], axis = 1)

    icc_avg_per_band = icc_avg_per_band[bands_names]
    icc_avg_per_band["total_avg"] = np.mean(icc_avg_per_band[bands_names], axis = 1)

    save_file = os.path.join(FOLDER_RESULTS, "ICC_avg_per_freq_band.csv")
    icc_avg_per_band.to_csv(save_file, index_label="ROI")

