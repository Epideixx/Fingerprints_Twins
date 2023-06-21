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
DATA_PATH = "Data/Schaefer"
FOLDER_RESULTS = "Results_Log_Schaefer_test"
ONLY_GT = True 
CORRELATION_TYPE = "icc" # "icc" or "pearson"

def main(data_path = DATA_PATH, main_folder_results = FOLDER_RESULTS, only_gt = ONLY_GT, correlation_type = CORRELATION_TYPE):

    # --- PATH TO SAVE THE RESULTS ---

    if not os.path.exists(main_folder_results):
        os.mkdir(main_folder_results)

    if only_gt:
        main_folder_results = os.path.join(main_folder_results, "Only_GT")
    else :
        main_folder_results = os.path.join(main_folder_results, "Include_SR")
    if not os.path.exists(main_folder_results):
        os.mkdir(main_folder_results)

    main_folder_results = os.path.join(main_folder_results, "PSD_Heritability")
    if not os.path.exists(main_folder_results):
        os.mkdir(main_folder_results)

    main_folder_results = os.path.join(main_folder_results, correlation_type)
    if not os.path.exists(main_folder_results):
        os.mkdir(main_folder_results)


    # ---   IMPORT DATA   ---

    record_1 = pd.read_csv(os.path.join(data_path, "record_1.csv"), index_col="Subject_ID")
    record_2 = pd.read_csv(os.path.join(data_path, "record_2.csv"), index_col="Subject_ID")
    record_3 = pd.read_csv(os.path.join(data_path, "record_3.csv"), index_col="Subject_ID")

    # Log the data

    record_1 = np.log(record_1)
    record_2 = np.log(record_2)
    record_3 = np.log(record_3)

    # ---   USEFULL ITEMS   ---

    ### Annotate Data such that twins are linked ###

    # Import extra data (confidential)
    all_data_restricted_filename = "Data/All_Data_RESTRICTED.csv"
    all_data_restricted = pd.read_csv(all_data_restricted_filename)

    # Remove subjects who don't have a MEG recording
    subjects_with_MEG = record_1.index
    mask = [True if subj_id in subjects_with_MEG else False for subj_id in all_data_restricted["Subject"]]
    all_data_restricted = all_data_restricted[mask]
    print("Number of subjects with MEG :", len(all_data_restricted))

    # Stock individuals in a dictionnary identified by the Family ID, and giving the pairs of twins as well as the type of twins
    twins_dict = {} # Items in twins_dict => familly ID : {type_twins (MZ/DZ/NT), list_ID_subject_same_family}
    if only_gt :
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

    # MZ Correlation

    save_file = os.path.join(main_folder_results, "corr_MZ.csv")

    mask_MZ = [True if subj_id in ids_MZ else False for subj_id in record_1.index]
    record_1_MZ = record_1[mask_MZ].reindex(ids_MZ).rename(index=rename_twins)
    zscore_1_MZ = zscore(record_1_MZ, axis = 1)
    record_2_MZ = record_2[mask_MZ].reindex(ids_MZ).rename(index=rename_twins)
    zscore_2_MZ = zscore(record_2_MZ, axis = 1)
    record_3_MZ = record_3[mask_MZ].reindex(ids_MZ).rename(index=rename_twins)
    zscore_3_MZ = zscore(record_3_MZ, axis = 1)

    corr_MZ = np.zeros(zscore_1_MZ.shape[1])
    corr_MZ_sub = np.zeros((9, zscore_1_MZ.shape[1]))

    n_subs = len(record_1_MZ)*9 # Number of pair of Twins !
    n_measurements = 2

    n = n_subs
    k = n_measurements

    for i, feature in tqdm(enumerate(record_1.columns), total = len(record_1.columns)):
        df1 = pd.DataFrame(np.array([[zscore_1_MZ.iloc[k][feature], zscore_1_MZ.iloc[k+1][feature]] for k in range(0, len(record_1_MZ), 2)]))
        df2 = pd.DataFrame(np.array([[zscore_2_MZ.iloc[k][feature], zscore_2_MZ.iloc[k+1][feature]] for k in range(0, len(record_1_MZ), 2)]))
        df3 = pd.DataFrame(np.array([[zscore_3_MZ.iloc[k][feature], zscore_3_MZ.iloc[k+1][feature]] for k in range(0, len(record_1_MZ), 2)]))
        df4 = pd.DataFrame(np.array([[zscore_1_MZ.iloc[k][feature], zscore_2_MZ.iloc[k+1][feature]] for k in range(0, len(record_1_MZ), 2)]))
        df5 = pd.DataFrame(np.array([[zscore_2_MZ.iloc[k][feature], zscore_1_MZ.iloc[k+1][feature]] for k in range(0, len(record_1_MZ), 2)]))
        df6 = pd.DataFrame(np.array([[zscore_1_MZ.iloc[k][feature], zscore_3_MZ.iloc[k+1][feature]] for k in range(0, len(record_1_MZ), 2)]))
        df7 = pd.DataFrame(np.array([[zscore_3_MZ.iloc[k][feature], zscore_1_MZ.iloc[k+1][feature]] for k in range(0, len(record_1_MZ), 2)]))
        df8 = pd.DataFrame(np.array([[zscore_2_MZ.iloc[k][feature], zscore_3_MZ.iloc[k+1][feature]] for k in range(0, len(record_1_MZ), 2)]))
        df9 = pd.DataFrame(np.array([[zscore_3_MZ.iloc[k][feature], zscore_2_MZ.iloc[k+1][feature]] for k in range(0, len(record_1_MZ), 2)]))

        for j, df_j in enumerate([df1, df2, df3, df4, df5, df6, df7, df8, df9]):
            if correlation_type == "icc":
                corr_MZ_sub[j][i] = compute_icc(df_j, n = len(record_1_MZ), k = n_measurements)
            elif correlation_type == "pearson" :
                corr_MZ_sub[j][i] = pearsonr(df_j[df_j.columns[0]], df_j[df_j.columns[1]])[0]
        df = pd.concat([df1, df2, df3, df4, df5, df6, df7, df8, df9], axis = 0).rename(columns={0:"Twin A", 1: "Twin B"})
        if correlation_type == "icc":
            corr_MZ[i] = compute_icc(df, n = n_subs, k = n_measurements)
        elif correlation_type == "pearson" :
            corr_MZ[i] = pearsonr(df[df.columns[0]], df[df.columns[1]])[0]

    corr_MZ = np.reshape(corr_MZ, (n_ROI, n_freqs))
    corr_MZ = pd.DataFrame(corr_MZ, columns=freqs, index=list_ROI)
    corr_MZ.to_csv(save_file, index_label="ROI")

    for j in range(9):
        save_file = "corr_MZ_partial_" + str(j) + ".csv"
        save_file = os.path.join(main_folder_results, save_file)
        corr_MZ = np.reshape(corr_MZ_sub[j], (n_ROI, n_freqs))
        corr_MZ = pd.DataFrame(corr_MZ, columns=freqs, index=list_ROI)
        corr_MZ.to_csv(save_file, index_label="ROI")

    # DZ Correlation
    save_file = os.path.join(main_folder_results, "corr_DZ.csv")

    mask_DZ = [True if subj_id in ids_DZ else False for subj_id in record_1.index]
    record_1_DZ = record_1[mask_DZ].reindex(ids_DZ).rename(index=rename_twins)
    zscore_1_DZ = zscore(record_1_DZ, axis = 1)
    record_2_DZ = record_2[mask_DZ].reindex(ids_DZ).rename(index=rename_twins)
    zscore_2_DZ = zscore(record_2_DZ, axis = 1)
    record_3_DZ = record_3[mask_DZ].reindex(ids_DZ).rename(index=rename_twins)
    zscore_3_DZ = zscore(record_3_DZ, axis = 1)

    corr_DZ = np.zeros(zscore_1_DZ.shape[1])
    corr_DZ_sub = np.zeros((9, zscore_1_MZ.shape[1]))

    n_subs = len(record_1_DZ)*9 # Number of pair of Twins !
    n_measurements = 2

    n = n_subs
    k = n_measurements

    for i, feature in tqdm(enumerate(record_1.columns), total = len(record_1.columns)):
        df1 = pd.DataFrame(np.array([[zscore_1_DZ.iloc[k][feature], zscore_1_DZ.iloc[k+1][feature]] for k in range(0, len(record_1_DZ), 2)]))
        df2 = pd.DataFrame(np.array([[zscore_2_DZ.iloc[k][feature], zscore_2_DZ.iloc[k+1][feature]] for k in range(0, len(record_1_DZ), 2)]))
        df3 = pd.DataFrame(np.array([[zscore_3_DZ.iloc[k][feature], zscore_3_DZ.iloc[k+1][feature]] for k in range(0, len(record_1_DZ), 2)]))
        df4 = pd.DataFrame(np.array([[zscore_1_DZ.iloc[k][feature], zscore_2_DZ.iloc[k+1][feature]] for k in range(0, len(record_1_DZ), 2)]))
        df5 = pd.DataFrame(np.array([[zscore_2_DZ.iloc[k][feature], zscore_1_DZ.iloc[k+1][feature]] for k in range(0, len(record_1_DZ), 2)]))
        df6 = pd.DataFrame(np.array([[zscore_1_DZ.iloc[k][feature], zscore_3_DZ.iloc[k+1][feature]] for k in range(0, len(record_1_DZ), 2)]))
        df7 = pd.DataFrame(np.array([[zscore_3_DZ.iloc[k][feature], zscore_1_DZ.iloc[k+1][feature]] for k in range(0, len(record_1_DZ), 2)]))
        df8 = pd.DataFrame(np.array([[zscore_2_DZ.iloc[k][feature], zscore_3_DZ.iloc[k+1][feature]] for k in range(0, len(record_1_DZ), 2)]))
        df9 = pd.DataFrame(np.array([[zscore_3_DZ.iloc[k][feature], zscore_2_DZ.iloc[k+1][feature]] for k in range(0, len(record_1_DZ), 2)]))
        for j, df_j in enumerate([df1, df2, df3, df4, df5, df6, df7, df8, df9]):
            if correlation_type == "icc":
                corr_DZ_sub[j][i] = compute_icc(df_j, n = len(record_1_MZ), k = n_measurements)
            elif correlation_type == "pearson" :
                corr_DZ_sub[j][i] = pearsonr(df_j[df_j.columns[0]], df_j[df_j.columns[1]])[0]
        df = pd.concat([df1, df2, df3, df4, df5, df6, df7, df8, df9], axis = 0).rename(columns={0:"Twin A", 1: "Twin B"})
        if correlation_type == "icc":
            corr_DZ[i] = compute_icc(df, n = n_subs, k = n_measurements)
        elif correlation_type == "pearson" :
            corr_DZ[i] = pearsonr(df[df.columns[0]], df[df.columns[1]])[0]

    corr_DZ = np.reshape(corr_DZ, (n_ROI, n_freqs))
    corr_DZ = pd.DataFrame(corr_DZ, columns=freqs, index=list_ROI)
    corr_DZ.to_csv(save_file, index_label="ROI")

    for j in range(9):
        save_file = "corr_DZ_partial_" + str(j) + ".csv"
        save_file = os.path.join(main_folder_results, save_file)
        corr_DZ = np.reshape(corr_DZ_sub[j], (n_ROI, n_freqs))
        corr_DZ = pd.DataFrame(corr_DZ, columns=freqs, index=list_ROI)
        corr_DZ.to_csv(save_file, index_label="ROI")


    # Heritability
    corr_MZ = pd.read_csv(os.path.join(main_folder_results, "corr_MZ.csv"), index_col="ROI")
    corr_DZ = pd.read_csv(os.path.join(main_folder_results, "corr_DZ.csv"), index_col="ROI")
    heritability = 2*(corr_MZ - corr_DZ)

    save_file = os.path.join(main_folder_results, "heritability.csv")
    heritability.to_csv(save_file, index_label="ROI")

    heritability_avg_per_band = heritability.copy()
    for ind, band in enumerate(bands):
        cols = [c for c in heritability.columns if float(c) >= band[0] and float(c) < band[1]]
        heritability_avg_per_band[bands_names[ind]] = np.mean(heritability[cols], axis = 1)

    heritability_avg_per_band = heritability_avg_per_band[bands_names]
    heritability_avg_per_band["total_avg"] = np.mean(heritability_avg_per_band[bands_names], axis = 1)

    save_file = os.path.join(main_folder_results, "heritability_avg_per_band.csv")
    heritability_avg_per_band.to_csv(save_file, index_label="ROI")

    ### COMPUTE HERITABILITY MEAN AND STD HERE !!!!
    corr_MZ = []
    corr_DZ = []

    for j in range(9):
        save_file = "corr_MZ_partial_" + str(j) + ".csv"
        save_file = os.path.join(main_folder_results, save_file)
        corr_MZ.append(pd.read_csv(save_file, index_col="ROI"))

        save_file = "corr_DZ_partial_" + str(j) + ".csv"
        save_file = os.path.join(main_folder_results, save_file)
        corr_DZ.append(pd.read_csv(save_file, index_col="ROI"))

    heritability = 2*(np.array(corr_MZ) - np.array(corr_DZ))

    heritability_mean = pd.DataFrame(np.mean(heritability, axis = 0), columns=freqs, index=list_ROI)
    heritability_std = pd.DataFrame(np.std(heritability, axis = 0), columns=freqs, index=list_ROI)

    save_file = os.path.join(main_folder_results, "heritability_mean.csv")
    heritability_mean.to_csv(save_file, index_label="ROI")

    save_file = os.path.join(main_folder_results, "heritability_std.csv")
    heritability_std.to_csv(save_file, index_label="ROI")

    heritability_avg_per_band = heritability_mean.copy()
    for ind, band in enumerate(bands):
        cols = [c for c in heritability_mean.columns if float(c) >= band[0] and float(c) < band[1]]
        heritability_avg_per_band[bands_names[ind]] = np.mean(heritability_mean[cols], axis = 1)

    heritability_avg_per_band = heritability_avg_per_band[bands_names]
    heritability_avg_per_band["total_avg"] = np.mean(heritability_avg_per_band[bands_names], axis = 1)

    save_file = os.path.join(main_folder_results, "Heritability_mean_avg_per_band.csv")
    heritability_avg_per_band.to_csv(save_file, index_label="ROI")


if __name__ == "__main__":
    main()