# ---   DEPENDENCIES   ---

import os
import re
import numpy as np
import pandas as pd
from scipy.stats import pearsonr
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
DATA_PATH = "Data/Schaefer_30_second"
FOLDER_RESULTS = "Results_Log_Schaefer_test"
ONLY_GT = True


# ---   MAIN FUNCTION   ---
def main(data_path = DATA_PATH, main_folder_results = FOLDER_RESULTS, only_gt = ONLY_GT):

    """
    Compute singleton and twin pairs differentiation accuracies using PSDs (see paper) computed from 30-second segments of recordings.
    It computes it separately for pairs of recordings from the same session and for pairs of recordings from different sessions.

    Inputs
    ------
    data_path : string
        Path of the folder containing the recordings (Scaefer_PSD_m1_30second_seg1.csv, Scaefer_PSD_m1_30second_seg2.csv, ...)
    main_folder_results : string
        Path folder to save results
    only_gt : boolean
        True if use only genetic test, False if also includes self-reported zygosity
    """

    # --- PATH TO SAVE THE RESULTS ---

    if not os.path.exists(main_folder_results):
        os.mkdir(main_folder_results)

    if only_gt:
        main_folder_results = os.path.join(main_folder_results, "Only_GT")
    else :
        main_folder_results = os.path.join(main_folder_results, "Include_SR")
    if not os.path.exists(main_folder_results):
        os.mkdir(main_folder_results)

    main_folder_results = os.path.join(main_folder_results, "PSD_Accuracy_30_sec")
    if not os.path.exists(main_folder_results):
        os.mkdir(main_folder_results)

    # ---   IMPORT DATA   ---

    # Import the csv and apply log to it
    record_1 = []
    for path_record_1 in [os.path.join(data_path, "Scaefer_PSD_m1_30second_seg" + str(i) + ".csv") for i in range(1, 9)] :
        record_1.append(np.log(pd.read_csv(path_record_1, index_col="Subject_ID")))

    print("Every PSD from session 1 has been loaded")

    record_2 = []
    for path_record_2 in [os.path.join(data_path, "Scaefer_PSD_m2_30second_seg" + str(i) + ".csv") for i in range(1, 9)] :
        record_2.append(np.log(pd.read_csv(path_record_2, index_col="Subject_ID")))

    print("Every PSD from session 2 has been loaded")

    record_3 = []
    for path_record_3 in [os.path.join(data_path, "Scaefer_PSD_m3_30second_seg" + str(i) + ".csv") for i in range(1, 9)] :
        record_3.append(np.log(pd.read_csv(path_record_3, index_col="Subject_ID")))

    print("Every PSD from session 3 has been loaded")



    # ---   USEFULL ITEMS   ---

    ### Annotate Data such that twins are linked ###

    # Import extra data (confidential)
    all_data_restricted_filename = "Data/All_Data_RESTRICTED.csv"
    all_data_restricted = pd.read_csv(all_data_restricted_filename)

    # Remove subjects who don't have a MEG recording
    subjects_with_MEG = record_1[0].index
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
    subjects_row_to_id = {i : id for i, id in enumerate(record_1[0].index)}
    subject_id_to_row = {id : i for i, id in enumerate(record_1[0].index)}

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

    # Summary 

    nb_pair = len([pair["subjects"] for pair in list(twins_dict.values()) if len(pair["subjects"]) >= 2])
    nb_pair_MZ = len(ids_MZ)//2
    nb_pair_DZ = len(ids_DZ)//2
    nb_not_twins = len(ids_NT)
    nb_pair_not_twins = len([pair["subjects"] for pair in list(twins_dict.values()) if len(pair["subjects"]) >= 2 and pair["type"] == "NT"])

    print("The total number of participants is : ", len(all_data_restricted))
    print("We will work with {} pairs of siblings : {} pairs of MZ twins, {} pairs of DZ twins and {} pairs of siblings with uncertain relation".format(nb_pair, nb_pair_MZ, nb_pair_DZ, nb_pair_not_twins))



    # ---   USEFULL FUNCTIONS   ---
    def corr_multi_subjects(dataset_1, dataset_2, plot = True, save_plot = None, save_csv = None):
        """
        Compute the pearson correlation of the dataset_1 with the dataset_2, subjects by subjects (= rows in datasets).
        Plot the correlation matrix if plot = True. Save the correlation matrix as a csv file if save_file = True.
        """

        # Compute the Pearson correlation between every pair of subjects
        corr = np.empty(shape=(dataset_1.shape[0], dataset_2.shape[0]))
        for i, subj_1 in enumerate(dataset_1.index):
            for j, subj_2 in enumerate(dataset_2.index):
                corr[i, j] = pearsonr(dataset_1.loc[subj_1], dataset_2.loc[subj_2])[0]

        # Convert to Dataframe 
        corr = pd.DataFrame(corr, index = [rename_twins[id] for id in dataset_1.index], columns= [rename_twins[id] for id in dataset_2.index])

        # Save csv if asked
        if save_csv:
            save_csv = save_csv + ".csv"
            corr.to_csv(save_csv, index = True, index_label="Subjects")

        # Plot if asked
        if plot : 
            f, ax = plt.subplots(figsize=(11, 9))
            cmap = sns.color_palette("coolwarm", as_cmap=True)
            sns.heatmap(corr, cmap=cmap,
                square=True, linewidths=.7, cbar_kws={"shrink": 1.0, "label": "Correlation"}) 
            if save_plot :
                save_plot = save_plot + ".pdf"
                plt.savefig(save_plot, format = "pdf", bbox_inches="tight")

        return corr


    def eval_accuracy(corr_df):
        """
        Using the correlation dataframe, it computes the accuracy for the fingerprint prediction, for the MZ twin prediction and for the DZ twin prediction.
        """

        # Autocorrelation accuracy
        nb_true_pred = corr_df.to_numpy().argmax(axis = 0) == np.arange(len(corr_df))
        nb_true_pred = np.count_nonzero(nb_true_pred)
        auto_accuracy = nb_true_pred/len(corr_df)


        # Crosscorrelation accuracy MZ
        mask_MZ = [True if "Twin_MZ" in name else False for name in corr_df.index]
        corr_df_rows_MZ = corr_df[mask_MZ]

        second_best_corr = corr_df_rows_MZ.to_numpy().argsort(axis = 1)[:, ::-1][:, 1]
        nb_true_pred = 0
        for i, num_column in enumerate(second_best_corr):
            match = re.search("Twin_MZ_[0-9]+(A|B)", corr_df_rows_MZ.index[i]).group(0)
            if corr_df.columns[num_column][:-1] == match[:-1] and match[-1] != corr_df.columns[num_column][-1]:
                nb_true_pred += 1
        cross_accuracy_MZ = nb_true_pred/len(corr_df_rows_MZ)


        # Crosscorrelation accuracy DZ
        mask_DZ = [True if "Twin_DZ" in name else False for name in corr_df.index]
        corr_df_rows_DZ = corr_df[mask_DZ]

        second_best_corr = corr_df_rows_DZ.to_numpy().argsort(axis = 1)[:, ::-1][:, 1]
        nb_true_pred = 0
        for i, num_column in enumerate(second_best_corr):
            match = re.search("Twin_DZ_[0-9]+(A|B)", corr_df_rows_DZ.index[i]).group(0)
            if corr_df.columns[num_column][:-1] == match[:-1] and match[-1] != corr_df.columns[num_column][-1]:
                nb_true_pred += 1
        cross_accuracy_DZ = nb_true_pred/len(corr_df_rows_DZ)
        
        return auto_accuracy, cross_accuracy_MZ, cross_accuracy_DZ


    # ---   MAIN    ---

    # ACCURACY WITHIN A RECORD ###


    group_of_records = {1 : record_1, 2 : record_2, 3 : record_3}
    folder_withing_record = os.path.join(main_folder_results, "Within_record")

    # Final dataframe containing all the results
    df_all_merge_within_records = []

    # For every session, we compute the accuracies 
    for id_record in range(1, 4):
        folder_ID_result = os.path.join(folder_withing_record, "Session" + str(id_record))

        if not os.path.exists(folder_ID_result):
            os.makedirs(folder_ID_result)

        # Dataframe containing all the accuracies for the session
        df_final = pd.DataFrame()
        save_final = "All_accuracies_every_freq.csv"
        save_final = os.path.join(folder_ID_result, save_final)

        # For every pair of recordings
        for i in range(1, 9):
            for j in range(1, 9):

                DATASET_1 = i 
                DATASET_2 = j

                record_A = group_of_records[id_record][DATASET_1 - 1]
                record_B = group_of_records[id_record][DATASET_2 - 1]

                file_result = "PSD_" + str(i) + "_VS_PSD_" + str(j) + "_accuracies.csv"
                print("Current match : ", "PSD_" + str(i) + "_VS_PSD_" + str(j))
                file_result = os.path.join(folder_ID_result, file_result)

                bands = [BROADBAND, DELTA, THETA, ALPHA, BETA, GAMMA, HIGH_GAMMA]
                bands_names = ["BROADBAND", "DELTA", "THETA", "ALPHA", "BETA", "GAMMA", "HIGH GAMMA"]

                df_every_freq = pd.DataFrame()

                # Compute the accuracy for every band (including BROADBAND)
                for k, band in tqdm(enumerate(bands), total = len(bands)):

                    # Compute correlation matrix
                    columns_band_freq = [c for c in record_A.columns if float(re.search("[0-9]+\.[0-9]*", c).group(0)) >= band[0] and float(re.search("[0-9]+\.[0-9]*", c).group(0)) < band[1]]
                    record_A_band_freq = record_A[columns_band_freq]
                    record_B_band_freq = record_B[columns_band_freq]

                    corr_subjects_band_freq = corr_multi_subjects(record_A_band_freq, record_B_band_freq, plot = False)  

                    # Compute the accuracies using the correlation matrix
                    corr = corr_subjects_band_freq.copy()
                    auto_acc, cross_MZ_acc, cross_DZ_acc = eval_accuracy(corr)

                    # If record 1 == record 2, we remove the singleton accuracy (=1)
                    if i == j :
                        auto_acc = np.nan

                    # Save the result in a dataframe
                    df = pd.DataFrame([[auto_acc, cross_MZ_acc, cross_DZ_acc]], columns=["Acc Autocorr", "Acc Crosscorr MZ", "Acc Crosscorr DZ"])
                    new_columns_map = {col_name : col_name + "_" + bands_names[k] for col_name in df.columns}
                    df.rename(columns=new_columns_map, inplace=True)
                    df_every_freq = pd.concat([df_every_freq, df], axis = 1)

                df_every_freq.to_csv(file_result)            

                # And to finish we merge that to the final df containing all the data
                df_final = pd.concat([df_final, df_every_freq], ignore_index = True, sort = False)

        df_final.to_csv(save_final)
        df_all_merge_within_records.append(df_final)

    # Save the final global result, with all the accuracies merged
    all_acc = pd.concat(df_all_merge_within_records, axis = 0)
    all_acc.to_csv(os.path.join(folder_withing_record, "All_accuracies_every_freq.csv"))

    # Create a stacked version to be used in R
    all_acc = pd.DataFrame(all_acc.stack(dropna = True)).reset_index().rename(columns={"level_1" : "columns", 0 : "values"})
    all_acc["columns"] = all_acc["columns"].apply(lambda x : re.split("_", x))
    all_acc["AccType"] = all_acc["columns"].apply(lambda x : x[0][4:])
    all_acc["FreqBand"] = all_acc["columns"].apply(lambda x : " ".join(x[1:]))
    all_acc.drop(columns=["level_0", "columns"], inplace = True)
    save_acc_path = os.path.join(folder_withing_record, "All_accuracies_every_freq_stacked.csv")
    all_acc.to_csv(save_acc_path)


    ### ACCURACY BETWEEN RECORDS ###


    group_of_records = {1 : record_1, 2 : record_2, 3 : record_3}
    folder_between_records = os.path.join(main_folder_results, "Between_records")

    # Final dataframe containing all the results
    df_all_merge_between_records = []

    # For every pair of different sessions, we compute the accuracies 
    for id_record_1 in range(1, 4):
        for id_record_2 in range(1, 4):
            if id_record_1 != id_record_2 :

                folder_ID_result = os.path.join(folder_between_records, "Session" + str(id_record_1) + "_VS_" + str(id_record_2))
                
                if not os.path.exists(folder_ID_result):
                    os.makedirs(folder_ID_result)

                # Dataframe containing all the accuracies for this pair of sessions
                df_final = pd.DataFrame()
                save_final = "All_accuracies_every_freq.csv"
                save_final = os.path.join(folder_ID_result, save_final)

                # For every pair of recordings
                for i in range(1, 9):
                    for j in range(1, 9):
                        DATASET_1 = i 
                        DATASET_2 = j

                        record_A = group_of_records[id_record_1][DATASET_1 - 1]
                        record_B = group_of_records[id_record_2][DATASET_2 - 1]

                        file_result = "PSD_" + str(i) + "_VS_PSD_" + str(j) + "_accuracies.csv"
                        print("Current match : ", "PSD_" + str(i) + "_VS_PSD_" + str(j))
                        file_result = os.path.join(folder_ID_result, file_result)


                        bands = [BROADBAND, DELTA, THETA, ALPHA, BETA, GAMMA, HIGH_GAMMA]
                        bands_names = ["BROADBAND", "DELTA", "THETA", "ALPHA", "BETA", "GAMMA", "HIGH GAMMA"]

                        df_every_freq = pd.DataFrame()

                        # Compute the accuracy for every band (including BROADBAND)
                        for k, band in tqdm(enumerate(bands), total = len(bands)):

                            df = []

                            # Compute correlation matrix
                            columns_band_freq = [c for c in record_A.columns if float(re.search("[0-9]+\.[0-9]*", c).group(0)) >= band[0] and float(re.search("[0-9]+\.[0-9]*", c).group(0)) < band[1]]
                            record_A_band_freq = record_A[columns_band_freq]
                            record_B_band_freq = record_B[columns_band_freq]

                            corr_subjects_band_freq = corr_multi_subjects(record_A_band_freq, record_B_band_freq, plot = False)
                            
                            # Compute the accuracies using the correlation matrix
                            corr = corr_subjects_band_freq.copy()
                            auto_acc, cross_MZ_acc, cross_DZ_acc = eval_accuracy(corr)
                            
                            # Save the result in a dataframe
                            df = pd.DataFrame([[auto_acc, cross_MZ_acc, cross_DZ_acc]], columns=["Acc Autocorr", "Acc Crosscorr MZ", "Acc Crosscorr DZ"])
                            new_columns_map = {col_name : col_name + "_" + bands_names[k] for col_name in df.columns}
                            df.rename(columns=new_columns_map, inplace=True)
                            df_every_freq = pd.concat([df_every_freq, df], axis = 1)

                        df_every_freq.to_csv(file_result)            

                        # And to finish we merge that to the final df containing all the data
                        df_final = pd.concat([df_final, df_every_freq], ignore_index = True, sort = False)

                df_final.to_csv(save_final)
                df_all_merge_between_records.append(df_final)

    # Save the final global result, with all the accuracies merged
    all_acc = pd.concat(df_all_merge_between_records, axis = 0)
    all_acc.to_csv(os.path.join(folder_between_records, "All_accuracies_every_freq.csv"))

    # Create a stacked version to be used in R
    all_acc = pd.DataFrame(all_acc.stack(dropna = True)).reset_index().rename(columns={"level_1" : "columns", 0 : "values"})
    all_acc["columns"] = all_acc["columns"].apply(lambda x : re.split("_", x))
    all_acc["AccType"] = all_acc["columns"].apply(lambda x : x[0][4:])
    all_acc["FreqBand"] = all_acc["columns"].apply(lambda x : " ".join(x[1:]))
    all_acc.drop(columns=["level_0", "columns"], inplace = True)
    save_acc_path = os.path.join(folder_between_records, "All_accuracies_every_freq_stacked.csv")
    all_acc.to_csv(save_acc_path)

if __name__ == "__main__":
    main()