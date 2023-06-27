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

def main(data_path = DATA_PATH, main_folder_results = FOLDER_RESULTS, only_gt = ONLY_GT):

    # --- PATH TO SAVE THE RESULTS ---

    if not os.path.exists(main_folder_results):
        os.mkdir(main_folder_results)

    if only_gt:
        main_folder_results = os.path.join(main_folder_results, "Only_GT")
    else :
        main_folder_results = os.path.join(main_folder_results, "Include_SR")
    if not os.path.exists(main_folder_results):
        os.mkdir(main_folder_results)

    main_folder_results = os.path.join(main_folder_results, "PSD_Correlations")
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
    ids_uncertain_SR = []
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

                ids_uncertain_SR.append(twins["subjects"][1])

    # Lists of IDs depending on the category of the subject
    ids_MZ = [k for k, v in rename_twins.items() if "Twin_MZ" in v]
    ids_DZ = [k for k, v in rename_twins.items() if "Twin_DZ" in v]
    ids_NT = [k for k, v in rename_twins.items() if "NotTwin" in v]
    ids_NT_only_GT = [id for id in ids_NT if id not in ids_uncertain_SR]

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
        Plot the correlation matrix if plot = True.
        """

        # Compute the Pearson correlation between every pair of subjects
        corr = np.empty(shape=(dataset_1.shape[0], dataset_2.shape[0]))
        for i, subj_1 in enumerate(dataset_1.index):
            for j, subj_2 in enumerate(dataset_2.index):
                corr[i, j] = pearsonr(dataset_1.loc[subj_1], dataset_2.loc[subj_2])[0]

        # Convert to Dataframe 
        corr = pd.DataFrame(corr, index = [rename_twins[id] for id in dataset_1.index], columns= [rename_twins[id] for id in dataset_2.index])
        if save_csv:
            save_csv = save_csv + ".csv"
            corr.to_csv(save_csv, index = True, index_label="Subjects")
        # Plot
        if plot : 
            f, ax = plt.subplots(figsize=(11, 9))
            cmap = sns.color_palette("coolwarm", as_cmap=True)
            sns.heatmap(corr, cmap=cmap,
                square=True, linewidths=.7, cbar_kws={"shrink": 1.0, "label": "Correlation"}) 
            if save_plot :
                save_plot = save_plot + ".pdf"
                plt.savefig(save_plot, format = "pdf", bbox_inches="tight")

        return corr


    def plot_autocorr_vs_crosscorr(corr_MZ, corr_DZ, corr_NT, corr_only_NT, save_plot = "test"):
        """
        Given the correlation matrices between:
        - MZ twins (ordered by pairs)
        - DZ twins (also ordered by pair)
        - NotTwins,
        it plots the histogram of the autocorrelation and crosscorrelation (twin A with twin B) for each category.
        The cross correlation for NT is just the correlation between two different subjects.

        ISSUE : Values being too close for the seaborn plot, correlation scores have to be normalized ...
        """
        # Autocorrelations
        autocorr_MZ = corr_MZ.to_numpy().diagonal()
        autocorr_DZ = corr_DZ.to_numpy().diagonal()
        autocorr_NT = corr_only_NT.to_numpy().diagonal()

        # Cross-correlations
        # [i, i+1] because pairs come together
        cross_corr_MZ = np.concatenate([[corr_MZ.to_numpy()[i, i+1], corr_MZ.to_numpy()[i+1, i]] for i in range(0, len(corr_MZ), 2)])
        cross_corr_DZ = np.concatenate([[corr_DZ.to_numpy()[i, i+1], corr_DZ.to_numpy()[i+1, i]] for i in range(0, len(corr_DZ), 2)])
        # For NT, as there are no correlation to find, we can just take the correlation with everyone, as a reference to the other ones
        cross_corr_NT = np.concatenate([np.concatenate([corr_NT.to_numpy()[i, :i], corr_NT.to_numpy()[i, i+1:]]) for i in range(0, len(corr_NT))]) 

        # Final Plot
        fig, (ax0, ax1) = plt.subplots(nrows=1, ncols= 2, figsize = (20, 9))

        ax0.hist([autocorr_MZ, autocorr_DZ, autocorr_NT], bins=30, alpha = 0.5, label=["MZ", "DZ", "NT"], histtype="stepfilled", density = True)
        ax0.legend(loc = "upper right")
        ax0.set_xlabel("Correlation")
        ax0.set_ylabel("Density") 
        ax0.set_title("Autocorrelation", fontweight = "bold", size = 16)   

        ax1.hist([cross_corr_MZ, cross_corr_DZ, cross_corr_NT], bins=30, alpha = 0.5, label=["MZ", "DZ", "NT"], histtype="stepfilled", density = True)
        ax1.legend(loc = "upper right")
        ax1.set_xlabel("Correlation")
        ax1.set_ylabel("Density")
        ax1.set_title("Cross-correlation", fontweight = "bold", size = 16)   

        save_plot = save_plot + ".pdf"
        plt.savefig(save_plot, format = "pdf", bbox_inches="tight")


    def eval_correlations(corr_MZ, corr_DZ, corr_NT, corr_only_NT, save_csv = None):
        """
        Note : corr_MZ and corr_DZ have to be ordered such that twins always are one next to each other
        """
        
        # Autocorrelation
        autocorr_MZ = corr_MZ.to_numpy().diagonal()
        autocorr_DZ = corr_DZ.to_numpy().diagonal()
        autocorr_NT = corr_only_NT.to_numpy().diagonal()

        autocorr = np.concatenate([autocorr_MZ, autocorr_DZ, autocorr_NT])
        autocorr_moy = np.mean(autocorr)
        autocorr_std = np.std(autocorr)

        # Cross correlation
        cross_corr_MZ = np.concatenate([[corr_MZ.to_numpy()[i, i+1], corr_MZ.to_numpy()[i+1, i]] for i in range(0, len(corr_MZ), 2)])
        cross_corr_DZ = np.concatenate([[corr_DZ.to_numpy()[i, i+1], corr_DZ.to_numpy()[i+1, i]] for i in range(0, len(corr_DZ), 2)])
        cross_corr_NT = np.concatenate([corr_NT.drop(columns = ["NotTwin" + str(i)]).iloc[i-1].to_numpy() for i in range(1, len(ids_NT) + 1)])

        # Save all the details
        if save_csv :
            df = pd.DataFrame({"Autocorr_MZ" : pd.Series(autocorr_MZ), "Autocorr_DZ" : pd.Series(autocorr_DZ), "Autocorr_NT" : pd.Series(autocorr_NT) , "Crosscorr MZ" : pd.Series(cross_corr_MZ), "Crosscorr DZ" : pd.Series(cross_corr_DZ), "Crosscorr NT" : pd.Series(cross_corr_NT)})
            save_csv = save_csv + ".csv"
            df.to_csv(save_csv)

        cross_corr_MZ_moy = np.mean(cross_corr_MZ)
        cross_corr_MZ_std = np.std(cross_corr_MZ)
        cross_corr_DZ_moy = np.mean(cross_corr_DZ)
        cross_corr_DZ_std = np.std(cross_corr_DZ)
        cross_corr_NT_moy = np.mean(cross_corr_NT)
        cross_corr_NT_std = np.std(cross_corr_NT)

        return autocorr_moy, cross_corr_MZ_moy, cross_corr_DZ_moy, cross_corr_NT_moy, autocorr_std, cross_corr_MZ_std, cross_corr_DZ_std, cross_corr_NT_std

    # ---   MAIN    ---

    df_final = pd.DataFrame()
    save_final = "All_correlation_every_freq.csv"
    save_final = os.path.join(main_folder_results, save_final)

    for i in range(1, 4):
        for j in range(1, 4):

            DATASET_1 = i 
            DATASET_2 = j

            records = {1 : record_1, 2 : record_2, 3 : record_3}
            record_A = records[DATASET_1]
            record_B = records[DATASET_2]

            folder_result = "Dataset_" + str(i) + "_VS_Dataset" + str(j)
            print("Current match : ", folder_result)
            folder_result = os.path.join(main_folder_results, folder_result)

            if not os.path.exists(folder_result):
                os.mkdir(folder_result)


            # Correlation every subject
            save_path = "Correlation_every_subject"
            save_path = os.path.join(folder_result, save_path)
            corr_subjects = corr_multi_subjects(record_A, record_B, save_plot=save_path, save_csv=save_path)

            # Compute Z-Score
            save_path = "Z_score_every_subject.csv"
            save_path = os.path.join(folder_result, save_path)
            df_zcore = pd.DataFrame(np.array([corr_subjects.index, zscore(corr_subjects, axis = 1).to_numpy().diagonal()]).T, columns = ["ID", "Z_score"]).set_index("ID")
            df_zcore.to_csv(save_path, index="ID")

            # Correlations between the MZ twins

            mask_MZ = [True if subj_id in ids_MZ else False for subj_id in record_A.index]
            record_A_MZ = record_A[mask_MZ].reindex(ids_MZ)
            record_B_MZ = record_B[mask_MZ].reindex(ids_MZ)


            save_path = "Correlation_MZ"
            save_path = os.path.join(folder_result, save_path)
            corr_MZ = corr_multi_subjects(record_A_MZ, record_B_MZ, save_plot=save_path, save_csv=save_path)


            # Correlations between the DZ twins

            mask_DZ = [True if subj_id in ids_DZ else False for subj_id in record_A.index]
            record_A_DZ = record_A[mask_DZ].reindex(ids_DZ)
            record_B_DZ = record_B[mask_DZ].reindex(ids_DZ)

            save_path = "Correlation_DZ"
            save_path = os.path.join(folder_result, save_path)
            corr_DZ = corr_multi_subjects(record_A_DZ, record_B_DZ, save_plot=save_path, save_csv=save_path)


            # Correlations between the other unrelated subjects
            if only_gt :
                ids_NT = ids_NT_only_GT

            mask_NT = [True if subj_id in ids_NT else False for subj_id in record_A.index]
            record_A_NT = record_A[mask_NT].reindex(ids_NT)
            record_B_NT = record_B[mask_NT].reindex(ids_NT)

            save_path = "Correlation_only_NT"
            save_path = os.path.join(folder_result, save_path)
            corr_only_NT = corr_multi_subjects(record_A_NT, record_B_NT, save_plot=save_path, save_csv=save_path)


            # Correlations between all unrelated subjects 
            mask_NT = [True if subj_id in ids_NT else False for subj_id in record_A.index]
            record_A_NT = record_A[mask_NT].reindex(ids_NT)

            if ONLY_GT :
                mask_uncertain_NT = [True if subj_id not in ids_uncertain_SR else False for subj_id in record_A.index]
                record_B_NT = record_B[mask_uncertain_NT]
            else :
                record_B_NT = record_B.copy()

            save_path = "Correlation_NT"
            save_path = os.path.join(folder_result, save_path)
            corr_NT = corr_multi_subjects(record_A_NT, record_B_NT, save_plot=save_path, save_csv=save_path)


            # Summary of the correlation distributions

            save_path = "Autocorr_and_Crosscorr_distributions"
            save_path = os.path.join(folder_result, save_path)
            plot_autocorr_vs_crosscorr(corr_MZ, corr_DZ, corr_NT, corr_only_NT, save_plot=save_path)


            # Correlations mean and std
            df = []
            bands = [BROADBAND, DELTA, THETA, ALPHA, BETA, GAMMA, HIGH_GAMMA]
            bands_names = ["BROADBAND", "DELTA", "THETA", "ALPHA", "BETA", "GAMMA", "HIGH GAMMA"]


            for k, band in tqdm(enumerate(bands), total = len(bands)) :

                columns_band_freq = [c for c in record_A.columns if float(re.search("[0-9]+\.[0-9]*", c).group(0)) >= band[0] and float(re.search("[0-9]+\.[0-9]*", c).group(0)) < band[1]]
                record_A_band_freq = record_A[columns_band_freq]
                record_B_band_freq = record_B[columns_band_freq]

                # MZ
                mask_MZ = [True if subj_id in ids_MZ else False for subj_id in record_A.index]
                record_A_MZ_band_freq = record_A_band_freq[mask_MZ].reindex(ids_MZ)
                record_B_MZ_band_freq = record_B_band_freq[mask_MZ].reindex(ids_MZ)

                corr_MZ_band_freq = corr_multi_subjects(record_A_MZ_band_freq, record_B_MZ_band_freq, plot = False)

                # DZ
                mask_DZ = [True if subj_id in ids_DZ else False for subj_id in record_A.index]
                record_A_DZ_band_freq = record_A_band_freq[mask_DZ].reindex(ids_DZ)
                record_B_DZ_band_freq = record_B_band_freq[mask_DZ].reindex(ids_DZ)

                corr_DZ_band_freq = corr_multi_subjects(record_A_DZ_band_freq, record_B_DZ_band_freq, plot = False)

                # NT
                mask_NT = [True if subj_id in ids_NT else False for subj_id in record_A.index]
                record_A_NT_band_freq = record_A_band_freq[mask_NT].reindex(ids_NT)
                record_B_NT_band_freq = record_B_band_freq[mask_NT].reindex(ids_NT)

                corr_only_NT_band_freq = corr_multi_subjects(record_A_NT_band_freq, record_B_NT_band_freq, plot = False)
                corr_NT_band_freq = corr_multi_subjects(record_A_NT_band_freq, record_B_band_freq, plot = False)

                # Final
                save_path = "All_correlations_per_class_band_" + bands_names[k]
                save_path = os.path.join(folder_result, save_path)
                autocorr_moy, cross_corr_MZ_moy, cross_corr_DZ_moy, cross_corr_NT_moy, autocorr_std, cross_corr_MZ_std, cross_corr_DZ_std, cross_corr_NT_std  = eval_correlations(corr_MZ_band_freq, corr_DZ_band_freq, corr_NT_band_freq, corr_only_NT_band_freq, save_csv=save_path)
                if i == j :
                    autocorr_moy, autocorr_std = np.nan, np.nan
                df.append([autocorr_moy, cross_corr_MZ_moy, cross_corr_DZ_moy, cross_corr_NT_moy, autocorr_std, cross_corr_MZ_std, cross_corr_DZ_std, cross_corr_NT_std])

            df = pd.DataFrame(df, index=bands_names, columns=["Autocorr", "Crosscorr MZ", "Crosscorr DZ", "Crosscorr NT", "Autocorr std", "Crosscorr MZ std", "Crosscorr DZ std", "Crosscorr NT std"])

            save_path = "Autocorr_Crosscorr_avg_std_per_freq_band.csv"
            save_path = os.path.join(folder_result, save_path)
            df.to_csv(save_path, index = True, index_label="Frequencies_Band")

            
            # Finally, we merge all the data collected for every band in one dataframe

            df_join = pd.DataFrame()
            for band in bands_names:
                name = "All_correlations_per_class_band_" + band + ".csv"
                name = os.path.join(folder_result, name)
                df = pd.read_csv(name, index_col=0)
                new_columns_map = {col_name : col_name + "_" + band for col_name in df.columns}
                df.rename(columns=new_columns_map, inplace=True)
                df_join = pd.concat([df_join, df], axis = 1)

            save_path = "All_correlations_per_class_band_merge.csv"
            save_path = os.path.join(folder_result, save_path)
            df_join.to_csv(save_path)            

            # And to finish we merge that to the final df containing all the data
            if i <= j : # Else results appear twice
                df_final = pd.concat([df_final, df_join], ignore_index = True, sort = False)

    df_final.to_csv(save_final)

    # Stack dataframe to plot in R
    all_corr = pd.DataFrame(df_final.stack(dropna=True)).reset_index().rename(columns={"level_1" : "columns", 0 : "values"})
    all_corr["columns"] = all_corr["columns"].apply(lambda x : re.split(" |_", x))
    all_corr["TypeCorr"] = all_corr["columns"].apply(lambda x : x[0])
    all_corr["TwinType"] = all_corr["columns"].apply(lambda x : x[1])
    all_corr["FreqBand"] = all_corr["columns"].apply(lambda x : " ".join(x[2:]))

    all_corr.drop(columns=["level_0", "columns"], inplace = True)

    save_corr_path = os.path.join(main_folder_results, "All_correlation_every_freq_stacked.csv")
    all_corr.to_csv(save_corr_path)

if __name__ == "__main__":
    main()