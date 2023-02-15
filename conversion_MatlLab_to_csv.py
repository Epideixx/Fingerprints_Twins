import os
import pandas as pd
import scipy.io
import re

import matplotlib.pyplot as plt
import numpy as np

PATH = "Data/HCP_Data"
SAVE_DIR = "Data"
LENGTH_TO_KEEP_ROI = 300

records = {1 : [], 2 : [], 3 : []}

for root, dirs, files in os.walk(PATH):
    for i, file in enumerate(files) :
        subject_id = int(re.search("[0-9]+_", file).group(0)[:-1])
        recording_id = re.search("recording_[1-3]", file).group(0)
        recording_id = int(recording_id[-1])
        path_file = os.path.join(PATH, file)
        
        data = scipy.io.loadmat(path_file)["TF"].squeeze()
        data = data[:, :LENGTH_TO_KEEP_ROI]
        data = data.reshape(-1)

        records[recording_id].append({"subject_id" : subject_id, "signal" : data})

        if i == 0 :
            freqs = scipy.io.loadmat(path_file)["Freqs"][0][:LENGTH_TO_KEEP_ROI]
            ROIs = scipy.io.loadmat(path_file)["RowNames"][0]
            columns = np.concatenate([[roi[0] + "_" + str(round(freq, 2)) for freq in freqs] for roi in ROIs])

record_1 = pd.DataFrame(data = [rec["signal"] for rec in records[1]], columns=columns, index=[rec["subject_id"] for rec in records[1]])
record_2 = pd.DataFrame(data = [rec["signal"] for rec in records[2]], columns=columns, index=[rec["subject_id"] for rec in records[2]])
record_3 = pd.DataFrame(data = [rec["signal"] for rec in records[3]], columns=columns, index=[rec["subject_id"] for rec in records[3]])

record_1.to_csv(os.path.join(SAVE_DIR, "record_1.csv"), index_label="Subject_ID")
record_2.to_csv(os.path.join(SAVE_DIR, "record_2.csv"), index_label="Subject_ID")
record_3.to_csv(os.path.join(SAVE_DIR, "record_3.csv"), index_label="Subject_ID")



