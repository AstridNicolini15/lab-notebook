# -*- coding: utf-8 -*-
"""
Created on Tue May  5 16:31:34 2026

@author: astrid.nicolini
"""
# %% 
import numpy as np
import os, sys , shutil 
sys.path += ['/home/user/lab-notebook/astrid/physion/src']

from physion.analysis.read_NWB import scan_folder_for_NWBfiles, Data
from physion.analysis.episodes.build import *
import matplotlib.pyplot as plt
import physion

#%% TEMPLATES

filename = '/home/user/DATA/Astrid/Cibele_data/PV-cells_WT_Adult_V1/NWBs/2024_03_22-14-35-26.nwb' #contrasts
filename = '/home/user/DATA/Astrid/Cibele_data/SST-cells_cond-GluN1-KO_Adult_V1/NWBs/2026_03_06-15-14-09.nwb' #orientations


filename = '/home/user/DATA/Astrid/Cibele_data/PYR-SynGCaMP_WT_V1/NWBs/2025_10_28-15-33-31.nwb'
data = physion.analysis.read_NWB.Data(filename, verbose=False) 
#data = Data(filename, verbose=False)
#protocol_name=[p for p in data.protocols if '8orientation' in p][0]
dFoF_parameters = dict(\
        roi_to_neuropil_fluo_inclusion_factor=1.15,
        neuropil_correction_factor = 0.7,
        method_for_F0 = 'sliding_percentile',
        percentile=5., # percent
        sliding_window = 5*60, # seconds
)
print(data.protocols)
data.build_dFoF(**dFoF_parameters, verbose=False)
data.build_Deconvolved(quantity = 'dFoF')
data.build_running()

Episodes = EpisodeData(data, 
                        quantities=['dFoF','running'], 
                        protocol_name = data.protocols[0],
                        verbose=False)



#%%

filename = '/home/user/DATA/Astrid/Taddy_data/OneDrive_1_9-3-2026/summary/Tunings_Wild-Type_contrast-1.0.npy'
Tunings = np.load(filename, allow_pickle = True)


#%%

folder_path = "/home/user/DATA/Astrid/Taddy_data/"

folders = [
    "SST_WT",
    "SST_GluN1",
    "SST_GluN3"
    ]

dates = [[],[],[]]
for i, folder in enumerate(folders) : 

        DATASET_online = scan_folder_for_NWBfiles(folder_path + folder)

        for filename in DATASET_online['files'] : 

                data = physion.analysis.read_NWB.Data(filename, verbose=False) 

                if data.protocols[0] == 'ff-gratings-8orientation-2contrasts-15repeats' : 

                        dates[i].append(data.nwbfile.identifier)

#%%

dFoF_parameters = dict(\
        roi_to_neuropil_fluo_inclusion_factor=1.15,
        neuropil_correction_factor = 0.7,
        method_for_F0 = 'sliding_percentile',
        percentile=5., # percent
        sliding_window = 5*60, # seconds
)

folder_path_online = "/home/user/DATA/Astrid/Taddy_data/OneDrive_1_9-3-2026/"

folders = ['Wild-Type',
           'GluN1-KO'] 

dates_online = [[],[]]
for i, folder in enumerate(folders) : 

        DATASET_online = scan_folder_for_NWBfiles(folder_path_online + folder)

        for filename in DATASET_online['files'] : 

                data = physion.analysis.read_NWB.Data(filename, verbose=False) 

                if data.protocols[0] == 'ff-gratings-8orientation-2contrasts-15repeats' : 

                        dates_online[i].append(data.nwbfile.identifier)

#%%