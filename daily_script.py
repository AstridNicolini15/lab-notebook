# -*- coding: utf-8 -*-
"""
Created on Tue May  5 16:31:34 2026

@author: astrid.nicolini
"""
# %% 
import numpy as np
import os, sys 
sys.path.insert(0,'/home/user/lab-notebook/astrid/physion/src')
#sys.path = [sys.path[i] for i in range(-1,len(sys.path)-1)]
from physion.analysis.read_NWB import scan_folder_for_NWBfiles, Data
from physion.analysis.episodes.build import *
import matplotlib.pyplot as plt
import physion

#%% TEMPLATES

filename = '/home/user/DATA/Astrid/Cibele_data/PV-cells_WT_Adult_V1/NWBs/2024_03_22-14-35-26.nwb' #contrasts
filename = '/home/user/DATA/Astrid/Cibele_data/SST-cells_cond-GluN1-KO_Adult_V1/NWBs/2026_03_06-15-14-09.nwb' #orientations

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
#data.build_Deconvolved(quantity = 'dFoF')
data.build_running()

Episodes = EpisodeData(data, 
                        quantities=['dFoF', 'running'], 
                        protocol_name = data.protocols[0],
                        verbose=False)

#%% allen


filename = Tunings_ko[11]['datafile']
allen_data = physion.analysis.read_NWB.Data(filename, verbose=False) 
#data = Data(filename, verbose=False)
#protocol_name=[p for p in data.protocols if '8orientation' in p][0]
dFoF_parameters = dict(\
        roi_to_neuropil_fluo_inclusion_factor=1.15,
        neuropil_correction_factor = 0.7,
        roi_specific_neuropil_correction_factor = True,
        method_for_F0 = 'sliding_percentile',
        percentile=5., # percent
        sliding_window = 5*60, # seconds
)

allen_data.build_dFoF(**dFoF_parameters, verbose=False)

Episodes = EpisodeData(allen_data, 
                        quantities=['dFoF'], 
                        protocol_name = allen_data.protocols[0],
                        verbose=False)

#%%
start = 5000
stop = 5100
validroi = 0
plt.figure(figsize = (7,7))
plt.plot(data.rawFluo[validroi,start:stop], color = 'red')
plt.plot(data.dFoF[validroi,start:stop])

plt.plot(allen_data.dFoF[data.valid_roiIndices[validroi],start:stop], color = 'black')


#%%

filename = '/home/user/DATA/Astrid/Cibele_data/summary/allen_neuropil_Tunings_SST-cells_WT_Adult_V1_contrast-1.0.npy'
allen_Tunings = np.load(filename, allow_pickle = True)


filename = '/home/user/DATA/Astrid/Cibele_data/summary/Tunings_SST-cells_WT_Adult_V1_contrast-1.0.npy'
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
summary_folder = '/home/user/DATA/Astrid/Cibele_data/summary/'

Tunings_wt_1 = np.load(summary_folder +  'Tunings_SST-cells_WT_Adult_V1_contrast-1.0.npy', allow_pickle = True)
Tunings_wt_05 = np.load(summary_folder +  'Tunings_SST-cells_WT_Adult_V1_contrast-0.5.npy', allow_pickle = True)
Tunings_ko_1 = np.load(summary_folder +  'Tunings_SST-cells_cond-GluN1-KO_Adult_V1_contrast-1.0.npy', allow_pickle = True)
Tunings_ko_05 = np.load(summary_folder +  'Tunings_SST-cells_cond-GluN1-KO_Adult_V1_contrast-0.5.npy', allow_pickle = True)
