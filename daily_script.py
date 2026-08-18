# -*- coding: utf-8 -*-
"""
Created on Tue May  5 16:31:34 2026

@author: astrid.nicolini
"""
# %% 
import numpy as np
import os, sys , shutil 
sys.path += ['./physion/src']

from physion.analysis.read_NWB import scan_folder_for_NWBfiles, Data
from physion.analysis.episodes.build import *
import matplotlib.pyplot as plt
import physion

#%%

filename = '/home/user/DATA/Astrid/Cibele_data/PV-cells_WT_Adult_V1/NWBs/2024_03_22-14-35-26.nwb' #contrasts
filename = '/home/user/DATA/Astrid/Cibele_data/SST-cells_cond-GluN1-KO_Adult_V1/NWBs/2026_03_06-15-14-09.nwb' #orientations

filename = Deconvolved_Tunings[5]['datafile']
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
data.build_Deconvolved()
data.build_running()

Episodes = EpisodeData(data, 
                        quantities=['Deconvolved','running'], 
                        protocol_name = data.protocols[0],
                        verbose=False)

#%% 

dFoF_parameters = dict(\
        roi_to_neuropil_fluo_inclusion_factor=1.15,
        neuropil_correction_factor = 0.7,
        method_for_F0 = 'sliding_percentile',
        percentile=5., # percent
        sliding_window = 5*60, # seconds
        with_correctedFluo_and_F0=True,
)
data.build_dFoF(**dFoF_parameters, verbose=False)


    

F = data.dFoF 

batch_size = 500

tau = 1.3 #in seconds

fs = data.CaImaging_dt

S = oasis(F, batch_size, tau, fs)

i_start = 1000
i_stop = 1040
plt.figure(figsize=(6,6))
plt.plot(data.t_dFoF[i_start:i_stop], S[0,i_start:i_stop])
plt.plot(data.t_dFoF[i_start:i_stop], data.dFoF[0,i_start:i_stop])

#%%


Tunings = np.load('/home/user/DATA/Astrid/run_rest_summary/Tunings_SST-cells_cond-GluN1-KO_Adult_V1_contrast-1.0.npy', allow_pickle=True)
T_Responses = get_tuning_responses(Tunings, average_by='sessions') 

Tunings = np.load('/home/user/DATA/Astrid/run_rest_summary/Tunings_SST-cells_cond-GluN1-KO_Adult_V1_contrast-0.5.npy', allow_pickle=True)
T_Responses0 = get_tuning_responses(Tunings, average_by='sessions') 

Deconvolved_Tunings = np.load('/home/user/DATA/Astrid/deconvolved_summary/Deconvolved_Tunings_SST-cells_cond-GluN1-KO_Adult_V1_contrast-1.0.npy', allow_pickle=True)
D_T_Responses = get_tuning_responses(Deconvolved_Tunings, average_by='sessions') 

Deconvolved_Tunings = np.load('/home/user/DATA/Astrid/deconvolved_summary/Deconvolved_Tunings_SST-cells_cond-GluN1-KO_Adult_V1_contrast-0.5.npy', allow_pickle=True)
D_T_Responses0 = get_tuning_responses(Deconvolved_Tunings, average_by='sessions') 

plt.figure(figsize=(6,6))
plt.plot(np.mean(T_Responses, axis=0)/np.max(np.mean(T_Responses, axis=0)), label='dFoF', color = 'orange')
plt.plot(np.mean(T_Responses0, axis=0)/np.max(np.mean(T_Responses0, axis=0)), label='dFoF', color = 'lightgrey')
plt.plot(np.nanmean(D_T_Responses, axis=0)/np.max(np.nanmean(D_T_Responses, axis=0)), label='deconvolved dFoF', color = 'red')
plt.plot(np.nanmean(D_T_Responses0, axis=0)/np.max(np.nanmean(D_T_Responses0, axis=0)), label='deconvolved dFoF', color = 'grey')


#%%
Sensitivities = np.load('/home/user/DATA/Astrid/Summaries_stats/Run_Sensitivities_PV-cells_WT_Adult_V1_angle-90.0.npy', allow_pickle=True)
Responses = np.array(get_responses(Sensitivities, average_by='sessions'))

#%%
folders = [#"PV-cells_WT_Adult_V1", 
    "SST-cells_WT_Adult_V1",
    "SST-cells_cond-GluN1-KO_Adult_V1"
    ]

ylims = [(-0.05,1),(-0.15,0.8),(-0.15,0.85)] #to automatize 

i=2
colors = [[pt.tab10(1), 'lightgrey'],[pt.tab10(2), 'lightgrey']]


#%%   

    


#Mean deltaF/F for run and stationary periods 

folders = [#"PV-cells_WT_Adult_V1", 
    "SST-cells_WT_Adult_V1",
    "SST-cells_cond-GluN1-KO_Adult_V1"
    ]

k= 2
folder = folders[k]
colors = [[pt.tab10(1),'lightgrey'], [pt.tab10(2),'lightgrey']]


summary = np.load('/home/user/DATA/Astrid/run_rest_summary/' + summary_protocol + '_' + folder + '_' + protocol_control_cond + '.npy', allow_pickle=True)  
#%%
sys.path += ['/home/user/lab-notebook/astrid']
from run_rest_responses.tuning_arousal_summary_functions import *
import physion.utils.plot_tools as pt
import matplotlib.pyplot as plt

#%%

istart = 1000
istop = 1300

plt.figure(figsize=(6,6))
plt.plot(data.dFoF[10,istart:istop])
plt.plot(data.Deconvolved[10,istart:istop])