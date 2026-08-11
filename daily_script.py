# -*- coding: utf-8 -*-
"""
Created on Tue May  5 16:31:34 2026

@author: astrid.nicolini
"""
# %% 
import numpy as np
import physion
from physion.analysis.read_NWB import scan_folder_for_NWBfiles, Data
import matplotlib.pyplot as plt

#%%

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
data.build_running()

Episodes = EpisodeData(data, 
                        quantities=['dFoF','running'], 
                        protocol_name = data.protocols[0],
                        verbose=False)
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


