#%%

import os, sys , shutil 
import multiprocessing
import numpy as np
import matplotlib.pyplot as pt


os.chdir('/home/user/lab-notebook/astrid/physion/src')
from physion.analysis.read_NWB import scan_folder_for_NWBfiles, Data
from physion.analysis.episodes.build import EpisodeData
from physion.analysis.protocols.orientation_tuning import *
from physion.assembling.dataset import read_spreadsheet

sys.path += ['/home/user/lab-notebook/astrid']
from neuropil_factors_and_SST_responsiveness.neuropil_factors_effect_functions import *
from run_rest_responses.tuning_arousal_summary_functions import *

#%% Plots (summaries with differents factors are necessary and can be created below)

ylims = [[-0.1,0.4], [-0.2,1.1],[-0.2,1.1]]
plot_tuning_responses_of_sign_and_non_sign_cells(folders, ylims, averaged_by_sessions = True) 

#effect of ROI_TO_NEUROPIL_INCLUSION_FACTOR on nb of valid rois 
neuropil_inclusion_factors = [0,0.5,0.85,1,1.15,1.3,1.5,1.75,2,2.5,3.5,6]
nb_valid_rois_by_inclusion_factor(folders, summary_folder, neuropil_inclusion_factors, neuropil_correction_factor = 0.7) 
    

#effect of NEUROPIL_CORRECTION_FACTOR on nb of valid and significant rois 
neuropil_correction_factors = [0,0.15,0.3,0.55,0.7,0.85,1]
nb_valid_and_sign_rois_by_correction_factor(folders, summary_folder, neuropil_correction_factors, neuropil_inclusion_factor = 1.15) 



for i, folder in enumerate(folders):
    plot_orientation_tuning_curve_multiple_corrfact_one_graph(['%s_contrast-1.0' % folder, 
                        '%s_contrast-0.5' % folder],
                        path=summary_folder,
                        average_by='sessions',
                        colors=[pt.tab10(i), 'lightgrey'],
                        neuropil_correction_factors = [],
                        neuropil_inclusion_factor = 1.15)
    


for i, folder in enumerate(folders):
    plot_orientation_tuning_variations_both_parameters(['%s_contrast-1.0' % folder, 
                        '%s_contrast-0.5' % folder],
                        path= summary_folder,
                        average_by='sessions',
                        colors=[pt.tab10(i), 'lightgrey'],
                        fig_args={'right':20, 'figsize' : (7,7)},
                        neuropil_inclusion_factors = [1.15,2,3],
                        neuropil_correction_factors = [0.7,0.85,1])
    