#%%
import sys 
import numpy as np
from scipy import stats
import pandas as pd 

sys.path += ['./physion/src']
from physion.utils import plot_tools as pt
from physion.analysis.read_NWB\
                         import scan_folder_for_NWBfiles, Data
from physion.analysis.episodes.build import EpisodeData
from physion.analysis.protocols.contrast_sensitivity\
                        import compute_sensitivity_per_cells

from physion.analysis. protocols.contrast_sensitivity import get_responses, get_gains


#%%
#------------Build functions----------------#


def fill_missing_contrast_values_with_nans(Sensitivity, Episodes): 
    corrected_sensitivity = {}
    contrast_values = np.unique(Episodes.contrast)
    nROIs = Episodes.data.nROIs

    for key in Sensitivity.keys() : 
        corrected_sensitivity [key] = []

    for contrast in contrast_values : 
        if contrast in Sensitivity['contrast'] : 
            ind_contrast = np.argwhere(Sensitivity['contrast'] == contrast)[0][0]
            for key in Sensitivity.keys() :
                if key != 'contrast' and key != 'ntrials': 
                    corrected_sensitivity[key].append(Sensitivity[key][:,ind_contrast])
                if key == 'ntrials' : 
                    corrected_sensitivity[key].append(Sensitivity[key][ind_contrast])
        else : 
            for key in Sensitivity.keys() : 
                if key != 'contrast': 
                    if 'significant' in key : 
                        corrected_sensitivity[key].append([False]*nROIs)
                    elif key == 'ntrials' : 
                        corrected_sensitivity[key].append(np.nan)
                    else : 
                        corrected_sensitivity[key].append([np.nan]*nROIs)

    #reshape correctly 
    for key in Sensitivity.keys() : 
        if key != 'contrast' : 
            corrected_sensitivity[key]  = np.array(corrected_sensitivity[key]).T
    corrected_sensitivity['contrast'] = contrast_values

    return corrected_sensitivity
