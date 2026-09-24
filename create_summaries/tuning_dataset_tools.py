#%% Tools necessary to create Tuning summaries with arousal cond, neuropil cond etc.

import os, sys , shutil 
import numpy as np
from scipy import stats
import pandas as pd 

os.chdir('/home/user/lab-notebook/astrid')
sys.path += ['./physion/src']
import physion.utils.plot_tools as pt
from physion.analysis.read_NWB\
                         import scan_folder_for_NWBfiles, Data
from physion.analysis.episodes.build import EpisodeData

from physion.analysis.protocols.orientation_tuning import *
from physion.analysis.protocols.contrast_sensitivity import *

#%%

def get_summary_prefix_name(quantity, filtering_cond_name, special_dFoF_params) :

    if quantity == 'dFoF' : 
        summary_prefix_name = ''
    elif quantity[:11] == 'Deconvolved' : 
        summary_prefix_name = quantity + '_'
    else : 
        raise ValueError("quantity should be either 'Deconvolved_something' or 'dFoF'")

    if filtering_cond_name is not None : 
        summary_prefix_name += filtering_cond_name  + '_' #Run or Rest for exemple

    if special_dFoF_params is not None : 

        summary_prefix_name += 'inclusion_factor-%s_correction_factor-%s_' % (special_dFoF_params['neuropil_inclusion_factor'], special_dFoF_params['neuropil_correction_factor'])

    return summary_prefix_name


def build_filtering_cond_quantities(filtering_cond_name, data, quantities):

    if filtering_cond_name == 'Run' or filtering_cond_name == 'Rest' :
        quantities += ['running']
        data.build_running(verbose=False)
        check_presence_locomotion_values(data)

    return data, quantities


def get_filtering_cond(filtering_cond_name, Episodes):

    if filtering_cond_name == 'Run' : 
        filtering_cond =  compute_arousal_mask(Episodes)[0]

    elif filtering_cond_name == 'Rest' : 
        filtering_cond =  compute_arousal_mask(Episodes)[1]

    return filtering_cond


def check_presence_locomotion_values(data):
    if np.sum(data.running) == 0 : 
        raise ValueError('No locomotion values present in the datafile')
    
def compute_arousal_mask(Episodes, speed_thr = 0.5) : 

    n_ep = len(Episodes.repeat)
    run_mask = np.zeros(n_ep, dtype=bool)
    rest_mask = np.zeros(n_ep, dtype=bool)
    temporal_mask = (0 <= Episodes.t) & (Episodes.t <= Episodes.time_duration[0]) #only consider locomotion during the visual stimuli

    for ep in range(n_ep):
        if np.mean(Episodes.running[ep][temporal_mask]) >= speed_thr :  
            run_mask[ep] = True
        else :
            rest_mask[ep] = True

    return run_mask, rest_mask 

def get_prefered_angles(prefered_angles_dataset, filename):

    short_filename = os.sep.join(filename.split('/')[-3:])
    prefered_angles = prefered_angles_dataset['prefered_angles'][(prefered_angles_dataset['short_filenames'] == short_filename)]

    if not prefered_angles.empty : 
        return prefered_angles.values[0]
    else : 
        raise ValueError('file %s not present in the prefered_angles_dataset' % short_filename)


def get_prefered_angles_dataset(data_folder, summary_folder, contrast) :

    folder = data_folder.split('/')[-2]
    prefered_angles_dataset = pd.DataFrame({})
    try : 
        Tunings = np.load(summary_folder + '/Tunings_%s_contrast-%s.npy' % (folder, contrast), allow_pickle=True)
        prefered_angles_dataset['short_filenames'] = [os.sep.join(Tuning['datafile'].split('/')[-3:]) for Tuning in Tunings] #ignore basepath
        prefered_angles_dataset['contrasts'] = [contrast]*len(prefered_angles_dataset['short_filenames'])
        prefered_angles_dataset['prefered_angles'] = [Tuning['prefered_angles'] for Tuning in Tunings]
        
    except: 
        print('No Tunings file found for folder %s and contrast %s' % (folder, contrast))
        print('You need to build the classical tunings summary first, before building one with a filtering condition')


    return prefered_angles_dataset

