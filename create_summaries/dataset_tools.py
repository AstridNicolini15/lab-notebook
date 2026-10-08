#%% 
import os, sys
import numpy as np
import pandas as pd 

os.chdir('/home/user/lab-notebook/astrid')
sys.path.insert(0,'/home/user/lab-notebook/astrid')
from general_tools.general_tools import *
sys.path += ['./physion/src']
import physion.utils.plot_tools as pt
from physion.analysis.read_NWB\
                         import scan_folder_for_NWBfiles, Data
from physion.analysis.episodes.build import EpisodeData

from physion.analysis.protocols.orientation_tuning import *
from physion.analysis.protocols.contrast_sensitivity import *

#%%

def get_summary_prefix_name(quantity, filtering_cond_name, dFoF_params, special_dFoF_params) :

    if quantity == 'dFoF' : 
        summary_prefix_name = ''
    elif quantity[:11] == 'Deconvolved' : 
        summary_prefix_name = quantity + '_'
    else : 
        raise ValueError("quantity should be either 'Deconvolved_something' or 'dFoF'")

    if filtering_cond_name is not None : 
        summary_prefix_name += filtering_cond_name  + '_' #Run or Rest for exemple

    if special_dFoF_params is not None: 
        summary_prefix_name += 'inclusion_factor-%s_correction_factor-%s_' % (dFoF_params['roi_to_neuropil_fluo_inclusion_factor'], dFoF_params['neuropil_correction_factor'])

    if dFoF_params['roi_specific_neuropil_correction_factor'] :
        summary_prefix_name += 'allen_neuropil_'


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
    temporal_mask = (0 <= Episodes.t) & (Episodes.t <= Episodes.time_duration[0]) #only consider locomotion during the visual stimuli

    run_mask = [np.mean(Episodes.running[ep][temporal_mask]) >= speed_thr for ep in range(n_ep)]
    run_mask = np.array(run_mask).astype(bool)
    rest_mask = ~run_mask

    return run_mask, rest_mask 


def get_perc_rest_from_rest_summary(rest_summary, return_mean = True) :

    n_trials = 15*8
    sessions_perc_rest = [np.sum(r_s['ntrials'])/n_trials for r_s in rest_summary]

    if return_mean : 
        return np.round(np.nanmean(sessions_perc_rest),3)
    else : 
        return sessions_perc_rest


def compute_perc_rest_from_nwb(folder, base_path, protocol = 'Tunings', compute = 'time', speed_thr = 0.5) :

    if protocol == 'Tunings' :
        control_cond = '8orientation'
    if protocol == 'Sensitivities' :
        control_cond = '8contrasts'

    DATASET = scan_folder_for_NWBfiles(base_path + '/' + folder + '/NWBs')
    filenames = [DATASET['files'][i] for i in range(len(DATASET['files'])) if control_cond in DATASET['protocols'][i][0]]

    perc_rest = []

    for filename in filenames :

        data = physion.analysis.read_NWB.Data(filename, verbose=False) 
        data.build_running()
        check_presence_locomotion_values(data)

        Episodes = EpisodeData(data, 
                                quantities=['running'], 
                                protocol_name = data.protocols[0],
                                verbose=False)

        if compute == 'time' :
            #percentage of time running during the entire session
            perc_rest.append(np.round((np.sum([data.running >= speed_thr]) * 100 )/ len(data.running),2))

        if compute == 'episodes' :
            #percentage of episodes considered runned
            perc_rest.append(np.round((np.sum(compute_arousal_mask(Episodes, speed_thr = speed_thr)[0]) *100) / len(Episodes.repeat),2))

    return perc_rest
