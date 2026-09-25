#%% 
import os, sys , shutil 
import numpy as np
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


#%%%%%%%%%%%%%%%old fcts, to integrate or deleted at some point%%%%%%%%%%%%%%%
def compute_perc_ep_run(folder, summary_path, return_mean = True) :
    
    PERC_EP = []
    keys =  ['%s_angle-0.0' % folder, 
                        '%s_angle-90.0' % folder]
    run0_Sensitivities = np.load(summary_path + '/Run_Sensitivities_%s.npy' % keys[0], 
            allow_pickle=True) 
    run90_Sensitivities = np.load(summary_path + '/Run_Sensitivities_%s.npy' % keys[1], 
        allow_pickle=True) 
    
    for j in range(len(run0_Sensitivities)) : 
        PERC_EP.append(((np.nansum(run0_Sensitivities[j]['ntrials']) + np.nansum(run90_Sensitivities[j]['ntrials'])) *100)/240)
    if return_mean == False : 
        return PERC_EP
    return np.round(np.nanmean(PERC_EP),2)



def compute_perc_time_run(base_path, folders, speed_thrs, return_mean = True):

    """
    Over do it a bit + necessitate the nwb files => if only perc_ep can be done more efficiently from Sensitivities
    If perc_ep the following code is necessary i believe 
    """

    raw_d = [[],[],[],[],[]]

    for folder in  folders : 

        DATASET = scan_folder_for_NWBfiles(base_path + '/' + folder + '/NWBs')
        files = [DATASET['files'][i] for i in range(len(DATASET['files'])) if '8contrasts' in DATASET['protocols'][i][0]]

        raw_d[0] = np.concatenate((raw_d[0],[folder]*len(files)*len(speed_thrs)))
        raw_d[1] = np.concatenate((raw_d[1],files*len(speed_thrs)))
        raw_d[2] = np.concatenate((raw_d[2],speed_thrs*len(files)))

        for file in files :

            data = physion.analysis.read_NWB.Data(file, verbose=False) 
            data.build_running()
            if np.sum(data.running) == 0 :
                raw_d[3].append(np.nan)
                raw_d[4].append(np.nan)

            else : 
                Episodes = EpisodeData(data, 
                                        quantities=['running'], 
                                        protocol_name = data.protocols[0],
                                        verbose=False)

                for speed_thr in speed_thrs :

                    #percentage of time running during the entire session
                    raw_d[3].append(np.round((np.sum([data.running >= speed_thr]) * 100 )/ len(data.running),2))

                    #percentage of episodes considered runned
                    raw_d[4].append(np.round((np.sum(compute_arousal_mask(Episodes, speed_thr = speed_thr)[0]) *100) / len(Episodes.repeat),2))



    d = pd.DataFrame(data = np.array(raw_d).T, columns= ['folder','file','speed_thr','perc_sess','perc_ep'])
    d['speed_thr'] = d['speed_thr'].astype(float)
    d['perc_sess'] = d['perc_sess'].astype(float)
    d['perc_ep'] = d['perc_ep'].astype(float)

    if return_mean == True : 

        raw_mean_d = [[],[],[],[]]

        for folder in folders :
            raw_mean_d[0] = np.concatenate((raw_mean_d[0],[folder]*len(speed_thrs)))
            raw_mean_d[1] = np.concatenate((raw_mean_d[1],speed_thrs))

            for speed_thr in speed_thrs :

                raw_mean_d[2].append(np.nanmean(d['perc_sess'].where((d['folder'] == folder) & (d['speed_thr'] == speed_thr)).values))
                raw_mean_d[3].append(np.nanmean(d['perc_ep'].where((d['folder'] == folder) & (d['speed_thr'] == speed_thr)).values))

        mean_d = pd.DataFrame(data = np.array(raw_mean_d).T, columns= ['folder','speed_thr','perc_sess','perc_ep'])
        mean_d['speed_thr'], mean_d['perc_sess'], mean_d['perc_ep'] = mean_d['speed_thr'].astype(float), mean_d['perc_sess'].astype(float), mean_d['perc_ep'].astype(float)
        #mean_d['perc_sess'] = mean_d['perc_sess'].astype(float)
        #mean_d['perc_ep'] = mean_d['perc_ep'].astype(float)
        return mean_d

    else : 
        return d 
