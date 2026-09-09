#%%
import os, sys , shutil 
import multiprocessing
import numpy as np
from scipy import stats
import math 
import pandas as pd 

sys.path += ['./physion/src']
from physion.utils import plot_tools as pt
from physion.analysis.read_NWB\
                         import scan_folder_for_NWBfiles, Data
from physion.analysis.episodes.build import EpisodeData
from physion.analysis.protocols.contrast_sensitivity\
                        import compute_sensitivity_per_cells


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






#%%
#------------Main Plot functions----------------#


def plot_contrast_sensitivity_with_uncertainty(keys,
                              arousal_keys,
                              summary_path='',
                              average_by='sessions',
                              uncertainties = ['std','session sem with propagation and independance hypothesis' ],
                              colors = None, 
                              plot_perc_run = False,
                              group_ROIs = False,
                              base_path = '',
                              ylims = [-0.05,1] ) : 

    folder = keys[0][:-10]
    for uncertainty in uncertainties : 

        fig_args={'right':25, 'ax_scale':(1.2, 1.7)}
        
        fig, ax = pt.figure(**fig_args)
        inset = pt.inset(ax, [2.5,0,1,1])
        inset2 = pt.inset(ax, [4.5,0,1,1])

        axes = [ax, inset, inset2]
        titles = [x[:-1] for x in arousal_keys]

        if type(keys)==str:
            keys, colors = [keys], [colors[0]]

        if plot_perc_run == True : 
            perc_ep = compute_perc_ep_run(folder, summary_path)
            fig.text(x=0, y=-0.2, s='%' + ' of episodes considered runned : ' + str(perc_ep))

        fig.text(x=0, y=1.1, s=folder)
        fig.text(x=0, y=-0.1, s= 'angle 0.0', color = colors[0])
        fig.text(x=0.1, y=-0.1, s= 'angle 90.0', color = colors[1])
        fig.text(x=0, y=-0.3, s= 'error bars = ' + uncertainty)


        for k, arousal_cond in enumerate(arousal_keys): 
            
            for (key, color) in zip(keys, colors):
                Sensitivities = \
                    np.load(summary_path + '/' + arousal_cond + 'Sensitivities_%s.npy' % key, allow_pickle=True)   

                if not group_ROIs :
                    Responses = get_responses(Sensitivities, average_by=average_by)
                else : 
                    Responses =  [s['Responses'][j] for s in Sensitivities for j in range(len(s['Responses']))]
                    fig.text(x=0.5, y = -0.3, s= 'group by ROIs : ' + str(group_ROIs))

                    
                if uncertainty == 'std' : 
                    uncertainty_sy = np.nanstd(Responses, axis=0, ddof = 1)

                elif uncertainty == 'sem' : 
                    uncertainty_sy = stats.sem(Responses, axis=0, nan_policy = 'omit', ddof = 1)

                elif uncertainty == 'session sem with propagation and independance hypothesis' : 
                    uncertainty_sy = session_sem_with_indepedance_hypothesis(Sensitivities)

                pt.plot(Sensitivities[0]['contrast'], 
                        np.nanmean(Responses, axis=0), 
                        sy=uncertainty_sy,
                        color=color,
                        ax=axes[k])
            
            pt.set_plot(axes[k], 
                title=titles[k],     
                ylabel='$\\delta$ $\\Delta$F/F',  
                xlabel='contrast',
                xticks=np.arange(3)*0.5)
            
        ax.set_ylim(ylims)
        inset.set_ylim(ylims)
        inset2.set_ylim(ylims)

    return fig, axes


#%%
#------------Additional Plot functions----------------#

#from physion.utils import plot_tools as pt
#from scipy import stats
#import math 
#import pandas as pd 

def get_responses(Sensitivities,
                  average_by='sessions'):

    if average_by=='sessions':
        # mean significant responses per session
        Responses = [np.mean(S['Responses'], axis=0) for S in Sensitivities]

    elif average_by=='subjects':
        subjects = np.array([Sensitivitie['subject']\
                                for Sensitivitie in Sensitivities])
        Responses = []
        # mean significant responses per session
        for subj in np.unique(subjects):
            sCond = (subjects==subj)
            Responses.append(\
                np.mean(\
                    np.concatenate([\
                        Sensitivities[i]['Responses']\
                          for i in np.arange(len(subjects))[sCond]]),
                    axis=0))

    elif average_by=='ROIs':
        # mean significant responses per session
        Responses = np.concatenate([\
                        S['Responses'] for S in Sensitivities])

    else:
        print()
        print(' choose average_by either "sessions" or "ROIs"  ')
        print()

    return Responses

def get_gains(Responses, contrast):
        """ gain from linear fit"""
        return np.array([np.polyfit(contrast, r, 1)[0]\
                        for r in Responses])

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
#%%
 
def compute_group_contrast(Episodes, grouped_contrast_values = [0.05,1]) : 
    #artificially modify episodes contrast refs as if there were only 2 contrasts
    #not very pretty but work.. 

    contrast_values = np.unique(Episodes.contrast)
    separator = contrast_values[4]

    for i in range(len(Episodes.contrast)) :


        if Episodes.contrast[i] < separator : 
            Episodes.contrast[i] = grouped_contrast_values[0]
        else : 
            Episodes.contrast[i] = grouped_contrast_values[1]

    return None
