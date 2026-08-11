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

def compute_sensitivity_per_cells(data, Episodes, 
                                  stat_test_props,
                                  response_significance_threshold = 0.05,
                                  filtering_cond=None,
                                  quantity='dFoF',
                                  angle=0.0,
                                  nMin_episodes = 2,
                                  verbose=False):
    """
    Modified only to pass : nMin_episodes arguments and add keys : std-values, ntrials
    note : 
    nMin_episodes is reduced to 2, because Run conditions have only a few if not none episodes for every contrast.
    """

    if filtering_cond is None:
        filtering_cond = Episodes.find_episode_cond() # True everywhere

    cond = Episodes.find_episode_cond(key='angle', 
                                        value=angle) &\
                                        filtering_cond
    
    # first for positive responses
    stat_test_props['sign'] = 'positive'
    summary_positive= Episodes.pre_post_statistics(\
                                stat_test_props = stat_test_props,
                                episode_cond=cond,
                                repetition_keys=['repeat', 'angle'],
                                nMin_episodes = nMin_episodes,
                                response_args=dict(quantity=quantity),
                                response_significance_threshold=response_significance_threshold,
                                multiple_comparison_correction=False,
                                loop_over_cells=True,
                                verbose=verbose)
    
    # second for negative responses
    stat_test_props['sign'] = 'negative'
    summary_negative = Episodes.pre_post_statistics(\
                            stat_test_props = stat_test_props,
                            episode_cond=cond,
                            repetition_keys=['repeat', 'angle'],
                            nMin_episodes = nMin_episodes,
                            response_args=dict(quantity=quantity),
                            response_significance_threshold=response_significance_threshold,
                            multiple_comparison_correction=False,
                            loop_over_cells=True,
                            verbose=verbose)

    semRESPONSES = np.array([
        summary_positive['value'][roi,:]/np.sqrt(summary_positive['ntrials'])
        for roi in range(data.nROIs)])

    output = {'Responses':summary_positive['value'],
              'semResponses':semRESPONSES,
              'contrast':summary_positive['contrast'],
              'significant_pos':summary_positive['significant'],
              'significant_neg':summary_negative['significant'],
              'std-values':summary_positive['std-value'], 
              'ntrials':summary_positive['ntrials']}

    return output

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


def correct_missing_responses_and_stds_with_nans(Sensitivity, Episodes): 
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

def correct_missing_responses_and_stds_with_nans_universal(summary, Episodes): 

    corrected_summary = {}

    if '8orientation' in Episodes.protocol_name : 
        stimulus = 'shifted_angle'
        stimulus_values = summary[stimulus] #in Tuning all values exists 

    elif '8contrasts' in Episodes.protocol_name : 
        stimulus = 'contrast'
        stimulus_values = np.unique(Episodes.contrast)

    nROIs = Episodes.data.nROIs

    for key in summary.keys() : 
        corrected_summary[key] = []

    for stim_val in stimulus_values : 
        print(stim_val)
        ind_stim = np.argwhere(summary[stimulus] == stim_val)[0][0]
        if (stimulus == 'shifted_angle' and np.sum(summary['Responses'][:,ind_stim]) != 0) or (stimulus == 'contrast' and stim_val in summary[stimulus]) : 
            print('a')
            for key in summary.keys() :
                print(key)
                if len(summary[key].shape) == 2: 
                    corrected_summary[key].append(summary[key][:,ind_stim])
                elif len(summary[key].shape) == 1 :
                    corrected_summary[key].append(summary[key][ind_stim])
        else : 
            print('b')
            for key in summary.keys() : 

                if key in ['significant_pos','significant_neg'] : 
                    corrected_summary[key].append([False]*nROIs)
                elif key == 'significant_ROIs' : 
                    corrected_summary[key].append([False])
                elif key == 'ntrials' : 
                    corrected_summary[key].append(np.nan)
                elif key != stimulus: 
                    corrected_summary[key].append([np.nan]*nROIs)

    #reshape correctly 
    for key in summary.keys() : 
        if key != 'contrast' : 
            corrected_summary[key]  = np.array(corrected_summary[key]).T
    corrected_summary['contrast'] = stimulus_values


    return corrected_summary


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


def session_sem_with_indepedance_hypothesis(Sensitivities) : #done very detailed because easily confusable


    #compute sessions means sems
    sessions_sems = []

    for s in Sensitivities :
    
        #compute sem of each ROI mean. ie ROI mean std due to eps. ie propagated eps std / square number of eps
        propag1_ROIs_sems = s['std-values'] / np.sqrt(s['ntrials'])

        #compute between-ROIs-variance-due-to-eps. ie the contributions of each ROIs-variance-due-to-eps to the between-ROI variance. is = to the mean of the variances due to eps. it is NOT the second propagation. bc we are in between ROIs not after ROIs mean 
        ROIs_variance_due_to_eps = np.nanmean(propag1_ROIs_sems**2, axis = 0) #mean of squared stds ie mean of variances ie is a variance

        #compute observed between ROIs variance
        observed_ROIs_variance = np.nanvar(s['Responses'], axis=0, ddof = 1)

        #compare observed-between-ROIs-variance to the mean of ROIs-variances-due-to-eps
        session_variance = np.nanmax((observed_ROIs_variance, ROIs_variance_due_to_eps), axis = 0)  #check for the nans propagation law
        #compute sem of the session mean 
        n_rois_per_contrast = [np.sum(np.isnan(s['Responses'][:,c]) == False) for c in range(s['contrast'].shape[0])] #to account for nans in the following divisions
        session_sem = np.sqrt(session_variance) / np.sqrt(n_rois_per_contrast) #std divided by sqrt number of ROIs ie sem of the session

        #more detailed version that is equivalent to the one above. 
        #session_sem = []
        #for i in range(len(s['contrast'])) : 
        #    if observed_ROIs_variance[i] > ROIs_variance_due_to_eps[i] :
        #        sem_i =  np.sqrt(observed_ROIs_variance[i]) / np.sqrt(n_rois_per_contrast[i]) 
        #        session_sem.append(sem_i)
        #    else :
        #        sem_i = np.sqrt(np.nansum(propag1_ROIs_sems[:,i]) / (n_rois_per_contrast[i]**2)) #here the second propagation
        #        session_sem.append(sem_i)
        
        sessions_sems.append(session_sem)

    sessions_sems = np.array(sessions_sems)


#compute population mean sem 

    #compute observed-between-sessions-variance 
    Responses = np.array(get_responses(Sensitivities, average_by='sessions'))
    observed_sessions_variance = np.nanvar(Responses, axis=0, ddof = 1)

    #compute between-sessions-variance-due-to-each-respective-sems ie the contributions of each session sem to the between-sessions variance. #is the mean of the previously computed squared sems (because sems is std and here we work with variance)
    sessions_variance_due_to_sems = np.nanmean(sessions_sems**2, axis = 0)

    #compare observed-between-sessions-variance to mean of sessions-variances computed before
    population_variance = np.nanmax((observed_sessions_variance, sessions_variance_due_to_sems), axis = 0)  #check for the nans propagation law

    #compute population sem
    n_sessions_per_contrast = [np.sum(np.isnan(Responses[:,c]) == False) for c in range(s['contrast'].shape[0])] #to account for nans in the following divisions
    population_sem = np.sqrt(population_variance) / np.sqrt(n_sessions_per_contrast) #std divided by sqrt number of sessions ie sem of the session

    return population_sem



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
 
 
def plot_contrast_sensitivity_with_uncertainty_all_pop_in_one_graph(folders,
                              arousal_keys,
                              summary_path='',
                              average_by='sessions',
                              uncertainties = ['std','session sem with propagation and independance hypothesis' ],
                              COLORS = None, 
                              plot_perc_run = False,
                              group_ROIs = False,
                              base_path = '',
                              YLIMS = [-0.05,1] ) : 

    fig_args={'right':25, 'ax_scale':(1.2, 1.7)}
            
    fig, ax = pt.figure(**fig_args)
    inset = pt.inset(ax, [2.5,0,1,1])
    inset2 = pt.inset(ax, [4.5,0,1,1])

    axes = [ax, inset, inset2]
    titles = [x[:-1] for x in arousal_keys]

    for i,folder in enumerate(folders) : 
        print(i)
        keys =  ['%s_angle-0.0' % folder, 
            '%s_angle-90.0' % folder]
        ylims = YLIMS[i]
        colors = COLORS[i]
        folder = keys[0][:-10]
        for uncertainty in uncertainties : 

            if type(keys)==str:
                keys, colors = [keys], [colors[0]]

            if plot_perc_run == True : 
                perc_ep = compute_perc_ep_run(folder, summary_path)
                fig.text(x=0, y=-0.2-(0.1*i), s='%' + ' of episodes considered runned : ')
                fig.text(x=0.3, y=-0.2-(0.1*i), s=  str(perc_ep), color = colors[0])

            fig.text(x=0, y=1.1-(0.1*i), s=folder, color = colors[0])
            fig.text(x=0+(0.2*i), y=-0.1, s= 'angle 0.0', color = colors[0])
            fig.text(x=0.1+(0.2*i), y=-0.1, s= 'angle 90.0', color = colors[1])
            fig.text(x=0, y=-0.4, s= 'error bars = ' + uncertainty)


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
if False : 
    
    #get percentage running 
    print('Computing percentage of time running...')
    mean_d = compute_perc_time_run(base_path, folders, [0.1,0.5], return_mean = True)

    perc_sess = [x for x in mean_d['perc_sess'].where((mean_d['folder'] == folder) & (mean_d['speed_thr'] == speed_thr)).values if np.isnan(x) == False][0]
    perc_ep = [x for x in mean_d['perc_ep'].where((mean_d['folder'] == folder) & (mean_d['speed_thr'] == speed_thr)).values if np.isnan(x) == False][0]

#%% old functions, are not the correct way to propagate uncertainty, i keep bc idk
if False : 

    def propagate_stds(Sensitivities) : 

        sessions_stdevs = [np.nanstd(s['Responses'], axis = 0, ddof = 1) for s in Sensitivities]
        #propagated_stds = np.nansum(sessions_stdevs, axis = 0) / len(sessions_stdevs)
        propagated_stds = np.nanmean(sessions_stdevs, axis =0) #denominator not always = to len sessions bc of nans

        return propagated_stds


    def propagate_sems(Sensitivities) : 

        sessions_sems = [stats.sem(s['Responses'], axis = 0, nan_policy = 'omit', ddof = 0) for s in Sensitivities]
        propagated_sems = np.nanmean(sessions_sems, axis = 0) 

        return propagated_sems

    def double_propagate_stds(Sensitivities) : 

        propagated_eps_std = np.array([np.nanmean(s['std-values'], axis = 0) for s in Sensitivities]) 
        ROIs_std = np.array([np.nanstd(s['Responses'], axis = 0, ddof = 1) for s in Sensitivities])
        #for every session, for every contrast take the max between the two values
        max_ROIs_stds = np.array([[np.max((ROIs_std[i][j], propagated_eps_std[i][j])) for j in range(ROIs_std.shape[1])] for i in range(ROIs_std.shape[0])])

        propagated_ROIs_stds = np.nanmean(max_ROIs_stds, axis = 0)
        sessions_stds = np.nanstd(np.array([np.nanmean(s['Responses'], axis = 0) for s in Sensitivities]), axis = 0, ddof = 1)
        #for every contrast take the max between sessions stds or the propagated ROI stds 
        max_sessions_stds = [np.max((sessions_stds[i], propagated_ROIs_stds[i])) for i in range(len(sessions_stds))]

        return max_sessions_stds

    def double_propagate_sems(Sensitivities) : 

        propagated_eps_sem = np.array([np.nanmean(s['std-values'], axis = 0)/np.sqrt(s['ntrials']) for s in Sensitivities]) #divide by sqrt of episodes numbers
        ROIs_sem = np.array([np.nanstd(s['Responses'], axis = 0, ddof = 1)/np.sqrt(np.sum(np.isnan(s['Responses'])==False, axis =0)) for s in Sensitivities]) #divide by sqrt of ROIs numbers + complicated formulation to take nans into account
        #for every session, for every contrast take the max between the two values
        max_ROIs_sems = np.array([[np.max((ROIs_sem[i][j], propagated_eps_sem[i][j])) for j in range(ROIs_sem.shape[1])] for i in range(ROIs_sem.shape[0])])

        propagated_ROIs_sems = np.nanmean(max_ROIs_sems, axis = 0)/np.sqrt(np.sum(np.isnan(max_ROIs_sems)==False, axis = 0)) #divide by sqrt of number of sessions + take nans into account 
        sessions_responses = np.array([np.nanmean(s['Responses'], axis = 0) for s in Sensitivities])
        sessions_sems = np.nanstd(sessions_responses, axis = 0, ddof = 1) / np.sqrt(np.sum(np.isnan(sessions_responses)==False, axis = 0))#divide by sqrt of number of sessions that are not nans 
        #for every contrast take the max between sessions sems or the propagated ROI sems 
        max_sessions_sems = [np.max((sessions_sems[i], propagated_ROIs_sems[i])) for i in range(len(sessions_sems))]

        return max_sessions_sems

