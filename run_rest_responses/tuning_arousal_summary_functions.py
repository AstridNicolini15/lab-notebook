#%%
import os, sys , shutil 
import numpy as np
from scipy import stats

os.chdir('/home/user/lab-notebook/astrid')
sys.path += ['./physion/src']
import physion.utils.plot_tools as pt
from physion.analysis.read_NWB\
                         import scan_folder_for_NWBfiles, Data
from physion.analysis.episodes.build import EpisodeData
from physion.analysis.protocols.orientation_tuning import *
from run_rest_responses.contrast_arousal_summary_functions import * 
#%%

def session_sem_with_indepedance_hypothesis_universal(Summaries, average_by = 'sessions') : #done very detailed because easily confusable


    #compute sessions means sems
    sessions_sems = []
    if 'shifted_angle' in Summaries[0].keys():
        visual_stim = 'shifted_angle'
    elif 'contrast' in Summaries[0].keys():
        visual_stim = 'contrast'

    for s in Summaries :

        #compute sem of each ROI mean. ie ROI mean std due to eps. ie propagated eps std / square number of eps
        ROIs_sems = s['std-values'] / np.sqrt(s['ntrials']) #does not perform the same way for both summaries bc in tuning ntrials is two dim while in sensitivity it is 1 dim 

        #compute between-ROIs-variance-due-to-eps. ie the contributions of each ROIs-variance-due-to-eps to the between-ROI variance. is = to the mean of the variances due to eps. it is NOT the second propagation. bc we are in between ROIs not after ROIs mean 
        ROIs_variance_due_to_eps = np.nanmean(ROIs_sems**2, axis = 0) #mean of squared stds ie mean of variances ie is a variance

        #compute observed between ROIs variance
        observed_ROIs_variance = np.nanvar(s['Responses'], axis=0, ddof = 1)

        #compare observed-between-ROIs-variance to the mean of ROIs-variances-due-to-eps
        session_variance = np.nanmax((observed_ROIs_variance, ROIs_variance_due_to_eps), axis = 0)  #check for the nans propagation law
        #compute sem of the session mean 
        n_rois_per_visual_stim = [np.sum(np.isnan(s['Responses'][:,c]) == False) for c in range(s[visual_stim].shape[0])] #to account for nans in the following divisions
        session_sem = np.sqrt(session_variance) / np.sqrt(n_rois_per_visual_stim) #sqrt of the std divided by sqrt number of ROIs ie sem of the session
        
        sessions_sems.append(session_sem)

    sessions_sems = np.array(sessions_sems)


#compute population mean sem 

    #compute observed-between-sessions-variance 
    if visual_stim == 'shifted_angle' : 
        Responses = np.array(get_tuning_responses(Summaries, average_by=average_by))
    elif visual_stim == 'contrast' : 
        Responses = np.array(get_responses(Summaries, average_by=average_by))

    observed_sessions_variance = np.nanvar(Responses, axis=0, ddof = 1)

    #compute between-sessions-variance-due-to-each-respective-sems ie the contributions of each session sem to the between-sessions variance. #is the mean of the previously computed squared sems (because sems is std and here we work with variance)
    sessions_variance_due_to_sems = np.nanmean(sessions_sems**2, axis = 0)

    #compare observed-between-sessions-variance to mean of sessions-variances computed before
    population_variance = np.nanmax((observed_sessions_variance, sessions_variance_due_to_sems), axis = 0)  #check for the nans propagation law

    #compute population sem
    n_sessions_per_visual_stim = [np.sum(np.isnan(Responses[:,c]) == False) for c in range(s[visual_stim].shape[0])] #to account for nans in the following divisions
    population_sem = np.sqrt(population_variance) / np.sqrt(n_sessions_per_visual_stim) #std divided by sqrt number of sessions ie sem of the session

    return population_sem
#%%


def get_tuning_responses(Tunings,
                         average_by='sessions'):
    """
    Change to include the normalisation by the value in prefered orientations
    """
    if average_by=='sessions':
        # mean significant responses per session
        Responses = [np.nanmean(Tuning['Responses'][Tuning['significant_ROIs'],:],
                        axis=0) for Tuning in Tunings]
        Responses = np.array([r/r[1] for r in Responses])

    elif average_by=='subjects':
        subjects = np.array([Tuning['subject']\
                                for Tuning in Tunings])
        Responses = []
        # mean significant responses per session
        for subj in np.unique(subjects):
            sCond = (subjects==subj)
            Responses.append(\
                np.mean(\
                    np.concatenate([\
                        Tunings[i]['Responses'][\
                            Tunings[i]['significant_ROIs'],:]\
                                 for i in np.arange(len(subjects))[sCond]]),
                    axis=0))
            
    elif average_by=='ROIs':
        # mean significant responses per session
        Responses = np.concatenate([\
                        Tuning['Responses'][Tuning['significant_ROIs'],:]\
                                                    for Tuning in Tunings])

    else:
        print()
        print(' choose average_by either "sessions", "subjects" or "ROIs"  ')
        print()

    return Responses

#%%

def compute_tuning_response_per_cells(data, Episodes,
                                      stat_test_props,
                                      response_significance_threshold = 0.05,
                                      filtering_cond=None,
                                      quantity='dFoF',
                                      contrast=1.0,
                                      nMin_episodes = 2,
                                      start_angle=-22.5, 
                                      angle_range=180,
                                      verbose=False):
    """

    All cells are considered in this analysis !!
      --> think about filtering them by resp['significant_ROIs'] when needed !!

    """

    shifted_angle = np.array(\
        [shift_orientation_according_to_pref(r, pref_angle=-start_angle,
                                             start_angle=start_angle,
                                             angle_range=angle_range)\
                    for r in Episodes.varied_parameters['angle']])

    if verbose:
        print('  - shifted_angle correspond to : ', shifted_angle)

    if filtering_cond is None:
        filtering_cond = Episodes.find_episode_cond() # True everywhere

    cond = Episodes.find_episode_cond(key='contrast', 
                                        value=contrast) &\
                                        filtering_cond
    
    summary = Episodes.pre_post_statistics(stat_test_props = stat_test_props,
                                                episode_cond=cond,
                                                repetition_keys=['repeat', 'contrast'],
                                                response_args=dict(quantity=quantity),
                                                response_significance_threshold=response_significance_threshold,
                                                multiple_comparison_correction=False,
                                                loop_over_cells=True,
                                                nMin_episodes = nMin_episodes,
                                                verbose=verbose)
        
    # if significant in at least one orientation
    significant = (np.sum(summary['significant'], axis=1)>0)

    # find preferred angle:
    ipref = np.argmax(summary['value'], axis=1).flatten()
    #print(ipref)

    prefered_angles = np.array(\
            [summary['angle'][i] for i in ipref])

    selectivities = np.array([\
        selectivity_index(summary['angle'],
                          summary['value'][roi, :])\
                            for roi in range(data.nROIs)])

    RESPONSES, semRESPONSES, Ntrials = [], [], []
    for roi in range(data.nROIs):

        RESPONSES.append(np.zeros(len(shifted_angle)))
        semRESPONSES.append(np.zeros(len(shifted_angle)))
        Ntrials.append(np.zeros(len(shifted_angle)))

        for angle, value, std, ntrials in zip(\
            summary['angle'],
            summary['value'][roi,:], 
            summary['std-value'][roi,:],
            summary['ntrials']):

            new_angle = shift_orientation_according_to_pref(angle,
                                                    pref_angle=prefered_angles[roi],
                                                    start_angle=start_angle,
                                                    angle_range=angle_range)
            iangle = np.flatnonzero(shifted_angle==new_angle)[0]

            RESPONSES[-1][iangle] = value
            semRESPONSES[-1][iangle] = std/np.sqrt(ntrials)
            Ntrials[-1][iangle] = ntrials


    return {'Responses':np.array(RESPONSES),
            'semResponses':np.array(semRESPONSES),
            'selectivities':np.array(selectivities),
            'shifted_angle':np.array(shifted_angle),
            'prefered_angles':np.array(prefered_angles),
            'significant_ROIs':np.array(significant),
            'std-values':summary['std-value'], 
            'ntrials': np.array(Ntrials[0])}

#--------------------------#

def compute_tuning_response_per_cells_with_arousal_cond(data, Episodes,
                                      arousal_cond,
                                      stat_test_props,
                                      response_significance_threshold = 0.05,
                                      filtering_cond=None,
                                      quantity='dFoF',
                                      contrast=1.0,
                                      nMin_episodes = 2,
                                      start_angle=-22.5, 
                                      angle_range=180,
                                      verbose=False) : 
    """
    /!\ Choice of significance based on summary without condition. 
    If statistical signficiance test were to be done on the actual run response => way less significant ROI ie cf commented lines 
    """
    Tuning_wo_cond = compute_tuning_response_per_cells(data, Episodes, 
                                                quantity='dFoF', 
                                                stat_test_props = stat_test_props, 
                                                response_significance_threshold = response_significance_threshold, 
                                                filtering_cond = None,
                                                contrast = contrast,
                                                nMin_episodes = nMin_episodes,
                                                start_angle=-22.5, 
                                                angle_range=180,
                                                verbose=False)


    if arousal_cond == '' : 
        Tuning = Tuning_wo_cond

    else : 
        Tuning  = {}
        for key in  Tuning_wo_cond.keys() : 
            if key == 'ntrials' : 
                temp_arr = np.zeros(Tuning_wo_cond['Responses'].shape)
            else : 
                temp_arr = np.zeros(Tuning_wo_cond[key].shape)
            temp_arr[:] = np.nan
            Tuning[key] = temp_arr
        
        Tuning['prefered_angles'] = Tuning_wo_cond['prefered_angles']
        Tuning['shifted_angle'] = Tuning_wo_cond['shifted_angle']
        Tuning['significant_ROIs'] = Tuning_wo_cond['significant_ROIs']


        if arousal_cond == 'Run_' : 
            filtering_arousal_cond = compute_arousal_mask(Episodes)[0]
        elif arousal_cond == 'Rest_' : 
            filtering_arousal_cond = compute_arousal_mask(Episodes)[1]


        cond = Episodes.find_episode_cond(key='contrast', 
                                            value=contrast) &\
                                            filtering_arousal_cond
        
        summary = Episodes.pre_post_statistics(stat_test_props = stat_test_props,
                                                    episode_cond=cond,
                                                    repetition_keys=['repeat', 'contrast'],
                                                    nMin_episodes = nMin_episodes,
                                                    response_args=dict(quantity=quantity),
                                                    response_significance_threshold=response_significance_threshold,
                                                    multiple_comparison_correction=False,
                                                    loop_over_cells=True,
                                                    verbose=verbose)
        #significant = np.zeros(data.nROIs)
        for i,angle in enumerate(summary['angle']) : 
            for roi in range(data.nROIs) : 
                new_angle = shift_orientation_according_to_pref(angle,
                                                            pref_angle=Tuning['prefered_angles'][roi],
                                                            start_angle=start_angle,
                                                            angle_range=angle_range)
                iangle = np.flatnonzero(Tuning['shifted_angle']==new_angle)[0]

                Tuning['Responses'][roi][iangle] = summary['value'][roi,i]
                Tuning['std-values'][roi][iangle] = summary['std-value'][roi,i]
                Tuning['ntrials'][roi][iangle] = summary['ntrials'][i]
                #significant[roi] += np.sum(summary['significant'][roi])

            #new_angle = shift_orientation_according_to_pref(summary['angle'][i],
            #                                        pref_angle=0,
            #                                        start_angle=start_angle,
            #                                        angle_range=angle_range)
            #iangle = np.flatnonzero(Tuning['shifted_angle']==new_angle)[0]

            #Tuning['ntrials'][iangle] = summary['ntrials'][i]
        
        #Tuning['significant_ROIs'] = significant > 0 
    return Tuning 

#####PLOT######


def plot_orientation_tuning_curve_with_uncertainty(keys,
                              arousal_keys,
                              summary_path='',
                              average_by='sessions',
                              uncertainties = ['std'],
                              colors = None, 
                              plot_perc_run = True,
                              group_ROIs = False,
                              gaussian_fit = True,
                              base_path = '',
                              ylims=[[0,1]]) : 
        
    if colors is None:
        colors = pt.plt.rcParams['axes.prop_cycle'].by_key()['color']
        
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
            perc_ep = compute_perc_ep_run(folder)
            fig.text(x=0, y=-0.2, s='%' + ' of episodes considered runned : ' + str(perc_ep))

        fig.text(x=0, y=1.1, s=folder)
        fig.text(x=0, y=-0.1, s= 'contrast 1.0', color = colors[0])
        fig.text(x=0.12, y=-0.1, s= 'contrast 0.5', color = colors[1])
        fig.text(x=0, y=-0.3, s= 'error bars = ' + uncertainty)


        x = np.linspace(-30, 180-30, 100)
        for k, arousal_cond in enumerate(arousal_keys): 

            for key, color in zip(keys, colors):

                # load data
                Tunings = np.load(summary_path + '/' + arousal_cond + 'Tunings_%s.npy' % key, allow_pickle=True)   

                Responses = get_tuning_responses(Tunings, average_by=average_by) #is normed
                # Gaussian Fit
                C, func = fit_gaussian(Tunings[0]['shifted_angle'],
                                        np.nanmean(Responses, axis=0))

                axes[k].plot(x, func(x), lw=2, alpha=.5, color=color)

                if uncertainty == 'std' : 
                    uncertainty_sy = np.nanstd(Responses, axis=0, ddof = 1)

                elif uncertainty == 'sem' : 
                    uncertainty_sy = stats.sem(Responses, axis=0, nan_policy = 'omit', ddof = 1) 

                elif uncertainty == 'session sem with propagation and independance hypothesis' : 
                    uncertainty_sy = session_sem_with_indepedance_hypothesis_universal(Tunings)

                pt.scatter(Tunings[0]['shifted_angle'], np.nanmean(Responses, axis=0), 
                sy=uncertainty_sy, 
                color=color, ax=axes[k], ms=2)
        

            pt.set_plot(axes[k], xticks=Tunings[0]['shifted_angle'], 
                        #yticks=np.arange(3)*0.5, 
                        ylim=ylims,
                        ylabel='norm. $\\delta$ $\\Delta$F/F',  
                        xlabel='angle ($^o$) from pref.',
                        title=titles[k], 
                        xticks_labels=['%i' % a if (a in [0, 90]) else '' for a in Tunings[0]['shifted_angle'] ])

    return fig, axes


def plot_orientation_tuning_curve_with_uncertainty_all_pop_in_one_graph(folders,
                              arousal_keys,
                              summary_path='',
                              average_by='sessions',
                              uncertainties = ['std'],
                              COLORS = None, 
                              plot_perc_run = True,
                              group_ROIs = False,
                              gaussian_fit = True,
                              base_path = '',
                              YLIMS=[[0,1]]) : 
    

    fig_args={'right':25, 'ax_scale':(1.2, 1.7)}
            
    fig, ax = pt.figure(**fig_args)
    inset = pt.inset(ax, [2.5,0,1,1])
    inset2 = pt.inset(ax, [4.5,0,1,1])

    axes = [ax, inset, inset2]
    titles = [x[:-1] for x in arousal_keys]

    for i,folder in enumerate(folders) : 
        print(i)
        keys =  ['%s_contrast-1.0' % folder, 
            '%s_contrast-0.5' % folder]
        ylims = YLIMS[i]
        colors = COLORS[i]
        folder = keys[0][:-13]

        for uncertainty in uncertainties : 

            if type(keys)==str:
                keys, colors = [keys], [colors[0]]

            if plot_perc_run == True : 
                perc_ep = compute_perc_ep_run(folder, summary_path)
                fig.text(x=0, y=-0.2, s='%' + ' of episodes considered runned : ' + str(perc_ep))


            fig.text(x=0, y=1.1-(0.1*i), s=folder, color = colors[0])
            fig.text(x=0+(0.3*i), y=-0.1, s= 'contrast 1.0', color = colors[0])
            fig.text(x=0.15+(0.3*i), y=-0.1, s= 'contrast 0.5', color = colors[1])
            fig.text(x=0, y=-0.2, s= 'error bars = ' + uncertainty)



            x = np.linspace(-30, 180-30, 100)
            for k, arousal_cond in enumerate(arousal_keys): 

                for key, color in zip(keys, colors):

                    # load data
                    Tunings = np.load(summary_path + '/' + arousal_cond + 'Tunings_%s.npy' % key, allow_pickle=True)   

                    Responses = get_tuning_responses(Tunings, average_by=average_by) #is normed
                    Responses = np.array([r/r[1] for r in Responses])
                    # Gaussian Fit
                    C, func = fit_gaussian(Tunings[0]['shifted_angle'],
                                            np.nanmean(Responses, axis=0))

                    axes[k].plot(x, func(x), lw=2, alpha=.5, color=color)

                    if uncertainty == 'std' : 
                        uncertainty_sy = np.nanstd(Responses, axis=0, ddof = 1)

                    elif uncertainty == 'sem' : 
                        uncertainty_sy = stats.sem(Responses, axis=0, nan_policy = 'omit', ddof = 1) 

                    elif uncertainty == 'session sem with propagation and independance hypothesis' : 
                        uncertainty_sy = session_sem_with_indepedance_hypothesis_universal(Tunings)

                    pt.scatter(Tunings[0]['shifted_angle'], np.nanmean(Responses, axis=0), 
                    sy=uncertainty_sy, 
                    color=color, ax=axes[k], ms=2)
            

                pt.set_plot(axes[k], xticks=Tunings[0]['shifted_angle'], 
                            #yticks=np.arange(3)*0.5, 
                            ylim=ylims,
                            ylabel='norm. $\\delta$ $\\Delta$F/F',  
                            xlabel='angle ($^o$) from pref.',
                            title=titles[k], 
                            xticks_labels=['%i' % a if (a in [0, 90]) else '' for a in Tunings[0]['shifted_angle'] ])

    return fig, axes