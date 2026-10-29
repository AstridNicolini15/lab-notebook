#%%
import os, sys 
os.chdir('/home/user/lab-notebook/astrid')
import numpy as np 
from scipy.stats import sem

from general_tools.propagated_uncertainty import *
from general_tools.general_tools import *
from create_summaries.Preprocessing_Settings import get_dFoF_params
sys.path.insert(0,'/home/user/lab-notebook/astrid/physion/src')

import physion.utils.plot_tools as pt
from physion.analysis.episodes.build import EpisodeData
from physion.analysis.protocols.orientation_tuning import get_tuning_responses, fit_gaussian 
from physion.analysis.protocols.contrast_sensitivity import *

import physion

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import to_rgb


#%%
#---------ORIENTATION TUNING AND CONTRAST SENSIVITY RESPONSES----------#

def draw_tuning_curve(Tunings, 
                    ax, 
                    color = None,
                    alpha = 1,
                    xy = (0,0),
                    uncertainty_sy = 'std',
                    graph_width_dict = {'lw' : 2, 'ms' : 3, 'fontsize' : 7},
                    clip_negative = False) : 
    

    Responses = get_tuning_responses(Tunings, average_by='sessions') #return only significant resp
    Responses = [r/r[1] for r in Responses]

    if clip_negative : 
        Responses = np.array([[max(x, 0) for x in r] for r in Responses])

    # Gaussian Fit
    C, func = fit_gaussian(Tunings[0]['shifted_angle'],
                            np.nanmean(Responses, axis=0))
    
    x_angles = Tunings[0]['shifted_angle']
    n_cells = np.sum([np.sum(t['significant_ROIs']) for t in Tunings])

    #draw
    lw = graph_width_dict['lw']
    ms = graph_width_dict['ms']
    fontsize = graph_width_dict['fontsize']

    x = np.linspace(-30, 180-30, 100)

    ax.scatter(x_angles, np.mean(Responses, axis=0), color=color, alpha = alpha)

    ax.plot(x, func(x), lw=lw, color=color, alpha = alpha)

    ax.annotate(text = 'N= ' + str(n_cells), xy = xy, color = color, fontsize = fontsize, alpha = alpha)

    if uncertainty_sy == 'std' : 
        uncertainty_sy = np.std(Responses, axis = 0, nan_policy = 'omit')
    elif uncertainty_sy == 'propagation' : 
        print('not implemented')
        #uncertainty_sy = session_sem_with_indepedance_hypothesis_universal(Tunings)

    ax.errorbar(Tunings[0]['shifted_angle'], 
            np.nanmean(Responses, axis=0),
            yerr=uncertainty_sy,
            alpha = alpha,
            elinewidth = 1,
            fmt = '.',
            color=color, 
            ms = ms)



def draw_sensitivitie_curve(Sensitivities,
                            ax, 
                            color = None,
                            alpha = 1,
                            xy = (0,0),
                            uncertainty_sy = 'std',
                            graph_width_dict = {'lw' : 2, 'ms' : 3, 'fontsize' : 4}
                            ):


    #get values
    #Responses = get_responses(Sensitivities)
    Responses = [np.mean(S['Responses'][np.sum(S['significant_pos'] + S['significant_neg'], axis = 1).astype(bool)], axis=0)\
                                for S in Sensitivities] #get significant responses

    n_cells = np.sum([np.sum(np.sum(S['significant_pos'] + S['significant_neg'], axis = 1).astype(bool)) for S in Sensitivities]) #TO CHECKs

    if uncertainty_sy == 'propagated' : 
        uncertainty_sy = session_sem_with_indepedance_hypothesis_universal(Sensitivities)
    elif uncertainty_sy == 'std' : 
        uncertainty_sy = np.std(Responses, axis = 0, nan_policy = 'omit')

    #draw
    lw = graph_width_dict['lw']
    ms = graph_width_dict['ms']
    fontsize = graph_width_dict['fontsize']

    pt.plot(Sensitivities[0]['contrast'], 
        np.nanmean(Responses, axis=0), 
        sy=uncertainty_sy,
        color=color,
        alpha = alpha,
        ax=ax)
                
    ax.annotate(text = 'N= ' + str(n_cells), xy = xy, color = color, alpha = alpha, fontsize = fontsize)

    
def plot_responses(ax,
                summary_filenames,
                colors,  
                alphas = None,
                summary_protocol = 'Tunings',
                clip_negative = True,
                special_dict = {}, 
                sensitivitie_uncertainty = 'propagated') :


    title, ylims = extract_special_dict(special_dict)
    if alphas == None : 
        alphas = [None for i in range(len(summary_filenames))]
    

    for i, summary_filename in enumerate(summary_filenames) : 
        
        keys = get_protocol_keys(summary_protocol = summary_protocol)

        for j, key in enumerate(keys) : 
            summaries = np.load(summary_filename + '%s.npy' % key, allow_pickle=True) 
            
            if summary_protocol == 'Tunings' : 

                xy =  (100+(30*j),1-(0.03*i))
                draw_tuning_curve(summaries,
                        ax = ax, 
                        color = colors[i][j],
                        alpha = alphas[i],
                        xy = xy, 
                        graph_width_dict = {'lw' : 3, 'ms' : 5, 'fontsize' : 10},
                        clip_negative = clip_negative) 


            elif summary_protocol == 'Sensitivities' : 

                xy = (0.05+(0.1*j),0.65-(0.04*i))
                draw_sensitivitie_curve(summaries,
                        ax = ax, 
                        color = colors[i][j],
                        alpha =  alphas[i],
                        xy = xy,
                        uncertainty_sy = sensitivitie_uncertainty,
                        graph_width_dict = {'lw' : 2, 'ms' : 3, 'fontsize' : 8})

    add_axis_labels_to_plot(ax, plot_type = summary_protocol)

    ax.set_title(title, fontsize = 13)
    ax.set_ylim(ylims)

#%%
#---------RESPONSIVENESS TO VISUAL STIM----------#
def get_responsiveness_to_visual_stim(summary, total_on):

    """
    The % of significant ROIs can be calculated on the total number of valid ROIs (total_on = 'nROIs_final')
    or on the total number of original ROIs (total_on = 'nROIs_original')
    """

    if 'significant_ROIs' in summary[0].keys(): 
        n_resp = np.sum([np.sum(Tuning['significant_ROIs']) for Tuning in summary])
        n_nonresp = np.sum([np.sum(~Tuning['significant_ROIs']) for Tuning in summary])

    elif 'significant_pos' in summary[0].keys(): 
        print('warning : positive and negative significants are considered')
        cells_responsiveness_by_contrasts = np.concatenate([np.sum(S['significant_pos'] + S['significant_neg'], axis = 1) for S in summary])

        n_resp = len([x for x in cells_responsiveness_by_contrasts if x!=0])
        n_nonresp = len([x for x in cells_responsiveness_by_contrasts if x==0])

    n_total = np.sum([s[total_on] for s in summary])

    perc_resp = (n_resp *100) / n_total
    perc_nonresp = (n_nonresp *100) / n_total

    return perc_resp, perc_nonresp


def draw_responsiveness_pie(ax, summary, color, total_on = 'nROIs_final') : 

    perc_resp, perc_nonresp = get_responsiveness_to_visual_stim(summary, total_on = total_on)
    ax.pie([perc_resp, perc_nonresp], 
                    colors = [color, 'white'],
                    startangle = 90,
                    wedgeprops={"edgecolor":"black",'linewidth': 1, 'width' : 0.6},
                    textprops = {"fontsize" : 10})
    
    ax.annotate(text = f'{perc_resp:.1f}% \n resp', xy= (-1.2,1), color = color, fontsize = 13)
    if total_on == 'nROIs_original' : 
        print('WARNING BUG : percentages are correct but pies are not coherent. idk why')
        ax.annotate(text = 'total on on nROIs_original',  xy= (0,0.7), fontsize = 13) #position to correct



def plot_responsiveness_pies(axes, summary_filenames, colors, total_on = 'nROIs_final', special_dict = None, summary_protocol = 'Tunings'):

    k = 0
    for i, summary_filename in enumerate(summary_filenames) : 

        keys = get_protocol_keys(summary_protocol = summary_protocol)
        
        for key in keys:

            summary = np.load(summary_filename + '%s.npy' % key, allow_pickle=True) 
            draw_responsiveness_pie(axes[k], summary, colors[i], total_on = total_on)
            k+=1

        axes[k-2].set_title('full contrast', color = 'black', loc = 'left', fontsize = 10, pad = 35)
        axes[k-1].set_title('half contrast', color = 'black', loc = 'left', fontsize = 10, pad = 35)



#%%
#---------MEAN FLUORESCENCE LEVEL----------#

def get_fluorescence_level(nwb_filenames, dFoF_parameters, F_value = 'correctedFluo0') : 

    '''
    F_value can be : 'dFoF', 'correctedFluo0','F0' or 'rawFluo'
    '''

    print("computing " + F_value)

    Fval = []

    for filename in nwb_filenames : 

        data = physion.analysis.read_NWB.Data(filename, verbose=False) 
        #getattr(data, 'build_' + F_value)
        data.build_dFoF(**dFoF_parameters, verbose=False)

        Fval.append(np.mean(getattr(data, F_value), axis = 1)) #one mean value per cell

    Fval = np.concatenate(Fval)

    return Fval


def plot_mean_F_val(ax, folders, colors, F_value, summary_paths, summary_protocol = 'Tunings', special_dict = None) :

    x = np.arange(0,len(folders))
    y = []
    y_err = []
    n_cells_pop = []

    for folder, summary_path in zip(folders, summary_paths) : 
        dFoF_parameters = get_dFoF_params(folder, None)
        dFoF_parameters['with_correctedFluo_and_F0'] = True
        nwb_filenames = get_nwb_filenames_from_summary(folder, summary_path, summary_protocol)

        Fval = get_fluorescence_level(nwb_filenames, dFoF_parameters = dFoF_parameters, F_value = F_value)
        y.append(np.mean(Fval))
        y_err.append(stats.sem(Fval))
        n_cells_pop.append(len(Fval))

    ax.errorbar(x = x, y = y, yerr = y_err, fmt='.', elinewidth=1, capthick=1,  capsize = 5, color = 'black')
    ax.bar(x = x, height = y, width = 0.7, color = colors)
    ax.set_ylabel('valid cells mean ' + F_value, fontsize = 13)
    ax.set_xticks([])

    ymax = ax.get_ylim()[1]
    for i in range(len(folders)) : 
        ax.text(s=' N=\n'+str(n_cells_pop[i]),  x = x[i] - 0.1, y = ymax/22, rotation = 'horizontal', fontsize = 10, color = 'white')

#%% 
#---------SUPPORT FUNCTIONS----------#

def add_axis_labels_to_plot(ax, plot_type = 'Tunings', fontsize = 13):

    if plot_type == 'Tunings' : 
        ticks = [-22.5,   0. ,  22.5,  45. ,  67.5,  90. , 112.5, 135. ]
        ax.set_xticks(ticks = ticks, labels = ['%i' % a if (a in [0, 90]) else '' for a in ticks], fontsize = fontsize)
        ax.set_ylabel('norm. $\\delta$ $\\Delta$F/F', fontsize = fontsize)
        ax.set_xlabel('angle ($^o$) from pref.', fontsize = fontsize)

    if plot_type == 'Sensitivities' : 
        ax.set_xticks(ticks = np.arange(3)*0.5)
        ax.set_ylabel('$\\delta$ $\\Delta$F/F', fontsize = fontsize)
        ax.set_xlabel('contrast', fontsize = fontsize)

          
def add_legend_to_plot(ax, labels, colors, loc = 'lower left', bbox_to_anchor = None, fontsize = 10) :  #ax can be fig 

    if len(colors[0]) == 2 and len(colors) == len(labels) : #include lightgrey control cond in legend
        labels = np.concatenate([[label]*2 for label in labels]).tolist()
        colors = [(c[0],c[1]) for c in colors]

    handles = []
    for i, label in enumerate(labels) : 
        patch = mpatches.Patch(color=colors[i], label=label)
        handles.append(patch)

    ax.legend(handles = handles, labels = labels, loc = loc, bbox_to_anchor = bbox_to_anchor, fontsize = fontsize)


