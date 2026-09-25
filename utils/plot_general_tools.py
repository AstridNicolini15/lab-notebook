#%%
import os ,sys
import numpy as np 
from scipy.stats import sem

os.chdir('/home/user/lab-notebook/astrid')
sys.path += ['./physion/src']
sys.path += ['./summary_plots']
import physion.utils.plot_tools as pt
from physion.analysis.episodes.build import EpisodeData
from physion.analysis.protocols.orientation_tuning import get_tuning_responses, fit_gaussian 
from physion.analysis.protocols.contrast_sensitivity import *
import physion
from utils.propagated_uncertainty import *

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import to_rgb



#%%

def draw_tuning_curve(key, 
                    special_dict = {}, 
                    summary_path = '/home/user/DATA/Astrid/summary',
                    ax = None, 
                    color = None,
                    alpha = 1,
                    xy = (0,0),
                    draw_uncertainty = False,
                    graph_width_dict = {'lw' : 2, 'ms' : 3, 'fontsize' : 7}) : 
    
    if 'name' in special_dict.keys() : 
        Tunings = np.load(summary_path + '/' + special_dict['name'] + 'Tunings_%s.npy' % key, allow_pickle=True) 
    else : 
        Tunings = np.load(summary_path + '/Tunings_%s.npy' % key, allow_pickle=True)  

    Responses = get_tuning_responses(Tunings, average_by='sessions') #return only significant resp
    Responses = [r/r[1] for r in Responses]

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

    ax.scatter(x_angles, np.mean(Responses, axis=0), 
            color=color, alpha = alpha)

    ax.plot(x, func(x), lw=lw, color=color, alpha = alpha)

    ax.annotate(text = 'N= ' + str(n_cells), xy = xy, color = color, fontsize = fontsize, alpha = alpha)

    if draw_uncertainty : 
        if 'std-values' not in Tunings[0].keys(): 
            print('no std values, uncertainty is not drawn')
        else : 
            #uncertainty_sy = session_sem_with_indepedance_hypothesis_universal(Tunings)
            uncertainty_sy = sem(Responses, axis = 0)
            ax.errorbar(Tunings[0]['shifted_angle'], 
                    np.nanmean(Responses, axis=0),
                    yerr=uncertainty_sy,
                    elinewidth = 1,
                    fmt = '.',
                    color=color, ms = ms)



def draw_sensitivitie_curve(key,
                            special_dict = {}, 
                            summary_path = '/home/user/DATA/Astrid/summary',
                            ax = None, 
                            color = None,
                            alpha = 1,
                            xy = (0,0),
                            draw_uncertainty = False,
                            graph_width_dict = {'lw' : 2, 'ms' : 3, 'fontsize' : 7}
                            ):

    #load
    if 'name' in special_dict.keys() : 
        Sensitivities = np.load(summary_path + '/' + special_dict['name'] + 'Sensitivities_%s.npy' % key, allow_pickle=True) 
    else : 
        Sensitivities = np.load(summary_path + '/Sensitivities_%s.npy' % key, allow_pickle=True)  

    #get values
    #Responses = get_responses(Sensitivities)
    Responses = [np.mean(S['Responses'][np.sum(S['significant_pos'] + S['significant_neg'], axis = 1).astype(bool)], axis=0)\
                                for S in Sensitivities] #get significant responses

    n_cells = np.sum([np.sum(np.sum(S['significant_pos'] + S['significant_neg'], axis = 1).astype(bool)) for S in Sensitivities]) #TO CHECKs
    if draw_uncertainty : 
        uncertainty_sy = session_sem_with_indepedance_hypothesis_universal(Sensitivities)
    else : 
        uncertainty_sy = None

    #draw
    lw = graph_width_dict['lw']
    ms = graph_width_dict['ms']
    fontsize = graph_width_dict['fontsize']

    pt.plot(Sensitivities[0]['contrast'], 
        np.nanmean(Responses, axis=0), 
        sy=uncertainty_sy,
        color=color,
        ax=ax)
                
    ax.annotate(text = 'N= ' + str(n_cells), xy = xy, color = color, fontsize = fontsize)
    

def add_axis_labels_to_plot(ax, plot_type = 'Tuning', fontsize = 13):

    if plot_type == 'Tuning' : 
        ticks = [-22.5,   0. ,  22.5,  45. ,  67.5,  90. , 112.5, 135. ]
        ax.set_xticks(ticks = ticks, labels = ['%i' % a if (a in [0, 90]) else '' for a in ticks], fontsize = fontsize)
        ax.set_ylabel('norm. $\\delta$ $\\Delta$F/F', fontsize = fontsize)
        ax.set_xlabel('angle ($^o$) from pref.', fontsize = fontsize)

    if plot_type == 'Sensitivities' : 
        ax.set_xticks(ticks = np.arange(3)*0.5)
        ax.set_ylabel('$\\delta$ $\\Delta$F/F', fontsize = fontsize)
        ax.set_xlabel('contrast', fontsize = fontsize)

          
def create_legend(folders, colors, with_grey = False) : 

    handles = []
    for i, folder in enumerate(folders) : 

        patch = mpatches.Patch(color=colors[i][0], label=folder)
        handles.append(patch)

        if with_grey == True : 
            patch = mpatches.Patch(color=colors[i][1], label=folder)
            handles.append(patch)

    return handles


def correct_summary_form(summary_path, folders) :
    if type(summary_path) == str :
        summary_path = [summary_path] * len(folders)
        print('summary_path is not a list, assuming same path for all folders')
    return summary_path

def plot_tuning_responses_many_pop(folders, 
                                   colors, 
                                   ax, 
                                   special_dict = {}, 
                                   summary_path =  [ "/home/user/DATA/Astrid/Cibele_data/summary"]) : 
    
    summary_path = correct_summary_form(summary_path, folders)

    ylims = [0,1.05]
    if bool(special_dict): 
        ax.set_title(special_dict['title'], fontsize = 13)
        ylims = special_dict['ylims']

    for i,folder in enumerate(folders) : 

        keys =  ['%s_contrast-1.0' % folder,
                '%s_contrast-0.5' % folder]
        
        for j, key in enumerate(keys) :

            xy = (115,1-(0.04*j)-(0.08*i))
            draw_tuning_curve(key, 
                    special_dict = special_dict, 
                    summary_path = summary_path[i],
                    ax = ax, 
                    color = colors[i][j],
                    alpha = 1,
                    xy = xy, 
                    draw_uncertainty = True,
                    graph_width_dict = {'lw' : 3, 'ms' : 5, 'fontsize' : 10}) 


    add_axis_labels_to_plot(ax = ax)
    ax.set_ylim(ylims)


def plot_contrast_responses_many_pop(folders, 
                                     colors, 
                                     ax, 
                                     special_dict = {}, 
                                     summary_path = '/home/user/DATA/Astrid/summary') : 

    summary_path = correct_summary_form(summary_path, folders)

    ylims = None
    if bool(special_dict): 
        ax.set_title(special_dict['title'], fontsize = 13)
        ylims = special_dict['ylims']

    for i,folder in enumerate(folders) : 

        keys = ['%s_angle-0.0' % folder, 
            '%s_angle-90.0' % folder]

        for j, key in enumerate(keys) :

            xy = (0.8,0.1-(0.04*j)-(0.08*i))
            draw_sensitivitie_curve(key,
                            special_dict = special_dict, 
                            summary_path = summary_path[i],
                            ax = ax, 
                            color = colors[i][j],
                            alpha = 1,
                            xy = xy,
                            draw_uncertainty = True,
                            graph_width_dict = {'lw' : 2, 'ms' : 3, 'fontsize' : 10})
            
    add_axis_labels_to_plot(ax, plot_type = 'Sensitivities')
    ax.set_ylim(ylims)

#%%
#---------RESPONSIVENESS TO VISUAL STIM----------#
     
def get_responsiveness_to_visual_stim(folder, summary_protocol = 'Tunings', summary_path =  '/home/user/DATA/Astrid/summary') : 

    perc_resp_to_visual = ()
    if summary_protocol == 'Tunings' :
        for key in ['%s_contrast-1.0' % folder, 
            '%s_contrast-0.5' % folder] : 

            Tunings = np.load(summary_path + '/' + summary_protocol + '_' + key + '.npy', allow_pickle=True) 

            n_cell_resp = np.sum([np.sum(Tuning['significant_ROIs']) for Tuning in Tunings])
            n_cell_nonresp = np.sum([np.sum(~Tuning['significant_ROIs']) for Tuning in Tunings])

            perc_resp_to_visual = (*perc_resp_to_visual, (n_cell_resp *100) / (n_cell_resp + n_cell_nonresp))
            perc_resp_to_visual = (*perc_resp_to_visual, (n_cell_nonresp *100) / (n_cell_resp + n_cell_nonresp))


    if summary_protocol == 'Sensitivities' :
        for key in ['%s_angle-90.0' % folder, 
            '%s_angle-0.0' % folder] : 

            Sensitivities = np.load(summary_path + summary_protocol + '_' + key + '.npy', allow_pickle=True) 
            cell_responsiveness = np.concatenate([np.sum(S['significant_pos'] + S['significant_neg'], axis = 1) for S in Sensitivities])

            n_cell_resp = len([x for x in cell_responsiveness if x !=0])
            n_cell_nonresp = len([x for x in cell_responsiveness if x ==0])

            perc_resp_to_visual = (*perc_resp_to_visual, (n_cell_resp *100) / (n_cell_resp + n_cell_nonresp))
            perc_resp_to_visual = (*perc_resp_to_visual, (n_cell_nonresp *100) / (n_cell_resp + n_cell_nonresp))

    return  perc_resp_to_visual


def plot_responsiveness_pie(folders, colors_list, axes, summary_protocol = 'Tunings', summary_path = '/home/user/DATA/Astrid/summary') : 

        
    summary_path = correct_summary_form(summary_path, folders)

    for i,folder in enumerate(folders) : 

        perc_resp_c1, perc_nonresp_c1, perc_resp_c05, perc_nonresp_c05 = np.round(get_responsiveness_to_visual_stim(folder, summary_protocol = summary_protocol, summary_path = summary_path[i]),3)
        
        for j, perc_resp, perc_nonresp in zip([0, 1], [perc_resp_c05, perc_resp_c1], [perc_nonresp_c05, perc_nonresp_c1]):
            axes[i*2+j].pie([perc_resp, perc_nonresp], 
                            colors = colors_list[i%2], 
                            startangle = 90,
                            wedgeprops={"edgecolor":"black",'linewidth': 1, 'width' : 0.6},
                            textprops = {"fontsize" : 10})
            
        axes[i*2].annotate(text = f'{perc_resp_c05:.1f}% \n resp', xy= (-1.2,1), color = colors_list[i%2][0], fontsize = 13)
        axes[i*2+1].annotate(text = f'{perc_resp_c1:.1f}% \n resp', xy= (-1.2,1), color = colors_list[i%2][0], fontsize = 13)
    axes[0].set_title('half contrast', color = 'black', loc = 'left', fontsize = 17, pad = 35)
    axes[1].set_title('full contrast', color = 'black', loc = 'left', fontsize = 17, pad = 35)
                