#%%
import os ,sys
os.chdir('/home/user/lab-notebook/astrid')
sys.path += ['./physion/src']
sys.path += ['./summary_plots']
import physion.utils.plot_tools as pt
from physion.analysis.read_NWB\
                         import scan_folder_for_NWBfiles, Data
from physion.analysis.episodes.build import EpisodeData
from physion.analysis.protocols.orientation_tuning import *
from run_rest_responses.tuning_arousal_summary_functions import *
import matplotlib.pyplot as plt
import physion

    #%%



def draw_tuning_curve(key, 
                    special_dict = {}, 
                    summary_path = '/home/user/DATA/Astrid/run_rest_summary',
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

    Responses = get_tuning_responses(Tunings, average_by='sessions')

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

    ax.scatter(x_angles, np.mean([r/r[1] for r in Responses], axis=0), 
            color=color, alpha = alpha)

    ax.plot(x, func(x), lw=lw, color=color, alpha = alpha)

    ax.annotate(text = 'N= ' + str(n_cells), xy = xy, color = color, alpha = alpha, fontsize = fontsize)

    if draw_uncertainty : 
        if 'std-values' not in Tunings[0].keys(): 
            print('no std values, uncertainty is not drawn')
        else : 
            uncertainty_sy = session_sem_with_indepedance_hypothesis_universal(Tunings)
            ax.errorbar(Tunings[0]['shifted_angle'], 
                    np.nanmean(Responses, axis=0),
                    yerr=uncertainty_sy,
                    elinewidth = 2,
                    fmt = '.',
                    color=color, ms = ms)




def add_axis_labels_to_plot(ax, plot_type = 'Tuning'):

    if plot_type == 'Tuning' : 
        ticks = [-22.5,   0. ,  22.5,  45. ,  67.5,  90. , 112.5, 135. ]
        ax.set_xticks(ticks = ticks, labels = ['%i' % a if (a in [0, 90]) else '' for a in ticks], fontsize = 13)
        ax.set_ylabel('norm. $\\delta$ $\\Delta$F/F', fontsize = 13)
        ax.set_xlabel('angle ($^o$) from pref.', fontsize = 13)


