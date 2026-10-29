#%% gOSI analysis
# cell 1 : cell and population gOSI functions  
# cell 2 : population gOSI at different inclusion factor
# cell 3 : cells's paired gOSI at contrast 1 and 0.5, distribution and wilcoxon stat test 
  
import sys
sys.path += ["/home/user/lab-notebook/astrid/physion/src"]
sys.path += ["/home/user/lab-notebook/astrid"]

from general_tools.cell_type_color_codes import get_working_data_defaults, reshape_colors_with

from physion.analysis.protocols.orientation_tuning import compute_selectivities
from physion.analysis.protocols.orientation_tuning import selectivity_index
import numpy as np
from math import radians
from scipy.stats import wilcoxon, sem 

from matplotlib.transforms import Bbox
import matplotlib.pyplot as plt
#from matplotlib.colors import to_rgb
#import physion.utils.plot_tools as pt

# %%

def get_cell_gOSI(Radians, Responses) : 
    """
    Responses and Radians are both 1D arrays of dim = (1*8)
    """
    cell_gOSI = np.sqrt(np.sum(Responses*np.cos(Radians))**2 + np.sum(Responses*np.sin(Radians))**2) / np.sum(Responses)

    return cell_gOSI


def get_OSIs_or_gOSIs(Tunings, osi_type = 'gOSI', only_significant = True, clip_negative = True) : 

    OSIs = []

    angles = Tunings[0]['shifted_angle']
    radians = [np.deg2rad(2*(angle + 22.5)) for angle in angles]

    for T in Tunings : 
        session_OSIs = []

        Responses = T['Responses']
        if only_significant : 
            Responses = Responses[T['significant_ROIs']]
        if clip_negative : 
            Responses = np.array([[max(x,0) for x in r ] for r in Responses]) 

        for R in Responses : #warning : do not work with summary with 2d significances
            if osi_type == 'gOSI' : 
                cell_OSI = get_cell_gOSI(radians, R)
            elif osi_type == 'OSI':
                cell_OSI = selectivity_index(angles, R)
            session_OSIs.append(cell_OSI)

        OSIs.append(session_OSIs)
    
    return OSIs


def plot_OSI_or_gOSI(ax, summary_filenames, colors, osi_type = 'gOSI', labels = None, x = None, mean_by = 'sessions', clip_negative = True): #with mean and std 

    height = []
    height_error = []
    for filename in summary_filenames : 
        Tunings = np.load(filename, allow_pickle = True)

        OSIs = get_OSIs_or_gOSIs(Tunings, osi_type = osi_type, clip_negative = clip_negative)

        if mean_by == 'cells': 
            OSIs = np.concatenate(OSIs)
            height.append(np.nanmean(OSIs))
            height_error.append(np.nanstd(OSIs))
            ax.set_title('mean by cells')

        elif mean_by == 'sessions':
            OSIs = [np.mean(session_gosi) for session_gosi in OSIs]
            height.append(np.nanmean(OSIs))
            height_error.append(sem(OSIs, nan_policy = 'omit'))

    if x is None : 
        x = np.arange(len(height))

    ax.grid(True, linewidth = 1, axis = 'y')
    ax.bar(x = x, height = height, width = 0.75, color = colors)
    ax.errorbar(x=x, y= height, yerr=height_error, fmt='.', elinewidth=1, capthick=1,  capsize = 5, color = 'black')
    if labels is not None : 
        ax.set_xticks(ticks = x, labels = labels, rotation = 45)
    ax.set_ylabel(osi_type)
    ax.set_xticks([])
    #ax.set_title('negative clipped')


#%% population gOSI at different inclusion factor
if False : 
    """
    In this cell : we plot tuning responses at differents inclusion factors 
        + the corresponding population gOSI (average over session, then average across sessions)

    """

    for experimenter in ['Cibele', 'Taddy'] : 

        summary_path, folders, colors = get_working_data_defaults(experimenters = experimenter,
                                                                folders_surnames = ['SST_WT', 'SST_GluN1_KO'],
                                                                color_with_grey = True)
        
        mosaic = """
            AABBCCDDD
            AABBCCDDD
            """

        fig = plt.figure(layout="constrained", figsize = (20,5))
        ax_dict = fig.subplot_mosaic(mosaic, gridspec_kw={"wspace": 0.3})


        for neuropil_inclusion_factor, ax in zip([1.15, 2, 3], [ax_dict["A"], ax_dict["B"], ax_dict["C"]]) : 
            special_dict = {'name' : 'inclusion_factor-' + str(neuropil_inclusion_factor) + '_correction_factor-0.7_' , 
                            'title' : 'neuropil inclusion \n factor = '+ str(neuropil_inclusion_factor),
                            'ylims' : (0,1.1)}
            
            plot_tuning_responses_many_pop(folders, colors, ax, special_dict = special_dict, summary_path = summary_path, clip_negative = False)

        colors = reshape_colors_with(colors, second_cond_low_alpha = True, second_cond_alpha = 0.3) 
        colors = np.concatenate([[c]*3 for c in np.concatenate(colors)])

        x = [0,1,2,4,5,6,10,11,12,14,15,16]
        height = []

        for i, folder in enumerate(folders) : 

            for key in ['%s_contrast-1.0' % folder,
                        '%s_contrast-0.5' % folder]:

                for inclusion_factor in [1.15, 2, 3] : 

                    prefixe = 'inclusion_factor-' + str(inclusion_factor) + '_correction_factor-0.7_'
                    Tunings = np.load(summary_path[i] + '/' + prefixe + 'Tunings_%s.npy' % key, allow_pickle=True)
                    gOSIs = get_gOSIs(Tunings)
                    height.append(np.mean([np.mean(session_gOSIs) for session_gOSIs in gOSIs])) #average of sessions averages

                    if False : #use_osi_instaed_of_gosi  
                        osi = []
                        for T in Tunings : 
                            SIs = compute_selectivities(T['Responses']/T['Responses'][1],
                                    using='orth-resp',#"fit",
                                    angles=np.linspace(-22.5, 135, 8))
                            osi.append(np.mean(np.array(SIs)[T['significant_ROIs']]))
                        height.append(np.mean(osi))

        ax_dict["D"].grid(True, linewidth = 1, axis = 'y')
        ax_dict["D"].bar(x = x, height = height, width = 0.75, color = colors)
        ax_dict["D"].set_xticks(ticks = x, labels = [1.15,2,3]*4, rotation = 45)
        ax_dict["D"].set_xlabel('Neuropil inclusion factor')
        ax_dict["D"].set_ylabel('gOSI')
        ax_dict["D"].set_ylim(0,0.45)
        ax_dict["D"].set_title('negative clipped')


    #%% cells's paired gOSI at contrast 1 and 0.5, distribution and wilcoxon stat test 

    """
    Hyp : SST cells orientation tuning is contrast dependant. ie. cells's gOSI decrease between contrast 1 and contrast 0.5; 

    In this cell : we plot the cells paired gOSIs at contrast 1 and 0.5 
        + the histogram of the paired difference between gOSIs
        + print the wilcoxon pval of the distribution of the paired difference between gOSIs

    Note : only cells significant at both contrast 1 and 0.5 are included
    """

    for experimenter in ['Cibele', 'Taddy'] : 

        summary_path, folders, colors = get_working_data_defaults(experimenters = experimenter,
                                                                folders_surnames = ['SST_WT', 'SST_GluN1_KO'])

        for i,folder in enumerate(folders) : 

            gOSIs = []
            concatenate_significant_cells = []

            for key in ['%s_contrast-1.0' % folder,
                        '%s_contrast-0.5' % folder]:

                Tunings = np.load(summary_path[i] + '/' + 'Tunings_%s.npy' % key, allow_pickle=True)

                gOSIs.append(np.concatenate(get_gOSIs(Tunings, only_significant = False)))
                concatenate_significant_cells.append(np.concatenate([T['significant_ROIs'] for T in Tunings]))

            significant_cells_in_both_cond = concatenate_significant_cells[0]*concatenate_significant_cells[1]

            diff_gOSI = gOSIs[0][significant_cells_in_both_cond] - gOSIs[1][significant_cells_in_both_cond]
            mean_diff_gOSI = np.mean(diff_gOSI)


            #plot
            fig,ax = plt.subplots(1,2,figsize = (11,7))
            ax[0].set_position(Bbox([[0, 0], [0.7, 1]]))
            ax[1].set_position(Bbox([[0.8, 0], [1, 0.3]]))

            #ax 0 : plot paired cells gOSIs
            for y1,y2 in zip(gOSIs[0][significant_cells_in_both_cond], gOSIs[1][significant_cells_in_both_cond]) : 
                ax[0].plot([0,1], [y1,y2], color = 'black', linewidth = 0.2)
                ax[0].scatter([0,1], [y1,y2], color = colors[i])

            ax[0].set_ylabel('gOSI of cells significant in both contrast')
            ax[0].set_xticks(ticks = [0,1], labels = ['contrast = 1', 'contrast = 0.5'])
            ax[0].set_xlim(-0.15,1.15)
            ax[0].set_ylim(-0.05,0.9)
            ax[0].annotate(text = 'mean diff gOSI \n =' + str(np.round(mean_diff_gOSI, 3)), xy = (1.05,0.06), fontsize = 9)
            ax[0].annotate(text = folder, xy = (1.05,0.80), color = colors[i])
            ax[0].annotate(text = 'N cell = ' + str(np.sum(significant_cells_in_both_cond)), xy = (1.05,0.13), color = colors[i])

            #ax 1 : plot hist and wilcoxon pval
            counts, bins = np.histogram(diff_gOSI, 20)
            ax[1].stairs(counts/len(diff_gOSI), bins, color = colors[i])

            ax[1].set_xlabel('diff gOSI values')
            ax[1].set_ylabel('%')
            ax[1].set_yticks(ticks = [0,0.125,0.25], labels = [0,0.125,0.25])
            #ax[1].set_xlim(-0.8,0.8)
            ax[1].annotate(text = 'Wilcoxon test \np value = \n' + str(wilcoxon(diff_gOSI)[1]), xy = (0,0.25), fontsize = 9)

            plt.show()

