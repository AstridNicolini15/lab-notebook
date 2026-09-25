#%% gOSI analysis
# cell 1 : cell and population gOSI functions  
# cell 2 : population gOSI at different inclusion factor
# cell 3 : cells's paired gOSI at contrast 1 and 0.5, distribution and wilcoxon stat test 
  
import sys
sys.path += ["/home/user/lab-notebook/astrid/physion/src"]

from physion.analysis.protocols.orientation_tuning import compute_selectivities

import numpy as np
from math import radians
from scipy import stats

from matplotlib.transforms import Bbox
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgb
import physion.utils.plot_tools as pt


from utils.plot_general_tools import plot_tuning_responses_many_pop
from utils.cell_type_color_codes import get_working_data_defaults, reshape_colors_with

# %%

def get_gOSIs(Tunings, only_significant = True) : 

    gOSIs = []

    angles = Tunings[0]['shifted_angle']
    Radians = [radians(2*(angle + 22.5)) for angle in angles]

    for T in Tunings : 
        session_gOSIs = []

        if only_significant : 
            Responses = T['Responses'][T['significant_ROIs']]
        else : 
            Responses = T['Responses']

        for R in Responses : #to note : do not work with summary with 2d significances

            cell_gOSI = get_cell_gOSI(R, Radians)
            session_gOSIs.append(cell_gOSI)

        gOSIs.append(session_gOSIs)
    
    return gOSIs


def get_cell_gOSI(Responses, Radians) : 

    """
    Responses and Radians are both 1D arrays of dim = (1*8)
    """

    cell_gOSI = np.sqrt(np.sum(Responses*np.cos(Radians))**2 + np.sum(Responses*np.sin(Radians))**2) / np.sum(np.abs(Responses))

    return cell_gOSI

#%% population gOSI at different inclusion factor

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
        
        plot_tuning_responses_many_pop(folders, colors, ax, special_dict = special_dict, summary_path = summary_path)

    colors = reshape_colors_with(colors, second_cond_low_alpha = True, second_cond_alpha = 0.3) 
    colors = np.concatenate([[c]*3 for c in np.concatenate(colors)])

    x = [0,1,2,4,5,6,10,11,12,14,15,16]
    height = []

    for folder in folders : 

        for key in ['%s_contrast-1.0' % folder,
                    '%s_contrast-0.5' % folder]:

            for inclusion_factor in [1.15, 2, 3] : 

                prefixe = 'inclusion_factor-' + str(inclusion_factor) + '_correction_factor-0.7_'
                Tunings = np.load(summary_path + '/' + prefixe + 'Tunings_%s.npy' % key, allow_pickle=True)
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


#%% cells's paired gOSI at contrast 1 and 0.5, distribution and wilcoxon stat test 

"""
Hyp : SST cells orientation tuning is contrast dependant. ie. cells's gOSI decrease between contrast 1 and contrast 0.5; 

In this cell : we plot the cells paired gOSIs at contrast 1 and 0.5 
    + the histogram of the paired difference between gOSIs
    + print the wilcoxon pval of the distribution of the paired difference between gOSIs

Note : only cells significant at contrast 1 and 0.5 are included
"""

for experimenter in ['Cibele', 'Taddy'] : 

    summary_path, folders, colors = get_working_data_defaults(experimenters = experimenter,
                                                            folders_surnames = ['SST_WT', 'SST_GluN1_KO'])

    for i,folder in enumerate(folders) : 

        gOSIs = []
        concatenate_significant_cells = []

        for key in ['%s_contrast-1.0' % folder,
                    '%s_contrast-0.5' % folder]:

            Tunings = np.load(summary_path + '/' + 'Tunings_%s.npy' % key, allow_pickle=True)

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
        ax[1].set_xlim(-0.8,0.8)
        ax[1].annotate(text = 'Wilcoxon test \np value = \n' + str(stats.wilcoxon(diff_gOSI)[1]), xy = (0,0.25), fontsize = 9)

        plt.show()
