#%%
import numpy as np
import sys, os 
os.chdir('/home/user/lab-notebook/astrid')
sys.path.insert(0,'/home/user/lab-notebook/astrid')
from general_tools.plot_general_tools import *
from general_tools.cell_type_color_codes import *
from create_summaries.Preprocessing_Settings import get_dFoF_params
import string
alphabet_upper = list(string.ascii_uppercase)

#%% CORE PLOTS TEMPLATES : also working with populations comparaisons
#%% choose populations

summary_paths, folders, colors= get_working_data_defaults(experimenters = ['Cibele'],
                                                        folders_surnames = ['SST_WT'], 
                                                        color_with_grey=True)


if False:
    add_legend_to_plot(fig, folders, colors, loc = 'outside lower left') 

#%%
def create_figure(mosaic, figsize, hspace = 0, wspace = 0): 
    fig= plt.figure(layout="constrained", figsize = figsize)
    fig.get_layout_engine().set(hspace = hspace, wspace = wspace)
    ax_dict = fig.subplot_mosaic(mosaic)
    return fig, ax_dict

#%%%%%%%%%%%%%%%% FLUORESCENCE LEVELS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mosaic = "A" 
fig, ax_dict = create_figure(mosaic, figsize = (2*len(folders),5))
colors = reshape_colors_with(colors)
plot_mean_F_val(ax_dict['A'], folders, colors, 'F0', summary_paths)

ymin = 0
ymax = 200
ax_dict['A'].set_ylim([ymin, ymax])


#%%%%%%%%%%%%%%%% TUNING CURVES OR CONTRAST RESPONSES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

summary_protocol = 'Tunings' #'Sensitivities
ylims = [0,1.05] #[-0.15,0.7]

mosaic = "A" 
fig, ax_dict = create_figure(mosaic, figsize = (5,5))
colors = reshape_colors_with(colors, second_cond_low_alpha = True, second_cond_alpha=0.3)
summary_filenames = get_summary_names(folders, summary_paths, summary_protocols = summary_protocol)

plot_responses(ax_dict['A'],
                summary_filenames,
                colors,  
                summary_protocol = summary_protocol,
                clip_negative = False,
                ylims = ylims)

#%%%%%%%%%%%%%%%%  RESPONSIVENESS PIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mosaic = "".join(alphabet_upper[:len(folders)*2])
fig, ax_dict = create_figure(mosaic, figsize = (2*len(folders),5))
#plot resp to visual pies
summary_filenames = get_summary_names(folders, summary_paths, summary_protocols = 'Tunings')
colors = reshape_colors_with(colors, second_cond_low_alpha=True)
plot_responsiveness_pies(list(ax_dict.values()), summary_filenames, colors, summary_protocol = 'Tunings', total_on = 'nROIs_final')



#%%%%%%%%%%%%%%%% OSIs and gOSIs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
from notebooks.gOSI import plot_OSI_or_gOSI 

mosaic =  """
    A
    """
fig, ax_dict = create_figure(mosaic, figsize = (5,5))

colors = np.concatenate(reshape_colors_with(colors, second_cond_low_alpha = True, second_cond_alpha=0.3))
summary_filenames = get_summary_names(folders, summary_paths, summary_protocols = 'Tunings')
plot_OSI_or_gOSI(ax_dict['A'], summary_filenames, colors, osi_type = 'gOSI', labels = None, x = None, mean_by = 'sessions')

