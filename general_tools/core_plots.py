#%%
import numpy as np
import sys, os 
os.chdir('/home/user/lab-notebook/astrid')
sys.path.insert(0,'/home/user/lab-notebook/astrid')
from general_tools.plot_general_tools import *
from general_tools.cell_type_color_codes import *
from create_summaries.Preprocessing_Settings import get_dFoF_params

#%% CORE PLOTS TEMPLATES : also working with populations comparaisons
#%% choose populations

summary_paths, folders, colors= get_working_data_defaults(experimenters = ['Cibele'],
                                                        folders_surnames = ['SST_WT'], 
                                                        color_with_grey=True)

add_legend = False
#if add_legend:
#    add_legend_to_plot(fig, folders, colors, loc = 'outside lower left') 


#%%%%%%%%%%%%%%%% FLUORESCENCE LEVELS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mosaic = "A" 
fig= plt.figure(layout="constrained", figsize = (2*len(folders),6))
fig.get_layout_engine().set(wspace = 0.15)
ax_dict = fig.subplot_mosaic(mosaic)

plot_mean_F_val(ax_dict['A'], folders, colors, 'correctedFluo0', summary_paths = summary_paths)

ymin = 0
ymax = 200
ax_dict['A'].set_ylim([ymin, ymax])


#%%%%%%%%%%%%%%%% TUNING CURVES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mosaic = "A" 
fig= plt.figure(layout="constrained", figsize = (6,6))
fig.get_layout_engine().set(wspace = 0.15)
ax_dict = fig.subplot_mosaic(mosaic)

colors = reshape_colors_with(colors, second_cond_low_alpha = True, second_cond_alpha=0.3)

plot_responses(ax_dict['A'],
                folders, 
                colors,  
                summary_protocol = 'Tunings',
                clip_negative = False,
                special_dict = {'ylims' : [0,1.05]}, 
                summary_paths = summary_paths)


#%%%%%%%%%%%%%%%% CONTRAST RESPONSES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

summary_paths, folders, colors= get_working_data_defaults(experimenters = ['Cibele'],
                                                        folders_surnames = ['SST_WT', 'SST_GluN1_KO'], 
                                                        color_with_grey=True)

mosaic = "A" 
fig= plt.figure(layout="constrained", figsize = (6,6))
fig.get_layout_engine().set(wspace = 0.15)
ax_dict = fig.subplot_mosaic(mosaic)

plot_responses(ax_dict['A'],
                folders, 
                colors,  
                summary_protocol = 'Sensitivities',
                clip_negative = True,
                special_dict = {'ylims' : [-0.15,0.7]}, 
                summary_paths = summary_paths)


#%%%%%%%%%%%%%%%%  RESPONSIVENESS PIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mosaic =  """
    ABCD
    """
fig= plt.figure(layout="constrained", figsize = (6,6))
fig.get_layout_engine().set(wspace = 0.15)
ax_dict = fig.subplot_mosaic(mosaic)

colors = reshape_colors_with(colors)
#plot resp to visual pies
plot_responsiveness_pies(list(ax_dict.values()), folders, colors,  summary_paths = summary_paths,  summary_protocol = 'Tunings', total_on = 'nROIs_final')



#%%%%%%%%%%%%%%%% OSIs and gOSIs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
from notebooks.gOSI import plot_OSI_or_gOSI

mosaic =  """
    A
    """
fig= plt.figure(layout="constrained", figsize = (6,6))
fig.get_layout_engine().set(wspace = 0.15)
ax_dict = fig.subplot_mosaic(mosaic)

colors = np.concatenate(reshape_colors_with(colors, second_cond_low_alpha = True, second_cond_alpha=0.3))
filenames = []
for folder, summary_path in zip(folders, summary_paths) :
    keys = get_folder_keys(folder, summary_protocol = 'Tunings') 
    for key in keys : 
        filenames.append(summary_path + '/Tunings_' + key + '.npy') 

plot_OSI_or_gOSI(ax_dict['A'], filenames, colors, osi_type = 'gOSI', labels = None, x = None, mean_by = 'sessions')



#%%%%%%%%%%%%%%%% MULTIPLE CORRECTION FACTORS TUNING OR SENSITIVITIE CURVES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
from notebooks.multiple_corrfacts import *

neuropil_correction_factors = [0,0.15,0.3,0.55,0.7,0.85,1]

fig = plt.figure(figsize = (7*len(folders),10))
for i, folder in enumerate(folders) : 
    fig.text(x= 0.37+0.23*i, y= 1, s = '%s' % folder, color = colors[i][0], fontsize = 15)

axes = get_multiple_corrfact_axes_and_implement_mosaic(folders, fig)


for i, (folder, summary_path, folder_colors) in enumerate(zip(folders, summary_paths, colors)) : 
    synthetic_folders, synthetic_colors, synthetic_summary_paths, alphas, special_dict  = get_input_multiple_corrfact_plot(folder, folder_colors, summary_path)
    plot_responses(axes[i][0],
                synthetic_folders, 
                synthetic_colors,  
                alphas = alphas,
                summary_protocol = 'Tunings', #'Sensitivities',
                clip_negative = False,
                special_dict = special_dict, 
                summary_paths = synthetic_summary_paths,
                #sensitivitie_uncertainty=None
                ) 

    axes[i][1].set_title('neuropil correction factors', fontsize = 13)
    plot_corrfact_legend_color_maps(colors = colors[i], fig = fig,  axes = [axes[i][1], axes[i][2]], neuropil_correction_factors = neuropil_correction_factors) 


if True : #plot associated gOSI (one folder at the time)

    summary_paths, folders, colors= get_working_data_defaults(experimenters = ['Cibele'],
                                                        folders_surnames = ['SST_WT'], 
                                                        color_with_grey=True)

    mosaic =  """
        A
        """
    fig= plt.figure(layout="constrained", figsize = (6,6))
    fig.get_layout_engine().set(wspace = 0.15)
    ax_dict = fig.subplot_mosaic(mosaic)

    summary_prefixe_names, alphas = get_multiple_corrfact_summary_prefixes_and_alphas()
    filenames = []
    for folder, summary_path in zip(folders, summary_paths) :
        keys = get_folder_keys(folder, summary_protocol = 'Tunings') 
        for key in keys : 
            for summary_prefixe_name in summary_prefixe_names : 
                filenames.append(summary_path + '/' + summary_prefixe_name + 'Tunings_' + key + '.npy') 
        
    x = [0,1,3,4,6,7,9,10,12,13,15,16,18,19]
    labels = np.concatenate([[str(val)]*2 for val in neuropil_correction_factors])
    colors = reshape_colors_with(colors, second_cond_low_alpha = True, second_cond_alpha=0.3)[0]
    colors = np.concatenate([colors for i in range(len(neuropil_correction_factors))])
    plot_OSI_or_gOSI(ax_dict['A'], filenames, colors, osi_type = 'gOSI', labels = labels, x = x, mean_by = 'sessions')
    
    ax_dict['A'].set_xlabel('neuropil correction factors')


#%%%%%%%%%%%%%%%% ALL STATES VS REST RESPONSES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#can be generalize to multiple folders ? 
mosaic = """
    AB
    CD
    """

fig= plt.figure(layout="constrained", figsize = (10,10))
fig.get_layout_engine().set(wspace = 0.15)
ax_dict = fig.subplot_mosaic(mosaic)

summary_prefixe_names = ['', 'Rest_', '', 'Rest_']
ylims = [[0,1.05],[0,1.05],[-0.1,0.6],[-0.1,0.6]]
for letter, summary_protocol, summary_prefixe_name, ylim in zip(['A','B','C','D'], 
                                                                ['Tunings']*2 + ['Sensitivities']*2,
                                                                summary_prefixe_names,
                                                                ylims) : 
    special_dict = {'summary_prefixe_name' : summary_prefixe_name, 'title' : '', 'ylims' : ylim}
    plot_responses(ax_dict[letter],
                    folders, 
                    colors,  
                    summary_protocol = summary_protocol,
                    clip_negative = True,
                    special_dict = special_dict, 
                    summary_paths = summary_paths)

for letter in ['A','B','C','D'] : 
    ax_dict[letter].grid(True, linewidth = 1, axis = 'y')

