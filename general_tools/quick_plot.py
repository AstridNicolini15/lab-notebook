#%%
import numpy as np
import sys, os 
os.chdir('/home/user/lab-notebook/astrid')
sys.path.insert(0,'/home/user/lab-notebook/astrid')
from general_tools.plot_general_tools import *
from general_tools.cell_type_color_codes import *

#%% choose populations

summary_path, folders, colors= get_working_data_defaults(experimenters = ['Cibele'],
                                                        folders_surnames = ['PYR_WT'], 
                                                        color_with_grey=True)

add_legend = True 

#%%%%%%%%%%%%%%%% TEMPLATE plot and compare FLUORESCENCE LEVEL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ymin = 0
ymax = 200

mosaic = "ABC" 
fig= plt.figure(layout="constrained", figsize = (12,6))
fig.get_layout_engine().set(wspace = 0.15)
ax_dict = fig.subplot_mosaic(mosaic)

if add_legend:
    handles = create_legend(folders, colors)
    fig.legend(handles = handles, loc = 'outside lower left', fontsize = 8)


for letter, inclufact in zip(['A', 'B', 'C'], [1.15,2,3]) : 
        
    dFoF_parameters = dict(\
        roi_to_neuropil_fluo_inclusion_factor = inclufact,
        neuropil_correction_factor = 0.7,
        #roi_specific_neuropil_correction_factor = False,
        method_for_F0 = 'sliding_percentile',
        percentile=5., # percent
        sliding_window = 5*60, # seconds
        with_correctedFluo_and_F0=True,
    )


    plot_mean_F_val(folders, dFoF_parameters, 'correctedFluo0', colors, ax_dict[letter], summary_path = summary_path)


    ax_dict[letter].set_ylim([ymin, ymax])
    ax_dict[letter].set_title('inclusion factor = ' + str(inclufact), fontsize = 10)


#%%%%%%%%%%%%%%%% TEMPLATE  plot and compare multiple TUNING CURVES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mosaic = "A" 
fig= plt.figure(layout="constrained", figsize = (6,6))
fig.get_layout_engine().set(wspace = 0.15)
ax_dict = fig.subplot_mosaic(mosaic)

if add_legend:
    handles = create_legend(folders, colors)
    fig.legend(handles = handles, loc = 'outside lower left', fontsize = 8)

colors = reshape_colors_with(colors, second_cond_grey = True)

plot_tuning_responses_many_pop(folders, colors, ax_dict["A"], summary_path = summary_path) 



#%%%%%%%%%%%%%%%% TEMPLATE plot and compare multiple CORRECTION FACTORS TUNING CURVES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
add_legend = False

colors = reshape_colors_with(colors, second_cond_grey = True)
neuropil_correction_factors = [0,0.15,0.3,0.55,0.7,0.85,1]

fig = plt.figure(figsize = (7*len(folders),10))

axes = get_multiple_corrfact_axes_and_implement_mosaic(folders, fig)


for i in range(len(folders)) : 
    special_dict = {'base_name' : 'inclusion_factor-1.15_correction_factor-'}
    plot_orientation_tuning_multiple_corrfact(folders[i],
                        colors=colors[i],
                        neuropil_correction_factors = neuropil_correction_factors,
                        ax = axes[i][0],
                        special_dict = special_dict,
                        summary_path = summary_path)

    axes[i][1].set_title('neuropil correction factors', fontsize = 13)
    #axes[i][0].set_title('neuropil inclusion factor = 4', fontsize = 8)
    plot_corrfact_legend_color_maps(colors = colors[i], fig = fig,  axes = [axes[i][1], axes[i][2]], neuropil_correction_factors = neuropil_correction_factors) 


colors = reshape_colors_with(colors)
if add_legend:
    handles = create_legend(folders, colors)
    fig.legend(handles = handles, loc = 'outside lower left', fontsize = 8)

#%%%%%%%%%%%%%%%% TEMPLATE plot and compare gOSI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#%% COMPOSED TEMPLATE profil response

mosaic = """
    AAAAAABBD
    AAAAAABBD
    AAAAAACCD
    AAAAAACCD
    """

fig = plt.figure(layout="constrained", figsize = (10,7))
fig.get_layout_engine().set(wspace = 0.15)
for i, folder in enumerate(folders) : 
    fig.text(x= 0, y= -0.05-0.03*i, s = '%s' % folder, color = colors[i][0], fontsize = 17)

ax_dict = fig.subplot_mosaic(mosaic)

#plot tuning resp
plot_tuning_responses_many_pop(folders, colors, ax_dict["A"], summary_path = summary_path) 


#plot resp to visual pies
plot_responsiveness_pie(folders, colors, [ax_dict["B"], ax_dict["C"]], summary_path = summary_path)


colors_without_grey = reshape_colors_with(colors)
plot_mean_F_val(folders, 'correctedFluo0', colors_without_grey, ax_dict["D"], summary_path = summary_path)



#%%  NOTEBOOK NEUROPIL AT HIGH FLUORESCENCE
#from notebooks.Tuning_multiple_corrfacts import get_multiple_corrfact_axes_and_implement_mosaic, plot_orientation_tuning_multiple_corrfact, plot_corrfact_legend_color_maps

colors = reshape_colors_with(colors, second_cond_grey = True)
neuropil_correction_factors = [0, 0.15,0.3,0.55,0.7,0.85,1]
inclusion_factors = [1.15,4,5]
summary_path = summary_path*len(inclusion_factors)
folders = folders*len(inclusion_factors)
colors = colors*len(inclusion_factors)

fig = plt.figure(figsize = (7*len(folders),10))

axes = get_multiple_corrfact_axes_and_implement_mosaic(folders, fig)


for i, inclufact in enumerate(inclusion_factors) : 
    special_dict = {'base_name' : 'inclusion_factor-' + str(inclufact) + '_correction_factor-'}
    plot_orientation_tuning_multiple_corrfact(folders[i],
                        colors=colors[i],
                        neuropil_correction_factors = neuropil_correction_factors,
                        ax = axes[i][0],
                        special_dict = special_dict,
                        summary_path = summary_path)
    axes[i][0].set_title('inclusion factor =' + str(inclufact))
    axes[i][1].set_title('neuropil correction factors', fontsize = 13)
    plot_corrfact_legend_color_maps(colors = colors[i], fig = fig,  axes = [axes[i][1], axes[i][2]], neuropil_correction_factors = neuropil_correction_factors) 


#%%

mosaic = "D" 
fig= plt.figure(layout="constrained", figsize = (6,6))
fig.get_layout_engine().set(wspace = 0.15)
ax_dict = fig.subplot_mosaic(mosaic)

if add_legend:
    handles = create_legend(folders, colors)
    fig.legend(handles = handles, loc = 'outside lower left', fontsize = 8)



inclusion_factors = [1.15, 2, 3, 4, 5]
correction_factors = [0,0.15,0.3,0.55,0.7,0.85,1] 

folder = folders[0]
summary_path = summary_path * len(inclusion_factors)
colors = reshape_colors_with(colors, second_cond_low_alpha = True, second_cond_alpha = 0.3) 
colors = np.concatenate([[c]*len(inclusion_factors) for c in colors[0] ])
colors = np.concatenate([colors]*len(correction_factors)) 

LEN = len(inclusion_factors)

x = []
xticks = [7+(i*15) for i in range(len(correction_factors))]


temp_x = [1,2,3,4,5,7,8,9,10,11]
x = []
last_x = -3
for i in range(len(correction_factors)) : 
    x.append([item + (last_x +4) for item in temp_x])
    last_x = x[-1][-1]
x = np.concatenate(x).tolist()


height = []

for corrfact in correction_factors : 
    
    for key in ['%s_contrast-1.0' % folder,
                '%s_contrast-0.5' % folder]:

        for i, inclusion_factor in enumerate(inclusion_factors) : 
            prefixe = 'inclusion_factor-' + str(inclusion_factor) + '_correction_factor-' + str(corrfact) + '_'
            Tunings = np.load(summary_path[i] + '/' + prefixe + 'Tunings_%s.npy' % key, allow_pickle=True)
            gOSIs = get_gOSIs(Tunings)
            height.append(np.mean([np.mean(session_gOSIs) for session_gOSIs in gOSIs])) #average of sessions averages


ax_dict["D"].grid(True, linewidth = 1, axis = 'y')
ax_dict["D"].bar(x = x, height = height, width = 0.75, color = colors)
ax_dict["D"].set_xticks(ticks = xticks, labels = correction_factors, rotation = 45)
ax_dict["D"].set_xlabel('Neuropil correction factor')
ax_dict["D"].set_ylabel('gOSI')
ax_dict["D"].set_ylim(0,0.40)
ax_dict["D"].set_title('negative clipped')
sec = ax_dict["D"].secondary_xaxis(location=0)
sec.set_xticks([2,3,4,5,6,8,9,10,11,12], labels=[str(item) for item in inclusion_factors*2], fontsize = 5, rotation = 90)
