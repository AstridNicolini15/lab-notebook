#%%
import os ,sys
import numpy as np 

sys.path.insert(0,'/home/user/lab-notebook/astrid/physion/src')
sys.path.insert(0,'/home/user/lab-notebook/astrid')
from general_tools.cell_type_color_codes import get_working_data_defaults
from general_tools.general_tools import *

import matplotlib.pyplot as plt
#%%

summary_path, folders, colors = get_working_data_defaults(experimenters = ['Cibele','Cibele'],
                                                        folders_surnames = ['SST_WT', "SST_GluN1_KO"],
                                                        color_with_grey=True)

summary_path, folders, colors = get_working_data_defaults(experimenters = ['Taddy', 'Taddy'],
                                                        folders_surnames = ['SST_WT', "SST_GluN1_KO"],
                                                        color_with_grey=True)

mosaic = """
    AAAAhBBBB
    AAAAhBBBB
    AAAAhBBBB
    AAAAhBBBB
    CCEEhGGII
    CCEEhGGII
    DDFFhHHJJ
    DDFFhHHJJ
    """


fig = plt.figure(layout="constrained", figsize = (14,12))

for i, folder in enumerate(folders) : #annotate with legend and mean neuropil correction factor
    fig.text(x = 0, y= -0.05-0.03*i, s = folder, color = colors[i][0], fontsize = 17)

    fig.text(x = 0.6, y= -0.02, s = 'mean neuropil correction factor : ', fontsize = 13)
    for j, key in enumerate(['%s_contrast-1.0' % folder, 
        '%s_contrast-0.5' % folder]) : 

        Tunings = np.load(summary_path[i] + '/allen_neuropil_' + 'Tunings_' + key + '.npy', allow_pickle=True) 
        mean_r =  np.mean(np.concatenate([T['neuropil_correction_factors'] for T in Tunings]))
        fig.text(x = 0.62+0.1*j, y= -0.05-0.03*i, s = str(np.round(mean_r,3)), color = colors[i][j], fontsize = 17)


ax_dict = fig.subplot_mosaic(mosaic)
ax_dict["h"].axis("off")
ax_dict["A"].grid(True, linewidth = 1, axis = 'y')
ax_dict["B"].grid(True, linewidth = 1, axis = 'y')
    
special_dict = {'name' : '',
                'title': '',
                'ylims' : [0,1.1]}
plot_tuning_responses_many_pop(folders, colors, ax_dict["A"], special_dict = special_dict, summary_path = summary_path) 
plot_responsiveness_pie(folders, colors, [ax_dict["C"], ax_dict["D"], ax_dict["E"], ax_dict["F"]], special_dict = special_dict, summary_path = summary_path)


special_dict = {'name' : 'allen_neuropil_',
                'title': 'allen neuropil',
                'ylims' : [0,1.1]}
plot_tuning_responses_many_pop(folders, colors, ax_dict["B"], special_dict = special_dict, summary_path = summary_path) 
plot_responsiveness_pie(folders, colors, [ax_dict["G"], ax_dict["H"],ax_dict["I"], ax_dict["J"]], special_dict = special_dict, summary_path = summary_path)

