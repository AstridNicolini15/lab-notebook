#%%import os, sys , shutil 
import sys
sys.path += ["/home/user/lab-notebook/astrid"]
from general_tools.cell_type_color_codes import get_working_data_defaults
from general_tools.plot_general_tools import plot_tuning_responses_many_pop, plot_responsiveness_pie, plot_mean_F_val
import matplotlib.pyplot as plt

#%%Cibele  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
summary_path, folders, colors= get_working_data_defaults(experimenters = 'Cibele',
                                                        folders_surnames = ['PV_WT'],
                                                        color_with_grey = True)

#%%Taddy  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
summary_path, folders, colors= get_working_data_defaults(experimenters = 'Taddy',
                                                        folders_surnames = ['SST_WT', 'SST_GluN1_KO'],
                                                        color_with_grey = True)

summary_path, folders, colors= get_working_data_defaults(experimenters = 'Cibele',
                                                        folders_surnames = ['PYR_WT'],
                                                        color_with_grey = True)

#%%
mosaic = """
    AAAAAABBDD
    AAAAAABBDD
    AAAAAACCEE
    AAAAAACCEE
    IIhhhhhhhh
    IIFFFFGGGG
    IIFFFFGGGG
    IIFFFFGGGG
    """

fig = plt.figure(layout="constrained", figsize = (13,13))
for i, folder in enumerate(folders) : 
    fig.text(x= 0, y= -0.05-0.03*i, s = '%s' % folder, color = colors[i][0], fontsize = 17)

ax_dict = fig.subplot_mosaic(mosaic)
ax_dict["h"].axis("off")

#plot tuning resp
plot_tuning_responses_many_pop(folders, colors, ax_dict["A"], summary_path = summary_path) 


#plot resp to visual pies
plot_responsiveness_pie(folders, colors, [ax_dict["B"], ax_dict["C"],ax_dict["D"], ax_dict["E"]], summary_path = summary_path)

colors_without_grey = reshape_colors_with(colors)
plot_mean_F_val(folders, 'correctedFluo0', colors_without_grey, ax_dict["I"], summary_path = summary_path)

for neuropil_inclusion_factor, ax in zip([2,3], [ax_dict["F"], ax_dict["G"]]) : 
    special_dict = {'name' : 'inclusion_factor-' + str(neuropil_inclusion_factor) + '_correction_factor-0.7_' , 
                    'title' : 'neuropil inclusion \n factor = '+ str(neuropil_inclusion_factor),
                    'ylims' : (0,1.1)}
    
    plot_tuning_responses_many_pop(folders, colors, ax, special_dict = special_dict, summary_path = summary_path)
