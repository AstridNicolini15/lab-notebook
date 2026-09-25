#%%
import sys
sys.path += ["/home/user/lab-notebook/astrid"]

from plot_general_tools import plot_tuning_responses_many_pop
from responsiveness_and_brightness import plot_responsiveness_pie, plot_mean_F_val
from utils.cell_type_color_codes import get_working_data_defaults
#%%
summary_path_cibele, folders_cibele, colors_cibele = get_working_data_defaults(experimenters = 'Cibele',
                                                        folders_surnames = ['SST_WT', 'SST_GluN1_KO'],
                                                        color_with_grey = True)




summary_path_taddy, folders_taddy, colors_taddy = get_working_data_defaults(experimenters = 'Taddy',
                                                        folders_surnames = ['SST_WT', 'SST_GluN1_KO'],
                                                        color_with_grey = True)


#%%

mosaic = """
    AAAAAAhBBBBBB
    AAAAAAhBBBBBB
    AAAAAAhBBBBBB
    AAAAAAhBBBBBB
    AAAAAAhBBBBBB
    ooooooooooooo
    EEFFjIIJJkCCC
    EEFFjIIJJkCCC
    KKLLjMMNNkCCC
    KKLLjMMNNkCCC
    """

fig = plt.figure(layout="constrained", figsize = (16,12))
ax_dict = fig.subplot_mosaic(mosaic)

ax_dict["h"].axis("off")
ax_dict["j"].axis("off")
ax_dict["o"].axis("off")
ax_dict["k"].axis("off")


for k, ax_letters, folders, summary_path, colors in zip([0,1],
                                                        [["A", "E", "F","K", "L"], ["B", "I", "J", "M", "N"]], 
                                                        [folders_cibele, folders_taddy], 
                                                        [summary_path_cibele, summary_path_taddy], 
                                                        [colors_cibele, colors_taddy]): 


    for i, folder in enumerate(folders) : 
        fig.text(x= 0+0.45*k, y= -0.02-0.04*i, s = '%s' % folder, color = colors[i][0], fontsize = 15)
    
    plot_tuning_responses_many_pop(folders, colors, ax_dict[ax_letters[0]], special_dict =  {'title' : '', 'ylims' : (0,1.05)}, summary_path = summary_path)

    axes = [ax_dict[ax_letters[1]], ax_dict[ax_letters[3]], ax_dict[ax_letters[2]], ax_dict[ax_letters[4]]]
    plot_responsiveness_pie(folders, colors, axes, summary_path = summary_path)


folders = folders_cibele + folders_taddy
summary_path = summary_path_cibele + summary_path_taddy
colors = colors_cibele + colors_taddy

colors_without_grey = []
for color in colors : 
    colors_without_grey.append(color[0])
plot_mean_F_val(folders, 'correctedFluo0', colors_without_grey, ax_dict["C"], summary_path = summary_path)

#%%Cibele vs Taddy  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

summary_path, folders, colors = get_working_data_defaults(experimenters = ['Cibele', 'Taddy'],
                                                        folders_surnames = ['SST_WT', 'SST_WT'],
                                                        color_with_grey = True)

mosaic = """
    AAAAAABBCC
    AAAAAADDEE
    AAAAAAFFGG

    """

fig = plt.figure(layout="constrained", figsize = (13,7))

ax_dict = fig.subplot_mosaic(mosaic)

plot_tuning_responses_many_pop(folders, colors, ax_dict["A"], summary_path = summary_path) 
plot_responsiveness_pie(folders, colors, [ax_dict["B"], ax_dict["C"],ax_dict["D"], ax_dict["E"], ax_dict["F"], ax_dict["G"]], summary_path = summary_path)

handles = create_legend(folders, colors, with_grey = False)
ax_dict["A"].legend(handles, folders, bbox_to_anchor = (0.7,0.5), fontsize = 10)
