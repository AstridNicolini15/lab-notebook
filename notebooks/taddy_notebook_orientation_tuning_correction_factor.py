#%% Orientation tuning responses, responsiveness and neuropil correction factor impact 
from notebooks.Tuning_multiple_corrfacts import * 

#%%

summary_path = "/home/user/DATA/Astrid/Taddy_data/summary"

folders = [
    "SST_WT",
    "SST_GluN1",
    "SST_GluN3"
    ]

colors = [[pt.tab10(1),'lightgrey'], [pt.tab10(4),'lightgrey'], [pt.tab10(5),'lightgrey']]

#%%

summary_path = "/home/user/DATA/Astrid/Taddy_data/photodiode_aligned/summary"

folders = [
    "SST_WT",
    "SST_GluN1",
    "SST_GluN3"
    ]

colors = [[pt.tab10(1),(*pt.tab10(1)[:3],0.3)], [pt.tab10(4),(*pt.tab10(4)[:3],0.3)], [pt.tab10(5),(*pt.tab10(5)[:3],0.3)]]


#%% plot orientation tuning and responsiveness

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

#%% plot neuropil correction factor impact 


colors = [[pt.tab10(1),'lightgrey'], [pt.tab10(4),'lightgrey'], [pt.tab10(5),'lightgrey']]
neuropil_correction_factors = [0,0.15,0.3,0.55,0.7,0.85,1]

fig = plt.figure(figsize = (7*len(folders),10))
for i, folder in enumerate(folders) : 
    fig.text(x= 0.37+0.23*i, y= 1, s = '%s' % folder, color = colors[i][0], fontsize = 15)

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
    plot_corrfact_legend_color_maps(colors = colors[i], fig = fig,  axes = [axes[i][1], axes[i][2]], neuropil_correction_factors = neuropil_correction_factors) 
