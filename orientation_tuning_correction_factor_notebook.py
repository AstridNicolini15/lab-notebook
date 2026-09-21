#%% Orientation tuning responses, responsiveness and neuropil correction factor impact 
from orientation_tuning_plot_tools import * 

#%%

summary_path = "/home/user/DATA/Astrid/Taddy_data/OneDrive_1_9-3-2026/summary"

folders = [
    "Wild-Type",
    "GluN1-KO"
    ]

colors = [[pt.tab10(1),'lightgrey'], [pt.tab10(4),'lightgrey']]

#%% plot orientation tuning and responsiveness

mosaic = """
    AAAAAABBCC
    AAAAAABBCC
    AAAAAADDEE
    AAAAAADDEE
    """

fig = plt.figure(layout="constrained", figsize = (13,7))

ax_dict = fig.subplot_mosaic(mosaic)

plot_tuning_responses_many_pop(folders, colors, ax_dict["A"], summary_path = summary_path) 
plot_responsiveness_pie(folders, colors, [ax_dict["B"], ax_dict["C"],ax_dict["D"], ax_dict["E"]], summary_path = summary_path)


#%% plot neuropil correction factor impact 
neuropil_correction_factors = [0,0.15,0.3,0.55,0.7,0.85,1]

mosaic = """
    AAAAAAAAAAAAhhBBBBBBBBBBBB
    AAAAAAAAAAAAhhBBBBBBBBBBBB
    AAAAAAAAAAAAhhBBBBBBBBBBBB
    AAAAAAAAAAAAhhBBBBBBBBBBBB
    AAAAAAAAAAAAhhBBBBBBBBBBBB
    AAAAAAAAAAAAhhBBBBBBBBBBBB
    AAAAAAAAAAAAhhBBBBBBBBBBBB
    AAAAAAAAAAAAhhBBBBBBBBBBBB
    AAAAAAAAAAAAhhBBBBBBBBBBBB
    AAAAAAAAAAAAhhBBBBBBBBBBBB
    HHHHHHHHHHHHHHHHHHHHHHHHHH
    HHHHHHHHHHHHHHHHHHHHHHHHHH
    cccccccccccciieeeeeeeeeeee
    ddddddddddddiiffffffffffff
    """


fig = plt.figure(figsize = (14,10))
for i, folder in enumerate(folders) : 
    fig.text(x= 0.37+0.35*i, y= 1, s = '%s' % folder, color = colors[i][0], fontsize = 15)

ax_dict = fig.subplot_mosaic(mosaic,  gridspec_kw={"wspace": 50,"hspace": 3})
ax_dict["H"].axis("off")
ax_dict["h"].axis("off")
ax_dict["i"].axis("off")

axes = [[ax_dict["A"], ax_dict['c'], ax_dict['d']],[ax_dict["B"], ax_dict['e'],ax_dict['f']]]


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
