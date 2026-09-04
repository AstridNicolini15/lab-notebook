#%%
from common_fcts import *
import matplotlib as matplotlib
import matplotlib.colors as mcolors

#%% Cibele %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

folders = ["PV-cells_WT_Adult_V1", 
    "SST-cells_WT_Adult_V1",
    #"SST-cells_cond-GluN1-KO_Adult_V1",
    ]

summary_path = ["/home/user/DATA/Astrid/Cibele_data/summary"]*len(folders)

colors = [[pt.tab10(0),'lightgrey'],  [(*to_rgb('#12522eff'), 1.0),'lightgrey']]
#colors = [[(*to_rgb('#12522eff'), 1.0),'lightgrey'], [(*to_rgb('#5d1490ff'), 1.0),'lightgrey']]

#%%Cibele vs Taddy  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

summary_path = ["/home/user/DATA/Astrid/Cibele_data/summary", "/home/user/DATA/Astrid/Taddy_data/OneDrive_1_9-3-2026/summary"]
#summary_path = ["/home/user/DATA/Astrid/Taddy_data/OneDrive_1_9-3-2026/summary",  "/home/user/DATA/Astrid/Taddy_data/OneDrive_1_9-3-2026/summary"]

folders = ["SST-cells_WT_Adult_V1",
    "Wild-Type",
    #"GluN1-KO",
    ]

colors = [ [(*to_rgb('#12522eff'), 1.0),'lightgrey'], [pt.tab10(1),'lightgrey']]

#colors = [[pt.tab10(1),'lightgrey'], [pt.tab10(4),'lightgrey']]

#%%
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
    fig.text(x= 0.37+0.4*i, y= 1, s = '%s' % folder, color = colors[i][0], fontsize = 15)

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
                        summary_path = summary_path[i])

    axes[i][1].set_title('neuropil correction factors', fontsize = 13)
    plot_corrfact_legend_color_maps(colors = colors[i], fig = fig,  axes = [axes[i][1], axes[i][2]], neuropil_correction_factors = neuropil_correction_factors) 


#%%

def plot_orientation_tuning_multiple_corrfact(folder,
                      colors=None,
                      neuropil_correction_factors =  [0,0.15,0.3,0.55,0.7,0.85,1],
                      ax = None,
                      ylims = (-0.05,1.05),
                      special_dict = None,
                      summary_path = '/home/user/DATA/Astrid/summary'):


    keys = ['%s_contrast-1.0' % folder, 
                '%s_contrast-0.5' % folder]

    for iter_x, (key, color) in enumerate(zip(keys, colors)):

        for iter_y,neuropil_correction_factor in enumerate(neuropil_correction_factors) :
            print(neuropil_correction_factor)
            special_dict['name'] = special_dict['base_name']  + str(neuropil_correction_factor) + '_'

            xy = (110+(25*iter_x),0.2-(0.03*iter_y))

            draw_tuning_curve(key, 
                    special_dict = special_dict, 
                    summary_path = summary_path,
                    ax = ax, 
                    color = color,
                    alpha = 0.3+ iter_y/10,
                    xy = xy)

    add_axis_labels_to_plot(ax = ax)
    ax.set_ylim(ylims)


def plot_corrfact_legend_color_maps(colors, fig, axes, neuropil_correction_factors =  [0,0.15,0.3,0.55,0.7,0.85,1]) :

    for ax,color in zip(axes,colors) :

        if type(color) == str : 
            color = (*mcolors.to_rgb(color),1.0)

        list_color_with_alphas = []
        color_as_list = list(color)

        for k,neuropil_correction_factor in enumerate(neuropil_correction_factors): 
            color_as_list[-1] = 0.3+ k/10
            list_color_with_alphas.append(tuple(color_as_list))

        cmap = mcolors.LinearSegmentedColormap.from_list("", list_color_with_alphas)
        norm = mcolors.Normalize(vmin=0, vmax=1)
        mappable = matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap)
        fig.colorbar(mappable = mappable,  cax=ax, orientation="horizontal", ticks = neuropil_correction_factors, shrink = 0.5 )

