#%%
import sys
sys.path += ["/home/user/lab-notebook/astrid"]

from plot_general_tools import *
import matplotlib as matplotlib
import matplotlib.colors as mcolors

#%% Cibele %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

folders = ["PV-cells_WT_Adult_V1", 
    "SST-cells_WT_Adult_V1",
    #"SST-cells_cond-GluN1-KO_Adult_V1",
    ]

summary_path = ["/home/user/DATA/Astrid/Cibele_data/summary"]*len(folders)

colors = [[(*to_rgb('#b30a7bff'), 1.0),'lightgrey'],  [(*to_rgb('#12522eff'), 1.0),'lightgrey']]
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

#%%

def plot_orientation_tuning_multiple_corrfact(folder,
                      colors=None,
                      neuropil_correction_factors =  [0,0.15,0.3,0.55,0.7,0.85,1],
                      ax = None,
                      ylims = (0,1.05),
                      special_dict = None,
                      summary_path = '/home/user/DATA/Astrid/summary'):


    keys = ['%s_contrast-1.0' % folder, 
                '%s_contrast-0.5' % folder]

    for iter_x, (key, color) in enumerate(zip(keys, colors)):

        for iter_y,neuropil_correction_factor in enumerate(neuropil_correction_factors) :

            special_dict['name'] = special_dict['base_name']  + str(neuropil_correction_factor) + '_'

            xy = (110+(25*iter_x),0.22-(0.03*iter_y))

            draw_tuning_curve(key, 
                    special_dict = special_dict, 
                    summary_path = summary_path,
                    ax = ax, 
                    color = color,
                    alpha = 0.3 + iter_y/10,
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


def create_multiple_corrfact_mosaic(N, alphabet_upper, alphabet_lower) : 

    """
    N = len(folders)
    """
    mosaic =  alphabet_upper[0]*12 
    for i in range(1,N) : 
        mosaic += alphabet_lower[25+1-i]*2
        mosaic += alphabet_upper[i]*12 
    mosaic = [mosaic + '\n   ']*10
    mosaic += [alphabet_lower[25+1-(N -1)-1]*(12*N + 2*(N -1)) + '\n   ']*2

    legend = alphabet_upper[N]*12
    for i in range(N+1, N*3) : 
        legend += alphabet_lower[i-1]*2
        legend += alphabet_upper[i]*12

    final_legend = [legend[:(len(legend)//2) -1] + '\n   ']
    final_legend += [legend[(len(legend)//2) +1:] + '\n   '] 

    mosaic += final_legend

    final_mosaic = ""
    final_mosaic = final_mosaic.join(mosaic)
    mosaic = final_mosaic

    return mosaic


def get_multiple_corrfact_axes_and_implement_mosaic(folders, fig): 

    import string
    alphabet_upper = list(string.ascii_uppercase)
    alphabet_lower = list(string.ascii_lowercase)

    mosaic = create_multiple_corrfact_mosaic(len(folders), alphabet_upper, alphabet_lower) 

    ax_dict = fig.subplot_mosaic(mosaic,  gridspec_kw={"wspace": 50,"hspace": 3})

    separators = np.unique([letter for letter in mosaic if letter.islower()])
    for letter in separators :
        ax_dict[letter].axis("off")

    axes = []
    for i in range(len(folders)) : 
        axes.append([ax_dict[alphabet_upper[i]], ax_dict[alphabet_upper[len(folders)+i]], ax_dict[alphabet_upper[len(folders)*2+i]]])

    return axes