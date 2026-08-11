#%%


# Source - https://stackoverflow.com/a/46778420
# Posted by ImportanceOfBeingErnest, modified by community. See post 'Timeline' for change history
# Retrieved 2026-07-23, License - CC BY-SA 4.0

list_color_with_alphas = []
color = list(pt.tab10(1))
for neuropil_fact in  [0,0.15,0.3,0.55,0.7,0.85,1] : 
    color[-1] = neuropil_fact
    list_color_with_alphas.append(tuple(color))

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors

x,y,c = zip(*np.random.rand(30,3)*4-2)

norm=plt.Normalize(-2,2)
cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", list_color_with_alphas)

plt.scatter(x,y,c=c, cmap=cmap, norm=norm)
plt.colorbar()
plt.show()

#%%
folders = ["PV-cells_WT_Adult_V1", 
    "SST-cells_WT_Adult_V1",
    #"SST-cells_cond-GluN1-KO_Adult_V1"
    ]

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
ax_dict = fig.subplot_mosaic(mosaic,  gridspec_kw={"wspace": 50,"hspace": 3})
ax_dict["H"].axis("off")
ax_dict["h"].axis("off")
ax_dict["i"].axis("off")


plot_orientation_tuning_curve_multiple_corrfact_one_graph(folders[0],
                      colors=[pt.tab10(0), 'lightgrey'],
                      neuropil_correction_factors =  [0,0.15,0.3,0.55,0.7,0.85,1],
                      ax = ax_dict["A"])

lightgrey = (*matplotlib.colors.to_rgb('lightgrey'),1.0)
ax_dict["c"].set_title('neuropil correction factors', fontsize = 13)
plot_color_maps(colors = [pt.tab10(0), lightgrey],fig = fig,  axes = [ax_dict["c"],ax_dict["d"]], neuropil_correction_factors =  [0,0.15,0.3,0.55,0.7,0.85,1]) 



plot_orientation_tuning_curve_multiple_corrfact_one_graph(folders[1],
                      colors=[pt.tab10(1), 'lightgrey'],
                      neuropil_correction_factors =  [0,0.15,0.3,0.55,0.7,0.85,1],
                      ax = ax_dict["B"])

lightgrey = (*matplotlib.colors.to_rgb('lightgrey'),1.0)
ax_dict["e"].set_title('neuropil correction factors', fontsize = 13)
plot_color_maps(colors = [pt.tab10(1), lightgrey],fig = fig,  axes = [ax_dict["e"],ax_dict["f"]], neuropil_correction_factors =  [0,0.15,0.3,0.55,0.7,0.85,1]) 

#%%
def plot_orientation_tuning_curve_multiple_corrfact_one_graph(folder,
                      colors=None,
                      neuropil_correction_factors =  [0,0.15,0.3,0.55,0.7,0.85,1],
                      ax = None,
                      ylims = (-0.05,1.05),
                      summary_path = '/home/user/DATA/Astrid/summary_neuropil_factor'):

    x = np.linspace(-30, 180-30, 100)
    keys = ['%s_contrast-1.0' % folder, 
                '%s_contrast-0.5' % folder]

    for i, (key, color) in enumerate(zip(keys, colors)):

        for k,neuropil_correction_factor in enumerate(neuropil_correction_factors) :
            # load data
            Tunings = np.load(summary_path + '/corrfact_' + str(neuropil_correction_factor) + 'Tunings_%s.npy' % key, 
                        allow_pickle=True)

            Responses = get_tuning_responses(Tunings,
                                           average_by='sessions')

            # Gaussian Fit
            C, func = fit_gaussian(Tunings[0]['shifted_angle'],
                                np.mean([r/r[1] for r in Responses], axis=0))


            ax.scatter(Tunings[0]['shifted_angle'], np.mean([r/r[1] for r in Responses], axis=0), 
                        color=color, alpha = 0.3+k/10)

            ax.plot(x, func(x), lw=2, color=color, alpha =0.3+ k/10)

            n_cells = np.sum([np.sum(t['significant_ROIs']) for t in Tunings])

            ax.annotate(text = 'N= ' + str(n_cells), xy = (110+(25*i),0.2-(0.03*k)), color = colors[i], alpha = 0.3+k/10, fontsize = 7)

    ax.set_xticks(ticks = Tunings[0]['shifted_angle'], labels = ['%i' % a if (a in [0, 90]) else '' for a in Tunings[0]['shifted_angle']], fontsize = 13)
    ax.set_ylabel('norm. $\\delta$ $\\Delta$F/F', fontsize = 13)
    ax.set_xlabel('angle ($^o$) from pref.', fontsize = 13)
    ax.set_ylim(ylims)



#HERE
def plot_color_maps(colors, fig, axes, neuropil_correction_factors =  [0,0.15,0.3,0.55,0.7,0.85,1]) :

    for ax,color in zip(axes,colors) :
        list_color_with_alphas = []
        color_as_list = list(color)

        for k,neuropil_correction_factor in enumerate(neuropil_correction_factors): 
            color_as_list[-1] = 0.3+ k/10
            list_color_with_alphas.append(tuple(color_as_list))

        cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", list_color_with_alphas)
        norm = matplotlib.colors.Normalize(vmin=0, vmax=1)
        mappable = matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap)
        fig.colorbar(mappable = mappable,  cax=ax, orientation="horizontal", ticks = neuropil_correction_factors, shrink = 0.5 )

