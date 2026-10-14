#%%
import os, sys
import numpy as np  
os.chdir('/home/user/lab-notebook/astrid')

import matplotlib as matplotlib
import matplotlib.colors as mcolors
import string

#%% coherence between all can be better : ie moins redondant 

def plot_corrfact_legend_color_maps(colors, fig, axes, neuropil_correction_factors =  [0,0.15,0.3,0.55,0.7,0.85,1]) :

    for ax,color in zip(axes,colors) :

        if type(color) == str : 
            color = (*mcolors.to_rgb(color),1.0)

        list_color_with_alphas = []
        color_as_list = list(color)

        for k in range(len(neuropil_correction_factors)): 
            color_as_list[-1] = 0.1 + (0.9*k/len(neuropil_correction_factors))
            list_color_with_alphas.append(tuple(color_as_list))

        cmap = mcolors.LinearSegmentedColormap.from_list("", list_color_with_alphas)
        norm = mcolors.Normalize(vmin=0, vmax=1)
        mappable = matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap)
        fig.colorbar(mappable = mappable,  cax=ax, orientation="horizontal", ticks = neuropil_correction_factors, shrink = 0.5 )


def create_multiple_corrfact_mosaic(N, alphabet_upper, alphabet_lower) : 

    """
    Choice of function to create an automatic iterative mosaic on N 
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

def get_input_multiple_corrfact_plot(folder, colors, summary_path, neuropil_correction_factors = [0,0.15,0.3,0.55,0.7,0.85,1], inclufact = '1.15', prefixe = ''): 

    summary_prefixe_names, alphas = get_multiple_corrfact_summary_prefixes_and_alphas(neuropil_correction_factors = neuropil_correction_factors, inclufact = inclufact, prefixe = prefixe)

    special_dict = {'summary_prefixe_name' : summary_prefixe_names, 'ylims': [0,1.05]}

    folders = [folder]*len(neuropil_correction_factors)
    colors = [colors]*len(neuropil_correction_factors)
    summary_paths = [summary_path]*len(neuropil_correction_factors)

    return folders, colors, summary_paths, alphas, special_dict


def get_multiple_corrfact_summary_prefixes_and_alphas(neuropil_correction_factors = [0,0.15,0.3,0.55,0.7,0.85,1], inclufact = '1.15', prefixe = ''): 

    summary_prefixe_names = []
    alphas = []
    for j, neuropil_correction_factor in enumerate(neuropil_correction_factors) :
        summary_prefixe_names.append(prefixe + 'roi_to_neuropil_fluo_inclusion_factor_' +inclufact+ '_neuropil_correction_factor_' + str(neuropil_correction_factor) + '_')
        alphas.append(0.1 + j*(0.9/(len(neuropil_correction_factors) -1)))

    return summary_prefixe_names, alphas

