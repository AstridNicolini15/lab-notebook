#%%
import numpy as np
import os, sys 
os.chdir('/home/user/lab-notebook/astrid')
sys.path.insert(0,'/home/user/lab-notebook/astrid/physion/src')

from matplotlib.colors import to_rgb
import pandas as pd 
import physion.utils.plot_tools as pt

#%%

data_from_list = ['Cibele']*13 + ['Taddy']*2

base_path = ["/home/user/DATA/Astrid/Cibele_data"]*13 + ["/home/user/DATA/Astrid/Taddy_data/OneDrive_1_9-3-2026"]*2

folders = [
    "SST-cells_WT_Adult_V1",
    "SST-cells_cond-GluN1-KO_Adult_V1",
    "PV-cells_WT_Adult_V1", 
    "PV-cells_cond-GluN1-KO_Adult_V1",
    "PYR-cells_WT_Adult_V1",
    "SST-cells_WT_P15-P19_V1",
    "SST-cells_WT_P20-P23_V1",
    "SST-cells_WT_P24-P27_V1",
    "PV-cells_WT_P20-P23_V1",
    "PV-cells_WT_P24-P27_V1",
    "PYR-cells_WT_P15-P19_V1",
    "PYR-cells_WT_P20-P23_V1",
    "PYR-cells_WT_P24-P27_V1",
    "Wild-Type",
    "GluN1-KO"
    ]

folders_surnames = [
    "SST_WT",
    "SST_GluN1_KO",
    "PV_WT", 
    "PV_GluN1_KO",
    "PYR_WT",
    "SST_WT_P15",
    "SST_WT_P20",
    "SST_WT_P24",
    "PV_WT_P20",
    "PV_WT_P24",
    "PYR_WT_P15",
    "PYR_WT_P20",
    "PYR_WT_P24",
    "SST_WT",
    "SST_GluN1_KO"
    ]


colors = [
    (*to_rgb('#12522eff'), 1),
    (*to_rgb('#5d1490ff'), 1),
    (*to_rgb('#b30a7bff'), 1),
    (*to_rgb('#9e0004ff'), 1),
    (*to_rgb('#005d8aff'), 1),
    (*to_rgb('#c6e7bbff'), 1),
    (*to_rgb('#99b993ff'), 1),
    (*to_rgb('#428e58ff'), 1),
    (*to_rgb('#ff9ab5ff'), 1),
    (*to_rgb('#e31593ff'), 1),
    (*to_rgb('#b8cfe0ff'), 1),
    (*to_rgb('#6caaccff'), 1),
    (*to_rgb('#008fbbff'), 1),
    pt.tab10(1), 
    pt.tab10(4), 
    ]



data = {'experimenter' : data_from_list, 'base_path' : base_path,'folder_name' : folders, 'folder_surname' : folders_surnames, 'color' : colors }
working_data_defaults = pd.DataFrame(data)



def get_working_data_defaults(experimenters = ['Cibele'], 
                              folders_surnames = ['SST_WT'], 
                              color_with_grey = False, 
                              alpha = 1,
                              summary_path_list = False) : 


    #correcting and checking input form 
    if type(experimenters) == str: 
        experimenters = [experimenters]

    if type(folders_surnames) == str: 
        folders_surnames = [folders_surnames]

    if len(experimenters) != len(folders_surnames) : 
        print("experimenters and folders surnames length do not match \n assuming %s experimenter for all folders " % experimenters[0])
        experimenters = [experimenters[0]] * len(folders_surnames)

    #collecting relevant values
    summary_folder = []
    folders = []
    colors = []
    for experimenter, folder_surname in zip(experimenters, folders_surnames) : 

        cond = (working_data_defaults['experimenter'] == experimenter) & \
                (working_data_defaults['folder_surname'] == folder_surname)
        
        summary_folder.append(working_data_defaults['base_path'][cond].values.item() + '/summary')
        folders.append(working_data_defaults['folder_name'][cond].values.item())
        colors.append(working_data_defaults['color'][cond].values.item())

    if alpha != 1 :
        colors = [(*c[:3], alpha) for c in colors]

    if color_with_grey : 
        colors = [[c, 'lightgrey'] for c in colors]

    return summary_folder, folders, colors 


def reshape_colors_with(colors,
                        alpha = 1,
                       second_cond_grey = False,
                       second_cond_low_alpha = False, 
                       second_cond_alpha = 0.3,
                       grey_to_rgb = False) : 

    #check and initialisation
    if second_cond_grey & second_cond_low_alpha :
        print('second_cond_grey and second_cond_low_alpha are mutually exclusive')
        return

    if len(colors[0]) == 2  : 
        colors = [c[0] for c in colors] #keep only cell type color, remove control cond color
        if (not second_cond_grey) & (not second_cond_low_alpha) & (not grey_to_rgb): 
            print('warning : second cond color removed')

    #do reshaping
    if alpha != 1 :
        colors = [(*c[:3], alpha) for c in colors]

    if second_cond_grey : 
        colors = [[c, 'lightgrey'] for c in colors]

    if second_cond_low_alpha : 
        colors = [[c, (*c[:3], second_cond_alpha)] for c in colors]

    if grey_to_rgb : 
        colors = [[c, (*to_rgb('lightgrey'),1)] for c in colors]

    return colors

