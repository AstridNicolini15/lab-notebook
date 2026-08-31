#%%

from summary_plots.WT_vs_KO_Tuning_summary_plot import plot_tuning_responses_many_pop, get_gaussian_fit_and_uncertainty

#%%


folders = [#"PV-cells_WT_Adult_V1", 
    "SST-cells_WT_Adult_V1",
    #"SST-cells_cond-GluN1-KO_Adult_V1"
    ]

mosaic = """
    AB
    """

fig = plt.figure(figsize = (16,9))
fig.subplots_adjust(wspace=0.3, hspace=0)
ax_dict = fig.subplot_mosaic(mosaic)


colors = [[pt.tab10(1),'lightgrey'], [pt.tab10(2),'lightgrey']]
special_dict = {'name' : '', 'title' : '', 'ylims' : (-0.15,1.1)}
plot_tuning_responses_many_pop(folders, colors, ax_dict["A"], special_dict, summary_path = '/home/user/DATA/Astrid/deconvolved_summary_test') 


#plot tuning resp
colors = [[pt.tab10(1),'lightgrey'], [pt.tab10(2),'lightgrey']]
special_dict = {'name' : 'Deconvolved_correctedFluo', 'title' : 'Deconvolved', 'ylims' : (-0.15,1.1)}
    
plot_tuning_responses_many_pop(folders, colors, ax_dict["B"], special_dict, summary_path = '/home/user/DATA/Astrid/deconvolved_summary') 

#%%


folders = [#"PV-cells_WT_Adult_V1", 
    #"SST-cells_WT_Adult_V1",
    "SST-cells_cond-GluN1-KO_Adult_V1"
    ]

mosaic = """
    ABCD
    """

fig = plt.figure(figsize = (16,9))
fig.subplots_adjust(wspace=0.3, hspace=0)
ax_dict = fig.subplot_mosaic(mosaic)

#plot tuning resp
colors = [[pt.tab10(2),'lightgrey']]
special_dict = {'name' : '', 'title' : 'Deconvolved', 'ylims' : (-0.15,1.1)}
    
plot_tuning_responses_many_pop(folders, colors, ax_dict["A"], special_dict, summary_path = '/home/user/DATA/Astrid/run_rest_summary') 

colors = [['forestgreen','lightgrey']]
special_dict = {'name' : 'Deconvolved_', 'title' : 'Deconvolved_dFoF', 'ylims' : (-0.15,1.1)}
    
plot_tuning_responses_many_pop(folders, colors, ax_dict["B"], special_dict, summary_path = '/home/user/DATA/Astrid/deconvolved_summary') 


colors = [['mediumseagreen','lightgrey']]
special_dict = {'name' : 'Deconvolved_correctedFluo', 'title' : 'Deconvolved_correctedFluo', 'ylims' : (-0.15,1.1)}
    
plot_tuning_responses_many_pop(folders, colors, ax_dict["C"], special_dict, summary_path = '/home/user/DATA/Astrid/deconvolved_summary_test') 



colors = [['lightseagreen','lightgrey']]
special_dict = {'name' : 'Deconvolved_rawFluo', 'title' : 'Deconvolved_rawFluo', 'ylims' : (-0.15,1.1)}
    
plot_tuning_responses_many_pop(folders, colors, ax_dict["D"], special_dict, summary_path = '/home/user/DATA/Astrid/deconvolved_summary_test') 

