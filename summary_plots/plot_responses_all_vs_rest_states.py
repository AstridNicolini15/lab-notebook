from common_fcts import *

#%%
summary_path =  "/home/user/DATA/Astrid/Cibele_data/summary"

folders = [#"PV-cells_WT_Adult_V1", 
    "SST-cells_WT_Adult_V1",
    "SST-cells_cond-GluN1-KO_Adult_V1"
    ]

colors = [[(to_rgb('#12522eff'), 0.8),'lightgrey'], [(to_rgb('#5d1490ff'), 0.8),'lightgrey']]

#%%

mosaic = """
    AB
    CD
    """


fig = plt.figure(layout="constrained", figsize = (10,10))
ax_dict = fig.subplot_mosaic(mosaic,  gridspec_kw={"wspace": 0.2,"hspace": 0.15})

fig.text(x = 0.2, y = 1.05, s = 'All states', fontsize = 16)
fig.text(x = 0.7, y = 1.05, s = 'Rest state only', fontsize = 16)


special_dict = {'name' : '', 'title' : '', 'ylims' : (0,1.1)}
plot_tuning_responses_many_pop(folders, colors,  ax_dict["A"], special_dict = special_dict, summary_path = summary_path)
special_dict = {'name' : 'Rest_', 'title' : '', 'ylims' : (0,1.1)}
plot_tuning_responses_many_pop(folders, colors,  ax_dict["B"], special_dict = special_dict, summary_path = summary_path)

special_dict = {'name' : '', 'ylabel' : '$\\delta$ $\\Delta$F/F', 'ylims' : (-0.15,0.55), 'title' : ''}
#plot_contrast_responses_many_pop(folders, colors, ax_dict["C"], plot_only_sign = True)
plot_contrast_responses_many_pop(folders, colors, ax_dict["C"], special_dict = special_dict, summary_path= summary_path) 


special_dict = {'name' : 'Rest_', 'ylabel' : '$\\delta$ $\\Delta$F/F', 'ylims' : (-0.15,0.55), 'title' : ''}
#plot_contrast_responses_many_pop(folders, colors, ax_dict["D"], plot_only_sign = True, special_dict = special_dict)
plot_contrast_responses_many_pop(folders, colors, ax_dict["D"], special_dict = special_dict, summary_path=summary_path) 

#%%

