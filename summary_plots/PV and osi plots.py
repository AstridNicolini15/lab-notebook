#%%

folders = [
    "PV-cells_WT_Adult_V1", 
     "PV-cells_cond-GluN1-KO_Adult_V1", 
]

colors = [[(to_rgb('#b30a7bff'), 0.8),'lightgrey'], 
          [(to_rgb('#9e0004ff'), 0.8),'lightgrey']
          ]

base_path =  '/home/user/DATA/Astrid/Cibele_data'

summary_folder = base_path + '/summary'

fig,ax = plt.subplots(figsize = (7,7))
special_dict = {'title' : 'interval post = [1,2]', 'ylims' : None}
#special_dict = {'name' : 'int_post_0_1_', 'title' : 'interval post = [0,1]', 'ylims' : None}

plot_tuning_responses_many_pop(folders, colors, ax, special_dict = special_dict, summary_path = summary_folder) 

#%%
folders = [
    "PV-cells_WT_Adult_V1", 
     "PV-cells_cond-GluN1-KO_Adult_V1", 
]

colors = [[(to_rgb('#b30a7bff'), 0.8),'lightgrey'], 
          [(to_rgb('#9e0004ff'), 0.8),'lightgrey']
          ]

base_path =  '/home/user/DATA/Astrid/Cibele_data'

summary_folder = base_path + '/summary'

mosaic = """
    ABC
    """

fig = plt.figure(layout="constrained", figsize = (15,6))
ax_dict = fig.subplot_mosaic(mosaic)
special_dict = {'name' : '', 'ylims' : None, 'title' : 'inclusion factor = 1.15 '}
plot_tuning_responses_many_pop(folders, colors, ax_dict["A"], special_dict = special_dict, summary_path = summary_folder) 


for i, inclusion_factor in enumerate([2,3]) : 
    ax = [ax_dict["B"], ax_dict["C"]][i]
    special_dict = {'name' : 'inclusion_factor-%i_correction_factor-0.7_' % inclusion_factor, 'ylims' : None, 'title' : 'inclusion factor = %i ' % inclusion_factor}
    plot_tuning_responses_many_pop(folders, colors, ax, special_dict = special_dict, summary_path = summary_folder) 




#%%

summary_folder = '/home/user/DATA/Astrid/Cibele_data/summary'

folders = [#"PV-cells_WT_Adult_V1", 
    "SST-cells_WT_Adult_V1",
    "SST-cells_cond-GluN1-KO_Adult_V1",
    ]

colors = [to_rgb('#12522eff')]*3 + [(to_rgb('#12522eff'), 0.3)]*3 + [to_rgb('#5d1490ff')]*3 + [(to_rgb('#5d1490ff'), 0.3)]*3
#%%

summary_folder = "/home/user/DATA/Astrid/Taddy_data/OneDrive_1_9-3-2026/summary"

folders = [
    "Wild-Type",
    "GluN1-KO"
    ]

colors = [pt.tab10(1)]*3 + [(pt.tab10(1)[:3], 0.3)]*3 + [pt.tab10(4)]*3 + [(pt.tab10(4)[:3], 0.3)]*3
alphas = [1]*6 + [0.2]*6 

#%%

x = [0,1,2,4,5,6,10,11,12,14,15,16]
height = []
labels = []
for folder in folders : 

    for key in ['%s_contrast-1.0' % folder,
                '%s_contrast-0.5' % folder]:

        for inclusion_factor in [1.15, 2, 3] : 
            labels.append(key)
            prefixe = 'inclusion_factor-' + str(inclusion_factor) + '_correction_factor-0.7_'
            Tunings = np.load(summary_folder + '/' + prefixe + 'Tunings_%s.npy' % key, allow_pickle=True)
            osi = [T['selectivities'][T['significant_ROIs']] for T in Tunings]
            height.append(np.mean(np.concatenate(osi)))

#%%

# Plot the legend with only those 4 items


fig, ax = plt.subplots(figsize = (7,7))
#plt.grid(True, linewidth = 1, axis = 'y')
ax.bar(x = x, height = height, width = 0.75, color = colors, label = labels )
ax.set_xticks(ticks = x, labels = [1.15,2,3,1.15,2,3,1.15,2,3,1.15,2,3], rotation = 45)
ax.set_xlabel('Neuropil inclusion factor')
ax.set_ylabel('OSI')

handles, labels = ax.get_legend_handles_labels()

# Choose 4 specific indices to show (e.g., indices 0, 3, 7, 11)
selected_indices = [0, 3, 6, 9]
selected_handles = [handles[i] for i in selected_indices]
selected_labels =  [labels[i] for i in selected_indices]

ax.legend(selected_handles, selected_labels, bbox_to_anchor=(0.4, -0.15))
