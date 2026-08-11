#%%

folders = [#"PV-cells_WT_Adult_V1", 
    "SST-cells_WT_Adult_V1",
    "SST-cells_cond-GluN1-KO_Adult_V1"
    ]


mosaic = """
    AB
    CD
    """


fig = plt.figure(layout="constrained", figsize = (10,10))
ax_dict = fig.subplot_mosaic(mosaic,  gridspec_kw={"wspace": 0.2,"hspace": 0.15})

fig.text(x = 0.2, y = 1.05, s = 'All states', fontsize = 16)
fig.text(x = 0.7, y = 1.05, s = 'Rest state only', fontsize = 16)

colors = [[pt.tab10(1),'lightgrey'], [pt.tab10(2),'lightgrey']]

special_dict = {'name' : '', 'title' : '', 'ylims' : (-0.125,1.1)}
plot_tuning_responses_many_pop(folders, colors,  ax_dict["A"], special_dict = special_dict)
special_dict = {'name' : 'Rest_', 'title' : '', 'ylims' : (-0.125,1.1)}
plot_tuning_responses_many_pop(folders, colors,  ax_dict["B"], special_dict = special_dict)

special_dict = {'name' : '', 'ylabel' : '$\\delta$ $\\Delta$F/F', 'ylims' : (-0.15,0.55), 'title' : ''}
plot_contrast_responses_many_pop(folders, colors, ax_dict["C"], plot_only_sign = True)

special_dict = {'name' : 'Rest', 'ylabel' : '$\\delta$ $\\Delta$F/F', 'ylims' : (-0.15,0.55), 'title' : ''}
#plot_contrast_responses_many_pop(folders, colors, ax_dict["D"], plot_only_sign = True, special_dict = special_dict)
plot_rest_contrast_responses_many_pop_significance_modified(folders, colors, ax_dict["D"], special_dict = special_dict) 

#%%


def plot_rest_contrast_responses_many_pop_significance_modified(folders, colors, ax, special_dict = None, summary_path = '/home/user/DATA/Astrid/run_rest_summary') : 
    
    for i,folder in enumerate(folders) : 
        for j,key in enumerate(['%s_angle-0.0' % folder, 
            '%s_angle-90.0' % folder]) :

            if special_dict is not None: 
                Sensitivities = np.load(summary_path + '/' + special_dict['name'] + '_Sensitivities_%s.npy' % key, allow_pickle=True)  
                original_Sensitivities = np.load(summary_path + '/Sensitivities_%s.npy' % key, allow_pickle=True) 
                ylabel = special_dict['ylabel']
                title = special_dict['name'] +  ' \n'
                ax.set_ylim(special_dict['ylims'])

                Responses = [np.mean(S['Responses'][np.sum(oS['significant_pos'] + oS['significant_neg'], axis = 1).astype(bool)], axis=0) for S, oS in zip(Sensitivities, original_Sensitivities)]
                n_cells = np.sum([np.sum(np.sum(s['significant_pos'] + s['significant_neg'], axis = 1).astype(bool)) for s in original_Sensitivities])
                #subtitle = 'Significant cells only'
                subtitle = ''

            uncertainty_sy = session_sem_with_indepedance_hypothesis(Sensitivities)
            
            pt.plot(Sensitivities[0]['contrast'], 
                np.nanmean(Responses, axis=0), 
                sy=uncertainty_sy,
                color=colors[i][j],
                ax=ax)
            
            ax.annotate(text = 'N= ' + str(n_cells), xy = (0.8,0.05-(0.05*j)-(0.1*i)), color = colors[i][j], fontsize = 13)
    ax.set_xticks(ticks = np.arange(3)*0.5)
    ax.set_ylabel(ylabel, fontsize = 13)
    ax.set_xlabel('contrast', fontsize = 13)
    if special_dict is not None :
        if 'title' in special_dict :  
            title = special_dict['title']
            subtitle = ''
    ax.set_title(label=title + subtitle, fontsize = 13)

