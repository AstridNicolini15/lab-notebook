#%% 
from scipy.stats import sem
'PV and SST orientation tuning comparaison'

summary_path, folders, colors= get_working_data_defaults(experimenters = 'Cibele',
                                                        folders_surnames = ['PV_WT', 'SST_WT'],
                                                        color_with_grey = True)


mosaic = """

    AAAAAABBBBBBCCCCC
    AAAAAABBBBBBCCCCC
    AAAAAABBBBBBCCCCC
    AAAAAABBBBBBCCCCC
    AAAAAABBBBBBCCCCC
    aaaaaabbbbbbCCCCC
    ccccccddddddCCCCC
    zzzzzzzzzzzzyyyyy
    EEEEEEFFFFFFyyyyy
    EEEEEEFFFFFFyyyyy
    EEEEEEFFFFFFyyyyy
    EEEEEEFFFFFFyyyyy
    EEEEEEFFFFFFyyyyy
    eeeeeeffffffyyyyy
    gggggghhhhhhyyyyy
    """



fig = plt.figure(layout="tight", figsize = (17,16))
#fig.get_layout_engine().set(wspace = 0.01)
ax_dict = fig.subplot_mosaic(mosaic)
ax_dict["z"].axis("off")
ax_dict["y"].axis("off")

neuropil_correction_factors = [0,0.15,0.3,0.55,0.7,0.85,1]
structure_letters = [[['A','a','c'],['B','b','d']],[['E','e','g'],['F','f','h']]]

for i, folder in enumerate(folders) : 

    for j, inclufact in enumerate([1.15,4]) : 

        letters = structure_letters[i][j]
        special_dict = {'base_name' : 'inclusion_factor-' + str(inclufact) + '_correction_factor-'}

        plot_orientation_tuning_multiple_corrfact(folder,
                            colors=colors[i],
                            neuropil_correction_factors = neuropil_correction_factors,
                            ax = ax_dict[letters[0]],
                            special_dict = special_dict,
                            summary_path = summary_path)
        ax_dict[letters[0]].set_title('neuropil inclusion factor = ' + str(inclufact), fontsize = 13)
        ax_dict[letters[1]].set_title('neuropil correction factors', fontsize = 13)
        plot_corrfact_legend_color_maps(colors = colors[i], fig = fig,  axes = [ax_dict[letters[1]], ax_dict[letters[2]]], neuropil_correction_factors = neuropil_correction_factors) 


    
height = []
height_error = []
gOSIs_list = []
for i, folder in enumerate(folders) : 
        
    for key in ['%s_contrast-1.0' % folder,
                '%s_contrast-0.5' % folder]:
        
        for inclusion_factor in [1.15,4] : 
        
            summary_prefixe_name = 'inclusion_factor-' + str(inclusion_factor) + '_correction_factor-0.7_'
            Tunings = np.load(summary_path[i] + '/' + summary_prefixe_name + 'Tunings_%s.npy' % key, allow_pickle=True)
            gOSIs = get_gOSIs(Tunings)
            gOSIs_list.append(np.concatenate(gOSIs))
            height.append(np.mean([np.mean(session_gOSIs) for session_gOSIs in gOSIs])) #average of sessions averages
            height_error.append(sem([np.mean(session_gOSIs) for session_gOSIs in gOSIs]))
            

x = [0,1,3,4,7,8,10,11]
ax_dict["C"].grid(True, linewidth = 1, axis = 'y')
colors = reshape_colors_with(colors, second_cond_low_alpha = True, second_cond_alpha = 0.3)
colors = np.concatenate([[colors[0][0]]*2, [colors[0][1]]*2, [colors[1][0]]*2, [colors[1][1]]*2])
ax_dict["C"].bar(x = x, height = height, width = 0.75, color = colors)
ax_dict["C"].errorbar(x=x, y= height, yerr=height_error, fmt='.', elinewidth=1, capthick=1,  capsize = 5, color = 'black')
ax_dict["C"].set_xticks(ticks = x, labels = [1.15,4]*4, rotation = 90, fontsize = 10)
ax_dict["C"].set_xlabel('Neuropil inclusion factor')
ax_dict["C"].set_ylabel('gOSI')
#ax_dict["D"].set_ylim(0,0.45)
ax_dict["C"].set_title('negative clipped')


ax_dict['y'].annotate(text = 'Mann-Whitney U tests : ', xy = (0,0.9), fontsize = 13)

r, pval = stats.mannwhitneyu(x = gOSIs_list[0], y = gOSIs_list[1])
ax_dict['y'].annotate(text = 'PV contrast 1 inclusion factor 1.15 vs 4 : \n           ' + str(np.round(pval,6)), xy = (0,0.8), fontsize = 13)
ax_dict['y'].annotate(text = '***', xy = (0,0.8), fontsize = 13)
r, pval = stats.mannwhitneyu(x = gOSIs_list[2], y = gOSIs_list[3])
ax_dict['y'].annotate(text = 'PV contrast 0.5 inclusion factor 1.15 vs 4 : \n           ' + str(np.round(pval,6)), xy = (0,0.7), fontsize = 13)

r, pval = stats.mannwhitneyu(x = gOSIs_list[4], y = gOSIs_list[5])
ax_dict['y'].annotate(text = 'SST contrast 1 inclusion factor 1.15 vs 4 : \n           ' + str(np.round(pval,6)), xy = (0,0.6), fontsize = 13)
ax_dict['y'].annotate(text = '*', xy = (0,0.6), fontsize = 13)
r, pval = stats.mannwhitneyu(x = gOSIs_list[6], y = gOSIs_list[7])
ax_dict['y'].annotate(text = 'SST contrast 0.5 inclusion factor 1.15 vs 4 : \n           ' + str(np.round(pval,6)), xy = (0,0.5), fontsize = 13)

r, pval = stats.mannwhitneyu(x = gOSIs_list[5], y = gOSIs_list[7])
ax_dict['y'].annotate(text = 'SST inclusion factor 4 contrast 1 vs 0.5 : \n           ' + str(np.round(pval,6)), xy = (0,0.4), fontsize = 13)
ax_dict['y'].annotate(text = '***', xy = (0,0.4), fontsize = 13)

r, pval = stats.mannwhitneyu(x = gOSIs_list[1], y = gOSIs_list[5])
ax_dict['y'].annotate(text = 'inclusion factor 4 contrast 1 PV vs SST : \n           ' + str(np.round(pval,6)), xy = (0,0.3), fontsize = 13)

r, pval = stats.mannwhitneyu(x = gOSIs_list[3], y = gOSIs_list[7])
ax_dict['y'].annotate(text = 'inclusion factor 4 contrast 0.5 PV vs SST : \n           ' + str(np.round(pval,6)), xy = (0,0.2), fontsize = 13)
ax_dict['y'].annotate(text = '***', xy = (0,0.2), fontsize = 13)


r, pval = stats.mannwhitneyu(x = gOSIs_list[0], y = gOSIs_list[4])
ax_dict['y'].annotate(text = 'inclusion factor 1.15 contrast 1 PV vs SST : \n           ' + str(np.round(pval,6)), xy = (0,0.1), fontsize = 13)
ax_dict['y'].annotate(text = '***', xy = (0,0.1), fontsize = 13)

r, pval = stats.mannwhitneyu(x = gOSIs_list[2], y = gOSIs_list[6])
ax_dict['y'].annotate(text = 'inclusion factor 1.15 contrast 0.5 PV vs SST : \n           ' + str(np.round(pval,6)), xy = (0,0), fontsize = 13)
ax_dict['y'].annotate(text = '***', xy = (0,0), fontsize = 13)