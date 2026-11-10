#%%



summary_paths, folders, colors= get_working_data_defaults(experimenters = ['Cibele'],
                                                        folders_surnames = ['SST_WT','PV_WT','PYR_WT','PYR_WT'])

prefixes = ['', '', '', 'Zscore_correctedFluo_roi_to_neuropil_fluo_inclusion_factor_0_neuropil_correction_factor_0.7_']
summary_filenames = get_summary_names(folders, summary_paths, summary_protocols = ['Sensitivities'])
summary_filenames = [f.split('Sensitivities')[0] + p + f.split('/')[-1] for p, f in zip(prefixes, summary_filenames)]


mosaic = "ABCD" 
fig, ax_dict = create_figure(mosaic, figsize = (14,4))
for summary_filename, color, letter in zip(summary_filenames, colors, ax_dict): 
    plot_contrast_amplitude_max_hist(ax_dict[letter], color, summary_filename, ylims = (0,0.6))

labels = ['SST_WT','PV_WT','PYR_WT','PYR_WT_Zscore_correctedFluo_corrfact07']
add_legend_to_plot(fig, labels, colors,bbox_to_anchor = (0,-0.3), fontsize = 8)  
fig.suptitle('significant pos only')


mosaic = "ABCD" 
fig, ax_dict = create_figure(mosaic, figsize = (14,4))
for summary_filename, color, letter in zip(summary_filenames, colors, ax_dict): 
    plot_cumulative_contrast_amplitude(ax_dict[letter], color, summary_filename, xlims=(-2,1))


#%% plot hist amplitude max resp 


def plot_contrast_amplitude_max_hist(ax, color, summary_filename, ylims = (0,0.6)):

    ARGMAXS = []
    for key in ['_angle-0.0','_angle-90.0'] : 
        Sensitivities =  np.load(summary_filename + key + '.npy', allow_pickle = True)
        contrasts = Sensitivities[0]['contrast']
        for Sensitivity in Sensitivities : 
            #sign = (np.sum(Sensitivity['significant_pos'], axis = 1)  + np.sum(Sensitivity['significant_neg'], axis = 1)).astype(bool)
            sign = (np.sum(Sensitivity['significant_pos'], axis = 1)).astype(bool)
                        
            responses = Sensitivity['Responses'][sign]
            ARGMAXS.append([np.argmax(r) for r in responses])

    ARGMAXS = np.concatenate(ARGMAXS)
    print(ARGMAXS)
    ax.hist(ARGMAXS, bins = np.arange(len(contrasts) +1), density = True, rwidth = 0.8, align = 'left', color = color)
    ax.set_xticks(np.arange(len(contrasts)), labels = np.round(contrasts,2))
    ax.set_ylim(ylims[0], ylims[1])
    ax.set_xlabel('contrast of max response')
    ax.set_ylabel('density significant cells')




def plot_cumulative_contrast_amplitude(ax, color, summary_filename, xlims=(0,1)):

    NORM_RESP = []
    for key in ['_angle-0.0','_angle-90.0'] : 
        Sensitivities =  np.load(summary_filename + key + '.npy', allow_pickle = True)
        for Sensitivity in Sensitivities : 
            #sign = (np.sum(Sensitivity['significant_pos'], axis = 1)  + np.sum(Sensitivity['significant_neg'], axis = 1)).astype(bool)
            sign = (np.sum(Sensitivity['significant_pos'], axis = 1)).astype(bool)
                        
            responses = Sensitivity['Responses'][sign]
            if len(responses) != 0 :
                NORM_RESP.append(np.concatenate([r/max(r) for r in responses]))
        
    
    NORM_RESP = np.concatenate(NORM_RESP)
    values, counts = np.unique(NORM_RESP, return_counts=True)
    probabilites = counts / len(NORM_RESP)
    cumsum = np.cumsum(probabilites)
    ax.plot(values,cumsum, color=color)
    ax.set_xlim(xlims[0], xlims[1])
    ax.set_xlabel('normed amplitude levels')
    ax.set_ylabel('cumulative probability of significant cells')



