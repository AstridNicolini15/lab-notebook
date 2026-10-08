#%%
#from scipy.stats import linregress

dFoF_parameters = dict(\
roi_to_neuropil_fluo_inclusion_factor=1.15,
neuropil_correction_factor = 0.7,
method_for_F0 = 'sliding_percentile',
percentile=5., # percent
sliding_window = 5*60, # seconds
)
#for all cell, ie before valid cell filtering (bc valid cells = 1.15 brighter than there neuropil)
#look results if done on rest times only?
#look briefly at contrast as well as tunings => should be the same

summary_path, folders, colors= get_working_data_defaults(experimenters = ['Cibele','Cibele','Taddy', 'Taddy','Cibele'],
                                                        folders_surnames = ['SST_WT', 'SST_GluN1_KO','SST_WT', 'SST_GluN1_KO','PV_WT'])

add_legend = True

mosaic = """
    ABCDE
    """

fig= plt.figure(layout="constrained", figsize = (25,7))
fig.get_layout_engine().set(wspace = 0.15)
ax_dict = fig.subplot_mosaic(mosaic)

ax_dict["A"].scatter([-10],[-10], color = colors[0], label = 'valid cells')
ax_dict["A"].scatter([-10],[-10], color = (1,0,0,1), label = 'invalid cells')
ax_dict["A"].legend(loc="upper left")

if add_legend:
    handles = create_legend(folders, colors)
    fig.legend(handles = handles, loc = 'outside lower left', fontsize = 8)

fig.suptitle('Tunings')

for i, letter in enumerate(['A', 'B', 'C','D','E']) : 
    print(i)
        
    summary = np.load(summary_path[i] + '/' + 'Tunings_' + folders[i] + '_contrast-1.0.npy', allow_pickle = True) #just one cond cause use only to get nwb filenames
    #summary = np.load(summary_path[i] + '/' + 'Sensitivities_' + folders[i] + '_angle-90.0.npy', allow_pickle = True) #just one cond cause use only to get nwb filenames

    ROI_fluo_values = []
    Neuropil_fluo_values =  []
    color_list_discriminate_valids = []

    filenames_list = [s['datafile'] for s in summary]
    for filename in filenames_list : 
        
        data = physion.analysis.read_NWB.Data(filename, verbose=False) 
        data.build_dFoF(**dFoF_parameters, verbose=False)
        color_list_discriminate_valids.append([colors[i] if j in data.valid_roiIndices else (1,0,0,1) for j in range(data.original_nROIs)])


        #run_trace = data.build_running(data.t_dFoF)
        #ind_rest = np.concatenate(np.argwhere(run_trace < 0.5))

        ROI_fluo_values.append(np.mean(data.Fluorescence.data, axis = 0))
        Neuropil_fluo_values.append(np.mean(data.Neuropil.data, axis = 0))

    ROI_fluo_values = np.concatenate(ROI_fluo_values)
    Neuropil_fluo_values = np.concatenate(Neuropil_fluo_values)
    color_list_discriminate_valids = np.concatenate(color_list_discriminate_valids)
    ax_dict[letter].scatter(x = ROI_fluo_values, y = Neuropil_fluo_values, color = color_list_discriminate_valids)
    ax_dict[letter].plot(np.linspace(min(ROI_fluo_values), max(ROI_fluo_values), 100), np.linspace(min(ROI_fluo_values), max(ROI_fluo_values), 100), color = 'black')
    ax_dict[letter].annotate(text = 'x = y', xy = (max(ROI_fluo_values), 0.95*max(ROI_fluo_values)),color = 'black')
    ax_dict[letter].set_xlabel('mean ROI fluorescence')
    ax_dict[letter].set_ylabel('mean Neuropil fluorescence')
    ax_dict[letter].annotate(text = 'N = ' +str(len(ROI_fluo_values)), xy = (max(ROI_fluo_values), max(ROI_fluo_values)), color = colors[i])
    ax_dict[letter].set_ylim(-0.05*min(Neuropil_fluo_values), max(ROI_fluo_values))

