#%% 
#explore PV young 


#plot evoked response by age 



def plot_evoked_response(filename, ax):

    dFoF_parameters = dict(\
        roi_to_neuropil_fluo_inclusion_factor=0,
        neuropil_correction_factor = 1,
        method_for_F0 = 'sliding_percentile',
        percentile=5., # percent
        sliding_window = 5*60, # seconds
        computed_neuropil_fact_method=None
    )

    #build
    data = Data(filename, verbose=False)
    data.build_rawFluo()
    #data.build_neuropil()
    quantities = ['rawFluo']
    #data.build_dFoF(**dFoF_parameters, verbose=False)
    #quantities = ['dFoF']
    Episodes = EpisodeData(data, 
                            quantities=quantities, 
                            protocol_name = data.protocols[0],
                            verbose=False)

    Episodes.dFoF = np.mean(Episodes.rawFluo, axis = 1)

    #compute values of interest
    period_evoked = np.mean(Episodes.dFoF , axis = 0)

    #plot
    ax.plot(Episodes.t, period_evoked, color = 'black')
    ax.set_xlabel('time')
    ax.set_ylabel('mean rawfluo')
    
    return ax 

DATASET = scan_folder_for_NWBfiles('/home/user/DATA/Astrid/Cibele_data/PV-cells_WT_Young_V1/NWBs')

#%%
import string
alphabet_upper = list(string.ascii_uppercase)
for age in np.unique(DATASET['ages']):  
    print("##############################")
    print("age: ", age)
    print("##############################")
    N = len(np.argwhere(DATASET['ages'] == age))
    mosaic = "".join(alphabet_upper[:N])
    fig = plt.figure(layout="constrained", figsize = (4*N,4))
    ax_dict = fig.subplot_mosaic(mosaic)
    fig.suptitle(str(age))
    k=0
    for i in np.argwhere(DATASET['ages'] == age) : 
        filename = DATASET['files'][i][0]
        ax = plot_evoked_response(filename, ax_dict[alphabet_upper[k]])
        ax_dict[alphabet_upper[k]].set_title(filename.split('/')[-1])
        k+=1
    plt.show()

#%%
filename = '/home/user/DATA/Astrid/Cibele_data/SST-cells_WT_Adult_V1/NWBs/2024_08_23-14-01-41.nwb'
mosaic = "A"
fig = plt.figure(layout="constrained", figsize = (4,4))
ax_dict = fig.subplot_mosaic(mosaic)
k=0

plot_evoked_response(filename, ax_dict["A"])
ax_dict['A'].axvspan(0,0.1, color = 'purple', alpha = 0.15)
#%%
dFoF_parameters = dict(\
    roi_to_neuropil_fluo_inclusion_factor=1.15,
    neuropil_correction_factor = 1,
    method_for_F0 = 'sliding_percentile',
    percentile=5., # percent
    sliding_window = 5*60, # seconds
    computed_neuropil_fact_method=None,
    verbose=False,
)

def plot_evoked_response_roi(val, ax, time,color):

    #compute values of interest
    period_evoked = np.mean(val, axis = 0)

    #plot
    ax.plot(time, period_evoked, color = color)
    ax.set_xlabel('time')
    ax.set_ylabel('mean rawfluo - 1*neuropil')
    
    return ax 
#%%

filename = '/home/user/DATA/Astrid/Cibele_data/PV-cells_WT_Young_V1/NWBs/2025_02_27-14-11-29.nwb'
#filename = '/home/user/DATA/Astrid/Cibele_data/PV-cells_WT_Adult_V1/NWBs/2024_04_29-15-10-12.nwb'
filename = '/home/user/DATA/Astrid/Cibele_data/newNWBs-Adult/PYR-cells_WT_Adult_V1/NWBs/2026_07_31-14-58-55.nwb'
data = Data(filename, verbose=False)
#data.build_dFoF(**dFoF_parameters)
data.build_correctedFluo(neuropil_correction_factor = 0.7)
quantities = ['correctedFluo']

Episodes = EpisodeData(data, 
                        quantities=quantities, 
                        protocol_name = data.protocols[0],
                        verbose=False)
N = data.nROIs
#%%
round_nb = 500
cols = 5
rows = 100//5

sign = (np.sum(Sensitivities[1]['significant_pos'], axis = 1) + np.sum(Sensitivities[1]['significant_neg'], axis = 1)).astype(bool)
fig, axes = plt.subplots(rows, cols, figsize=(cols * 3, rows * 3), layout = 'constrained')

axes_flat = axes.flatten()
for roi in range(0+round_nb ,100+round_nb ) : 
    if sign[roi] : 
        val =  Episodes.correctedFluo[:,roi,:] #- 0.7* Episodes.neuropil[:,roi,:]
        plot_evoked_response_roi(val, axes_flat[roi-round_nb], Episodes.t, color = 'green')
        axes_flat[roi - round_nb ].set_title('n° ' + str(roi))

plt.show()

#%%


for i in range(len(Sensitivities)): 
    sign = (np.sum(Sensitivities[i]['significant_pos'], axis = 1) + np.sum(Sensitivities[i]['significant_neg'], axis = 1)).astype(bool)
    responses = Sensitivities[i]['Responses'][sign]
    argmaxs = [np.argmax(r) for r in responses]


    mosaic = "A" 
    fig, ax_dict = create_figure(mosaic, figsize = (3,3))


    ax_dict['A'].hist(argmaxs, bins = np.array([0,1,2,3,4,5,6,7,8]), density = True, rwidth = 0.8, align = 'left')
    ax_dict['A'].set_xticks(np.array([0,1,2,3,4,5,6,7]), labels = np.round(Episodes.varied_parameters['contrast'],2))
    ax_dict['A'].set_ylim(0,0.6)
    ax_dict['A'].set_xlabel('contrast of max response')
    ax_dict['A'].set_ylabel('density' + 'significant cells')
    ax_dict['A'].set_title(Sensitivities[i]['datafile'].split('/')[-1])

#%%
data.build_photodiode()
plot_trace_and_photodiode(data, 149, quantity = None)
#%%

for age in np.unique(DATASET['ages']):  
    evoked_response = []
    for i in np.argwhere(DATASET['ages'] == age) :
        i = i[0]
        if DATASET['protocols'][i][0] == 'ff-gratings-2orientations-8contrasts-15repeats': #'ff-gratings-8orientation-2contrasts-15repeats' : 
            filename = DATASET['files'][i]
            data = Data(filename, verbose=False)
            data.build_rawFluo()
            data.build_neuropil()
            quantities = ['rawFluo', 'neuropil']

            Episodes = EpisodeData(data, 
                                    quantities=quantities, 
                                    protocol_name = data.protocols[0],
                                    verbose=False)
            val = Episodes.rawFluo - 1*Episodes.neuropil
            val = np.mean(val, axis = 0)
            val = np.mean(val, axis = 0)
            evoked_response.append(val)

    mean_evoked_response = np.mean(evoked_response, axis = 0)
    #plot
    fig = plt.figure(figsize = (4,4))
    plt.plot(Episodes.t, mean_evoked_response, color = 'green')
    plt.xlabel('time')
    plt.ylabel('mean rawfluo - 1*neuropil')
    plt.title('age: ' + str(age) + '\n protocol = contrast')
    plt.show()

#%%
summary_paths, folders, colors= get_working_data_defaults(experimenters = ['Cibele'],
                                                        folders_surnames = ['PYR_WT'], 
                                                        color_with_grey=True)
summary_protocol = 'Sensitivities'
summary_filename = '/home/user/DATA/Astrid/Cibele_data/newNWBs-Adult/summary/Zscore_correctedFluo_roi_to_neuropil_fluo_inclusion_factor_0_neuropil_correction_factor_0.7_Sensitivities_PYR-cells_WT_Adult_V1'
session_number = -2 #6
ylims = (-1,1)
plot_cells_response(summary_filename, summary_protocol, session_number, colors, ylims)



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# RESPONSES OF EACH SESSIONS

summary_protocol = 'Sensitivities' #'Tunings'
ylims = [-0.05,0.2] #[0,1.05] 
summary_filename = '/home/user/DATA/Astrid/Cibele_data/newNWBs-Adult/summary/Zscore_correctedFluo_roi_to_neuropil_fluo_inclusion_factor_0_neuropil_correction_factor_0.7_Sensitivities_PYR-cells_WT_Adult_V1'

fig = plot_individual_responses(summary_filename, summary_protocol = summary_protocol, colors = colors[0], ylims = ylims)
