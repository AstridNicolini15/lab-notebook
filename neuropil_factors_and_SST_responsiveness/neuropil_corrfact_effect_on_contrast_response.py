#%%

neuropil_inclusion_factors = [1.15]
neuropil_correction_factors = [0,0.15,0.3,0.55,0.7,0.85,1] 


arousal_keys = ['']
folders = ["PV-cells_WT_Adult_V1", 
    "SST-cells_WT_Adult_V1",
    #"SST-cells_cond-GluN1-KO_Adult_V1"
    ]

summary_path = '/home/user/DATA/Astrid/run_rest_summary'
base_path = os.path.expanduser('~/DATA/Astrid/Cibele_data')

nMin_episodes = 2
nMIN_DATAFILES = 2


# age intervals in Young
AGE_INTERVALS = [\
    (15,19), (20,23), (24,27), (16,21), (22,27)]
save_summary = True


#%%

folders = ["PV-cells_WT_Adult_V1", 
    "SST-cells_WT_Adult_V1",
    "SST-cells_cond-GluN1-KO_Adult_V1"
    ]

mosaic = """
    ABC
    """

plot_only_sign = False
fig = plt.figure(layout="constrained", figsize = (18,6))
ax_dict = fig.subplot_mosaic(mosaic,  gridspec_kw={"wspace": 0.2,"hspace": 0.05})

if plot_only_sign :
    fig.suptitle('Signficant cells only', fontsize = 15)
plot_contrast_response_multiple_corrfact_one_graph(folder = folders[0], ax = ax_dict["A"], colors=[pt.tab10(0), 'lightgrey'],
                                neuropil_correction_factors = neuropil_correction_factors, ylims = [0,0.22], plot_only_sign = plot_only_sign)

plot_contrast_response_multiple_corrfact_one_graph(folder = folders[1], ax = ax_dict["B"], colors=[pt.tab10(1), 'lightgrey'],
                                neuropil_correction_factors = neuropil_correction_factors, ylims = [-0.1,0.5], plot_only_sign = plot_only_sign)

plot_contrast_response_multiple_corrfact_one_graph(folder = folders[2], ax = ax_dict["C"], colors=[pt.tab10(2), 'lightgrey'],
                                neuropil_correction_factors = neuropil_correction_factors, ylims = [-0.1,0.5], plot_only_sign = plot_only_sign)



def plot_contrast_response_multiple_corrfact_one_graph(folder,
                                ax,
                                colors=None,
                                neuropil_correction_factors = [],
                                neuropil_inclusion_factor = 1.15,
                                plot_only_sign = False,
                                ylims = None, 
                                summary_path =  '/home/user/DATA/Astrid/run_rest_summary'):

    keys = ['%s_angle-0.0' % folder, 
                '%s_angle-90.0' % folder]
    for i, (key, color) in enumerate(zip(keys, colors)):
        for k,neuropil_correction_factor in enumerate(neuropil_correction_factors) :

            Sensitivities = np.load(summary_path + '/' +  str('corrfact_') + str(neuropil_correction_factor) + 'inclufact_' + str(neuropil_inclusion_factor) + 'Sensitivities_%s.npy' % key, allow_pickle=True)   
            
            if plot_only_sign : 
                Responses = [np.mean(S['Responses'][np.sum(S['significant_pos'] + S['significant_neg'], axis = 1).astype(bool)], axis=0) for S in Sensitivities]
                n_cell = np.sum([np.sum(np.sum(S['significant_pos'] + S['significant_neg'], axis = 1).astype(bool)) for S in Sensitivities])
            else : 
                Responses = [np.mean(S['Responses'], axis=0) for S in Sensitivities]
                n_cell = np.sum([len(S['Responses']) for S in Sensitivities])


            uncertainty_sy = session_sem_with_indepedance_hypothesis(Sensitivities)

            ax.plot(Sensitivities[0]['contrast'], 
                        np.nanmean(Responses, axis=0), 
                        color=color, 
                        alpha =0.3+ k/10, lw = 3)
            if 'PV' in folder : 
                ax.text(x = 0.75 + (0.15*i), y= ylims[1]/3.5-(k*ylims[1]/25), s = 'N =' + str(n_cell), color=color, alpha =0.3+ k/10, fontsize = 9)
            else : 
                ax.text(x = 0.75 + (0.15*i), y= ylims[1]/8-(k*ylims[1]/25), s = 'N =' + str(n_cell), color=color, alpha =0.3+ k/10, fontsize = 9)

    ax.set_ylabel('$\\delta$ $\\Delta$F/F', fontsize = 13)
    ax.set_xlabel('contrast', fontsize = 13)
    ax.set_xticks(np.arange(3)*0.5)
    ax.set_ylim(ylims)

