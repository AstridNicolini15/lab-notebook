#%%
from common_fcts import *
#%%
folders = [#"PV-cells_WT_Adult_V1", 
    "SST-cells_WT_Adult_V1",
    "SST-cells_cond-GluN1-KO_Adult_V1"
    ]

mosaic = """
    AAAAAABB
    AAAAAABB
    AAAAAACC
    AAAAAACC
    EEHHHHHH
    EEFFFGGG
    EEFFFGGG
    EEFFFGGG
    """

fig = plt.figure(layout="constrained", figsize = (10,11))
ax_dict = fig.subplot_mosaic(mosaic)
ax_dict["H"].axis("off")

#plot tuning resp
colors = [[pt.tab10(1),'lightgrey'], [pt.tab10(2),'lightgrey']]
plot_tuning_responses_many_pop(folders, colors, ax_dict["A"]) 


#plot resp to visual pies
colors_list = [[pt.tab10(1), 'white'] , [pt.tab10(2),'white']]
plot_responsiveness_pie(folders, colors_list, [ax_dict["B"], ax_dict["C"]])


colors = [pt.tab10(1), pt.tab10(2)]
#plot_mean_F_val(folders, 'correctedFluo0', colors, ax_dict["E"])

colors = [[pt.tab10(1),'lightgrey'], [pt.tab10(2),'lightgrey']]
for neuropil_inclusion_factor, ax in zip([2,3], [ax_dict["F"], ax_dict["G"]]) : 
    special_dict = {'name' : 'withstd_corrfact_0.7inclufact_' + str(neuropil_inclusion_factor), 
                    'title' : 'neuropil inclusion \n factor = '+ str(neuropil_inclusion_factor),
                    'ylims' : None}
    
    plot_tuning_responses_many_pop(folders, colors, ax, special_dict = special_dict, summary_path = '/home/user/DATA/Astrid/summary_neuropil_factor')

#%%
#---------TUNING RESPS----------#
def plot_tuning_responses_many_pop(folders, 
                                   colors, 
                                   ax, 
                                   special_dict = {}, 
                                   summary_path = '/home/user/DATA/Astrid/run_rest_summary') : 
    

    ylims = None
    if bool(special_dict): 
        ax.set_title(special_dict['title'], fontsize = 13)
        ylims = special_dict['ylims']

    for i,folder in enumerate(folders) : 

        keys =  ['%s_contrast-1.0' % folder, 
                '%s_contrast-0.5' % folder]
        
        for j, key in enumerate(keys) :

            xy = (100,0.8-(0.06*j)-(0.12*i))
            draw_tuning_curve(key, 
                    special_dict = special_dict, 
                    summary_path = summary_path,
                    ax = ax, 
                    color = colors[i][j],
                    alpha = 1,
                    xy = xy, 
                    draw_uncertainty = True,
                    graph_width_dict = {'lw' : 3, 'ms' : 5, 'fontsize' : 10}) 


    add_axis_labels_to_plot(ax = ax)
    ax.set_ylim(ylims)


def get_defaults() : 
    summary_path = '/home/user/DATA/Astrid/run_rest_summary'
    return summary_path 

#%%
#---------RESPONSIVENESS TO VISUAL STIM----------#
     
def get_responsiveness_to_visual_stim(folder, summary_protocol = 'Tunings') : 

    perc_resp_to_visual = ()
    if summary_protocol == 'Tunings' :
        for key in ['%s_contrast-1.0' % folder, 
            '%s_contrast-0.5' % folder] : 

            Tunings = np.load('/home/user/DATA/Astrid/run_rest_summary/' + summary_protocol + '_' + key + '.npy', allow_pickle=True) 

            n_cell_resp = np.sum([np.sum(Tuning['significant_ROIs']) for Tuning in Tunings])
            n_cell_nonresp = np.sum([np.sum(~Tuning['significant_ROIs']) for Tuning in Tunings])

            perc_resp_to_visual = (*perc_resp_to_visual, (n_cell_resp *100) / (n_cell_resp + n_cell_nonresp))
            perc_resp_to_visual = (*perc_resp_to_visual, (n_cell_nonresp *100) / (n_cell_resp + n_cell_nonresp))


    if summary_protocol == 'Sensitivities' :
        for key in ['%s_angle-90.0' % folder, 
            '%s_angle-0.0' % folder] : 

            Sensitivities = np.load('/home/user/DATA/Astrid/run_rest_summary/' + summary_protocol + '_' + key + '.npy', allow_pickle=True) 
            cell_responsiveness = np.concatenate([np.sum(S['significant_pos'] + S['significant_neg'], axis = 1) for S in Sensitivities])

            n_cell_resp = len([x for x in cell_responsiveness if x !=0])
            n_cell_nonresp = len([x for x in cell_responsiveness if x ==0])

            perc_resp_to_visual = (*perc_resp_to_visual, (n_cell_resp *100) / (n_cell_resp + n_cell_nonresp))
            perc_resp_to_visual = (*perc_resp_to_visual, (n_cell_nonresp *100) / (n_cell_resp + n_cell_nonresp))

    return  perc_resp_to_visual


def plot_responsiveness_pie(folders, colors_list, axes, summary_protocol = 'Tunings') : 

    for i,folder in enumerate(folders) : 
        perc_resp_c1, perc_nonresp_c1, perc_resp_c05, perc_nonresp_c05 = np.round(get_responsiveness_to_visual_stim(folder, summary_protocol = summary_protocol),3)
        wedges, texts = axes[i].pie([perc_resp_c1, perc_nonresp_c1], 
                        colors = colors_list[i], 
                        startangle = 90,
                        wedgeprops={"edgecolor":"black",'linewidth': 1, 'width' : 0.5},
                        textprops = {"fontsize" : 10})
        axes[i].annotate(text = f'{perc_resp_c1:.1f}% \n resp', xy= (-1.2,1), color = colors_list[i][0], fontsize = 13)
    #ax.set_axis_off()

#%%
#---------MEAN RAW FLUO BAR----------#

def get_unprocessed_F(folder, F_value = 'correctedFluo0', summary_protocol = 'Tunings') : 

    print("computing RawFluo")
    
    if summary_protocol == 'Tunings' :
        summary = np.load('/home/user/DATA/Astrid/run_rest_summary/' + summary_protocol + '_' + folder + '_contrast-1.0.npy', allow_pickle=True)  #only one control cond is sufficient

    elif summary_protocol == 'Sensitivities' :
        summary = np.load('/home/user/DATA/Astrid/run_rest_summary/' + summary_protocol + '_' + folder + '_angle-90.0.npy', allow_pickle=True)  #only one control cond is sufficient

    dFoF_parameters = dict(\
            roi_to_neuropil_fluo_inclusion_factor=1.15,
            neuropil_correction_factor = 0.7,
            method_for_F0 = 'sliding_percentile',
            percentile=5., # percent
            sliding_window = 5*60, # seconds
            with_correctedFluo_and_F0=True,
    )

    #load summary to get nwb filenames
    Fval = []

    for ses_number in range(len(summary)) : 

        filename = summary[ses_number]['datafile']
        data = physion.analysis.read_NWB.Data(filename, verbose=False) 
        data.build_dFoF(**dFoF_parameters, verbose=False)

        if F_value == 'correctedFluo0' : 
            Fval.append(np.mean(data.correctedFluo0, axis = 1))
        elif F_value == 'dFoF' : 
            Fval.append(np.mean(data.dFoF, axis = 1))

    Fval = np.concatenate(Fval)

    return Fval

def plot_mean_F_val(folders, F_value, colors, ax, summary_protocol = 'Tunings') :

    x = np.arange(0,len(folders))
    y = []
    y_err = []
    n_cells_pop = []

    for i,folder in enumerate(folders) : 

        Fval = get_unprocessed_F(folder, F_value = F_value, summary_protocol = summary_protocol)
        #Fval = [1,2,3]
        y.append(np.mean(Fval))
        y_err.append(stats.sem(Fval))
        n_cells_pop.append(len(Fval))

    ax.errorbar(x = x, y = y, yerr = y_err, fmt='.', elinewidth=1, capthick=1,  capsize = 5, color = 'black')
    ax.bar(x = x, height = y, width = 0.7, color = colors)
    ax.set_ylabel(F_value, fontsize = 13)
    ax.set_xticks([])

    ymax = ax.get_ylim()[1]
    for i in range(len(folders)) : 
        ax.text(s=' N=\n'+str(n_cells_pop[i]),  x = x[i] - 0.1, y = ymax/22, rotation = 'horizontal', fontsize = 10)







# %%
