#%%
n_sign_list = ['','sign3dir_','sign5dir_','sign7dir_']

plot_orientation_tuning_curve_with_uncertainty_all_pop_in_one_graph(folders,
                              n_sign_list,
                              summary_path=summary_path,
                              average_by='sessions',
                              uncertainties = ['std'],
                              COLORS = [[pt.tab10(1),'lightgrey'],[pt.tab10(2),'lightgrey']], 
                              plot_perc_run = False,
                              group_ROIs = False,
                              gaussian_fit = True,
                              base_path = base_path,
                              YLIMS=[[-0.05,1.05],[-0.05,1.05]]) 




c12 = list(pt.tab10(1))
c12[-1] = 0.35
c12 = tuple(c12)
c22 = list(pt.tab10(2))
c22[-1] = 0.35 
c22 = tuple(c22)



plot_orientation_tuning_curve_diff_significatif_many_pop(folders,
                              n_sign_list,
                              summary_path=summary_path,
                              average_by='sessions',
                              uncertainties = ['sem'],
                              COLORS = [[pt.tab10(1),c12],[pt.tab10(2),c22]], 
                              plot_perc_run = False,
                              group_ROIs = False,
                              gaussian_fit = True,
                              base_path = base_path,
                              YLIMS=[[-0.05,1.1],[-0.05,1.1]]) 
    
#%%
import matplotlib.pyplot as plt
def plot_orientation_tuning_curve_diff_significatif_many_pop(folders,
                              n_sign_list,
                              summary_path='',
                              average_by='sessions',
                              uncertainties = ['std'],
                              COLORS = None, 
                              plot_perc_run = True,
                              group_ROIs = False,
                              gaussian_fit = True,
                              base_path = '',
                              YLIMS=[[0,1]]) : 
    
    fig, axes = plt.subplots(ncols=4, nrows=1, figsize=(6,4))
    plt.subplots_adjust(right=2.5, hspace =0)
    for i,folder in enumerate(folders) : 
        print(i)
        keys =  ['%s_contrast-1.0' % folder, 
            '%s_contrast-0.5' % folder]
        ylims = YLIMS[i]
        colors = COLORS[i]
        folder = keys[0][:-13]

        for uncertainty in uncertainties : 

            if type(keys)==str:
                keys, colors = [keys], [colors[0]]

            if plot_perc_run == True : 
                perc_ep = compute_perc_ep_run(folder, summary_path)
                fig.text(x=0, y=-0.2, s='%' + ' of episodes considered runned : ' + str(perc_ep))


            fig.text(x=0.3, y=1.2-(0.1*i), s=folder, color = colors[0])
            fig.text(x=0.3+(0.3*i), y=-0.1, s= 'contrast 1.0', color = colors[0])
            fig.text(x=0.45+(0.3*i), y=-0.1, s= 'contrast 0.5', color = colors[1])
            fig.text(x=0.3, y=-0.2, s= 'error bars = ' + uncertainty)
            fig.text(x=0.3, y=0.15, s='n significant cells = ', color = 'black')
               


            x = np.linspace(-30, 180-30, 100)
            for k, prefixe in enumerate(n_sign_list):
                m = 0 
                for key, color in zip(keys, colors):

                    # load data
                    Tunings = np.load(summary_path + '/' + prefixe + 'Tunings_%s.npy' % key, allow_pickle=True)   

                    Responses = get_tuning_responses(Tunings, average_by=average_by) #is normed
                    # Gaussian Fit
                    C, func = fit_gaussian(Tunings[0]['shifted_angle'],
                                            np.nanmean(Responses, axis=0))

                    axes[k].plot(x, func(x), lw=2, color=color)

                    if uncertainty == 'std' : 
                        uncertainty_sy = np.nanstd(Responses, axis=0, ddof = 1)

                    elif uncertainty == 'sem' : 
                        uncertainty_sy = stats.sem(Responses, axis=0, nan_policy = 'omit', ddof = 1) 

                    elif uncertainty == 'session sem with propagation and independance hypothesis' : 
                        uncertainty_sy = session_sem_with_indepedance_hypothesis_universal(Tunings)

                    axes[k].scatter(Tunings[0]['shifted_angle'], np.nanmean(Responses, axis=0), color=color)
                    axes[k].errorbar(x=Tunings[0]['shifted_angle'],y = np.nanmean(Responses, axis=0), yerr = uncertainty_sy, fmt='.',color=color)
                    n_sign_cell = np.sum([np.sum(t['significant_ROIs']) for t in Tunings])
                    axes[k].text(x=20+(50*i), y=-0.4-m*0.1, s=str(n_sign_cell), color = color,fontsize = 12)
                    m+=1
                axes[k].set_xticks(ticks = Tunings[0]['shifted_angle'], labels = ['%i' % a if (a in [0, 90]) else '' for a in Tunings[0]['shifted_angle']])
                axes[k].set_ylabel('norm. $\\delta$ $\\Delta$F/F')
                axes[k].set_ylim([-0.05,1.05])
                axes[k].set_xlabel('angle ($^o$) from pref.')
                axes[k].set_title(['significant in 1 dir','significant in 3 dir','significant in 5 dir','significant in 7 dir'][k])
     

    return fig, axes


#%%

values_of_interest = ['Fluorescence','Neuropil']

def plot_bar_mean_many_pop(values_of_interest, folders, summary_protocol, protocol_control_cond) : 

    VALUES = np.empty((len(values_of_interest), 0)).tolist()

    for i, folder in enumerate(folders) :

        #load summary to get nwb filenames
        summary = np.load('/home/user/DATA/Astrid/run_rest_summary/' + summary_protocol + '_' + folder + '_' + protocol_control_cond + '.npy', allow_pickle=True)  


        dFoF_parameters = dict(\
                roi_to_neuropil_fluo_inclusion_factor=1.15,
                neuropil_correction_factor = 0.7,
                method_for_F0 = 'sliding_percentile',
                percentile=5., # percent
                sliding_window = 5*60, # seconds
                with_correctedFluo_and_F0=True,
        )


        values = np.empty((len(values_of_interest), 0)).tolist()

        for ses_number in range(len(summary)) : 
            print(ses_number)

            filename = summary[ses_number]['datafile']
            data = physion.analysis.read_NWB.Data(filename, verbose=False) 
            data.build_dFoF(**dFoF_parameters, verbose=False)

            for j, value_name in enumerate(values_of_interest) : 
                values[j].append(np.mean(data.__dict__[value_name].data, axis = 0))


        values = [np.concatenate(v) for v in values]

        for k in range(len(values)) : 
            VALUES[k].append(values[k])

    y = []
    y_err = []
    x = [0,1,3,4]
    for k in range(len(VALUES)) : 
        for i in range(len(folders)) : 
            y.append(np.mean(VALUES[k][i]))
            y_err.append(stats.sem(VALUES[k][i]))

    print(y_err)
     #plot 
    plt.figure(figsize = (7,4))
    plt.errorbar(x = x, y = y, yerr =  y_err, fmt='.', elinewidth=1, capthick=1,  capsize = 5, color = 'black')
    plt.bar(x = x, height = y, width = 0.7, color = [pt.tab10(1),pt.tab10(2)]*len(values_of_interest), alpha = 0.5)
    x1 = 0.7*np.arange(len(VALUES[0][0]))/len(VALUES[0][0])-0.35
    x2 = 0.7*np.arange(len(VALUES[0][1]))/len(VALUES[0][1])-0.35
    plt.scatter(x = [x1, x1+3], y = [VALUES[0][0], VALUES[1][0]], color = pt.tab10(1))
    plt.scatter(x = [x2+1, x2+4], y = [VALUES[0][1],VALUES[1][1]], color = pt.tab10(2))
    plt.ylim(0,250)
    plt.ylabel('F (cells average)')
    #plt.ylim(0,2.1)
    plt.xticks( [0.5,3.5], values_of_interest, rotation = 45)
    plt.text(s=folders[0], x = -1, y = -75, color = pt.tab10(1))
    plt.text(s=folders[1], x = -1, y = -90, color = pt.tab10(2))
    plt.text(s='error bar = sem', x = -1, y = -105)
    plt.text(s=summary_protocol,  x = -1, y = -120)
    plt.show()



    plt.figure(figsize = (7,4))
    #plt.errorbar(x = x, y = y, yerr =  y_err, fmt='.', elinewidth=1, capthick=1,  capsize = 5, color = 'black')
    plt.bar(x = [0,1], height = [np.mean(VALUES[0][0]/VALUES[1][0]), np.mean(VALUES[0][1]/VALUES[1][1])], width = 0.7, color = [pt.tab10(1),pt.tab10(2)], alpha = 0.5)
    x1 = 0.7*np.arange(len(VALUES[0][0]))/len(VALUES[0][0])-0.35
    x2 = 0.7*np.arange(len(VALUES[0][1]))/len(VALUES[0][1])-0.35
    plt.scatter(x = [x1], y = [VALUES[0][0]/ VALUES[1][0]], color = pt.tab10(1))
    plt.scatter(x = [x2+1], y = [VALUES[0][1]/VALUES[1][1]], color = pt.tab10(2))
    plt.ylim(0,12)
    plt.ylabel('F (cells average)')
    #plt.ylim(0,2.1)
    plt.text(s=folders[0], x = -1, y = -75, color = pt.tab10(1))
    plt.text(s=folders[1], x = -1, y = -90, color = pt.tab10(2))
    plt.text(s='error bar = sem', x = -1, y = -105)
    plt.text(s=summary_protocol,  x = -1, y = -120)
    plt.show()

#%%

    plt.figure(figsize = (3,6))
    plt.errorbar(x = [0,1], y = [np.mean(DFOF[0]), np.mean(DFOF[1])], yerr =  [stats.sem(DFOF[0]), stats.sem(DFOF[1])], fmt='.', elinewidth=1, capthick=1,  capsize = 5, color = 'black')
    plt.bar(x = [0,1], height = [np.mean(DFOF[0]), np.mean(DFOF[1])], width = 0.7, color = [pt.tab10(1),pt.tab10(2)], alpha = 0.7)
    plt.scatter([x1], DFOF[0], color = pt.tab10(1))
    plt.scatter([x2+1], DFOF[1], color = pt.tab10(2))
    plt.ylabel('$\\Delta$F/F')
    #plt.text(s=folders[0], x = -1, y = -0.3, color = pt.tab10(1))
    #plt.text(s=folders[1], x = -1, y = -0.4, color = pt.tab10(2))
    #plt.text(s='error bar = sem', x = -1, y = -0.5)
    #plt.text(s=summary_protocol,  x = -1, y = -0.6)
    plt.xticks( [0.5],  ['dFoF'], rotation = 45)
    
    plt.show()



    plt.figure(figsize = (3,6))
    y0 = CORRECTED[0] / F0[0]
    y1 = CORRECTED[1] / F0[1]
    plt.errorbar(x = [0,1], y = [np.mean(y0), np.mean(y1)], yerr =  [stats.sem(y0), stats.sem(y1)], fmt='.', elinewidth=1, capthick=1,  capsize = 5, color = 'black')
    plt.bar(x = [0,1], height = [np.mean(y0), np.mean(y1)], width = 0.7, color = [pt.tab10(1),pt.tab10(2)], alpha = 0.7)
    plt.scatter([x1], y0, color = pt.tab10(1))
    plt.scatter([x2+1],y1, color = pt.tab10(2))
    plt.ylabel('$\\Delta$F/F')
    #plt.text(s=folders[0], x = -1, y = -0.3, color = pt.tab10(1))
    #plt.text(s=folders[1], x = -1, y = -0.4, color = pt.tab10(2))
    #plt.text(s='error bar = sem', x = -1, y = -0.5)
    #plt.text(s=summary_protocol,  x = -1, y = -0.6)
    plt.xticks( [0.5],  ['correctedFluo / \n correctedFluo0'], rotation = 45)
    
    plt.show()

#%%