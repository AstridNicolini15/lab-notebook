#%%




def plot_bar_mean_raw_corrected_fluo(folders, summary_protocol, protocol_control_cond) : 

    RAW = []
    CORRECTED = []
    F0 = []
    DFOF = []

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

        raw = []
        corrected = []
        f0 = []
        dfof = []
        for ses_number in range(len(summary)) : 
            print(ses_number)

            filename = summary[ses_number]['datafile']
            data = physion.analysis.read_NWB.Data(filename, verbose=False) 
            data.build_dFoF(**dFoF_parameters, verbose=False)

            raw.append(np.mean(data.rawFluo, axis = 1))
            corrected.append(np.mean(data.correctedFluo, axis = 1))
            f0.append(np.mean(data.correctedFluo0, axis = 1))
            dfof.append(np.mean(data.Neuropil.data, axis = 1))

        raw = np.concatenate(raw)
        corrected = np.concatenate(corrected)
        f0 = np.concatenate(f0)
        dfof = np.concatenate(dfof)
    
        RAW.append(raw)
        CORRECTED.append(corrected)
        F0.append(f0)
        DFOF.append(dfof)

    if False : 
        y = [np.mean(RAW[0]), np.mean(RAW[1]), np.mean(CORRECTED[0]), np.mean(CORRECTED[1]), np.mean(F0[0]), np.mean(F0[1])]
        yerr = [stats.sem(RAW[0]), stats.sem(RAW[1]), stats.sem(CORRECTED[0]), stats.sem(CORRECTED[1]), stats.sem(F0[0]), stats.sem(F0[1])]
        #plot 
        plt.figure(figsize = (7,4))
        plt.errorbar(x = [0,1,3,4,6,7], y = y, yerr =  yerr, fmt='.', elinewidth=1, capthick=1,  capsize = 5, color = 'black')
        plt.bar(x = [0,1,3,4,6,7], height = y, width = 0.7, color = [pt.tab10(1),pt.tab10(2)]*3, alpha = 1)
        x1 = 0.7*np.arange(len(RAW[0]))/len(RAW[0])-0.35
        x2 = 0.7*np.arange(len(RAW[1]))/len(RAW[1])-0.35
        #plt.scatter(x = [x2+1, x2+4, x2+7], y = [RAW[1], CORRECTED[1], F0[1]], color = pt.tab10(2))
        #plt.scatter(x = [x1, x1+3,x1+6], y = [RAW[0], CORRECTED[0], F0[0]], color = pt.tab10(1))
        plt.ylabel('F (cells average)')
        #plt.ylim(0,2.1)
        plt.xticks( [0.5,3.5,6.5],  ['rawFluo','correctedFluo', 'correctedFluo0'], rotation = 45)
        plt.text(s=folders[0], x = -1, y = -75, color = pt.tab10(1))
        plt.text(s=folders[1], x = -1, y = -90, color = pt.tab10(2))
        plt.text(s='error bar = sem', x = -1, y = -105)
        plt.text(s=summary_protocol,  x = -1, y = -120)

    x1 = 0.7*np.arange(len(RAW[0]))/len(RAW[0])-0.35
    x2 = 0.7*np.arange(len(RAW[1]))/len(RAW[1])-0.35

    plt.figure(figsize = (3,6))
    plt.errorbar(x = [0,1], y = [np.mean(CORRECTED[0]/F0[0]), np.mean(CORRECTED[1]/F0[1])], yerr =  [stats.sem(CORRECTED[0]/F0[0]), stats.sem(CORRECTED[1]/F0[1])], fmt='.', elinewidth=1, capthick=1,  capsize = 5, color = 'black')
    plt.bar(x = [0,1], height = [np.mean(CORRECTED[0]/F0[0]), np.mean(CORRECTED[1]/F0[1])], width = 0.7, color = [pt.tab10(1),pt.tab10(2)], alpha = 0.7)
    plt.scatter([x1], CORRECTED[0]/F0[0], color = pt.tab10(1))
    plt.scatter([x2+1], CORRECTED[1]/F0[1], color = pt.tab10(2))
    plt.ylabel('$\\Delta$F/F')
    #plt.text(s=folders[0], x = -1, y = -0.3, color = pt.tab10(1))
    #plt.text(s=folders[1], x = -1, y = -0.4, color = pt.tab10(2))
    #plt.text(s='error bar = sem', x = -1, y = -0.5)
    #plt.text(s=summary_protocol,  x = -1, y = -0.6)
    plt.xticks( [0.5],  ['correctedFluo / \n correctedFluo0'], rotation = 45)
    
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



def plot_bar_comparaison_average_run_triggered_dfof(folders) : 


    plt.figure(figsize = (5,7))

    time_from_other_events_cond = 0
    for i,folder in enumerate(folders[1:]) : 
        summary_protocol = 'Tunings'
        protocol_control_cond = 'contrast-1.0'
        summary = np.load('/home/user/DATA/Astrid/run_rest_summary/' + summary_protocol + '_' + folder + '_' + protocol_control_cond + '.npy', allow_pickle=True)  

        summary_portocol, protocol_control_cond = get_summary_info(summary)

        RUN, TRACE, ses_included = get_run_triggered_dfof_and_run(folder, summary, time_from_other_events_cond = time_from_other_events_cond) 
        sign_ROIs = np.concatenate([summary[ses_number]['significant_ROIs'] for ses_number in ses_included])
        TRACE = np.array(np.concatenate(TRACE))
        mean_TRACE = np.mean(TRACE, axis = 0)
        n_cell = len(TRACE)

        
        ind_onset = (10*TRACE.shape[1])//20 - 2
        ind_onset_m2 = (8*TRACE.shape[1])//20
        ind_onset_p2 = (12*TRACE.shape[1])//20

        var_pre_post_sign = [np.mean(tr[ind_onset:ind_onset_p2]) - np.mean(tr[ind_onset_m2:ind_onset]) for tr in TRACE[sign_ROIs]]

        color = pt.tab10(i+1)
        plt.errorbar(x = [i], y = [np.mean(var_pre_post_sign)], yerr =  [sem(var_pre_post_sign)], fmt='.', elinewidth=1, capthick=1,  capsize = 5, color = 'black')
        plt.bar(x = [i], height = [np.mean(var_pre_post_sign)], width = 0.7, color = color, alpha = 0.5)
        plt.scatter(0.7*(np.arange(0,len(var_pre_post_sign))/len(var_pre_post_sign)) - 0.35 + i, var_pre_post_sign, s= 2, color = color)
        plt.text(s='n_cells = '+ str(np.sum(sign_ROIs)) + ' / mean = ' + str(np.round(np.mean(var_pre_post_sign),3)),  x = -0.5 , y = -0.-0.1*i, color = color)
    
    plt.ylabel('Variation $\\Delta$F/F (post-pre)')
    plt.xticks( [0,1],  ['responsive cells \n ' + str(folders[0]),'responsive cells \n ' + str(folders[1])])
    #plt.text(s=folder, x = -0.5, y = -1.1)
    plt.text(s='error bar = sem', x = -0.5, y = -0.6)
    plt.text(s=summary_protocol,  x = -0.5, y = -0.7)
    plt.text(s='time_from_other_events_cond = ' + str(time_from_other_events_cond),  x = -0.5, y = -0.8)
    #plt.text(s='error bar = sem', x = -0.5, y = -1.1)
    #plt.text(s=summary_protocol,  x = -0.5, y = -1.2)
    #plt.text(s='time_from_other_events_cond = ' + str(time_from_other_events_cond),  x = -0.5, y = -1.3)
    #plt.text(s='RESP : n_cells = '+ str(np.sum(sign_ROIs)) + ' / mean = ' + str(np.round(np.mean(var_pre_post_sign),3)),  x = -0.5, y = -1.5)
    #plt.text(s='NON RESP : n_cells = '+ str(np.sum(~sign_ROIs)) + ' / mean = ' + str(np.round(np.mean(var_pre_post_non_sign),3)),  x = -0.5, y = -1.6)
    plt.tight_layout()
    plt.show()

#%%

folders = [#"PV-cells_WT_Adult_V1", 
    "SST-cells_WT_Adult_V1",
    "SST-cells_cond-GluN1-KO_Adult_V1"
    ]
folder = folders[0]
keys =  ['%s_angle-0.0' % folder, 
                        '%s_angle-90.0' % folder]
arousal_cond = ''
key = keys[0]
Sensitivities = \
    np.load(summary_path + '/' + arousal_cond + 'Sensitivities_%s.npy' % key, allow_pickle=True)   

Responses = get_responses(Sensitivities, average_by='sessions')

for i in range(len(Sensitivities)) : 
    plt.figure(figsize = (7,7))
    plt.plot(Responses[i], color = 'black')

    plt.plot(Sensitivities[i]['Responses'].T, color = 'orange')
    plt.show()



for i,folder in enumerate(folders) : 
    plt.figure(figsize=(4,4))
    keys =  ['%s_angle-0.0' % folder, 
                        '%s_angle-90.0' % folder]
    arousal_cond = ''
    key = keys[0]
    Sensitivities = \
        np.load(summary_path + '/' + arousal_cond + 'Sensitivities_%s.npy' % key, allow_pickle=True)   

    Responses = get_responses(Sensitivities, average_by='sessions')

    plt.plot(np.array(Responses).T, color = ['orange','lightgreen'][i])
    plt.ylim(-0.7,1.4)




for i,folder in enumerate(folders) : 
    plt.figure(figsize=(4,4))
    keys =  ['%s_contrast-1.0' % folder, 
            '%s_contrast-0.5' % folder]
        
    arousal_cond = ''
    key = keys[0]
    Tunings = np.load(summary_path + '/' + '' + 'Tunings_%s.npy' % key, allow_pickle=True)   

    #Responses = get_tuning_responses(Tunings, average_by='sessions')
    Responses = [np.nanmean(Tuning['Responses'][Tuning['significant_ROIs'],:],
                        axis=0) for Tuning in Tunings]

    plt.plot(np.array(Responses).T, color = ['orange','lightgreen'][i])
    plt.ylim(-0.7,2.5)

for i,folder in enumerate(folders) : 
    m_resp = []
    plt.figure(figsize = (7,7))
    keys =  ['%s_contrast-1.0' % folder, 
            '%s_contrast-0.5' % folder]
        
    arousal_cond = ''
    key = keys[0]
    Tunings = np.load(summary_path + '/' + '' + 'Tunings_%s.npy' % key, allow_pickle=True)   
    Responses = [np.nanmean(Tuning['Responses'][Tuning['significant_ROIs'],:],
                        axis=0) for Tuning in Tunings]
    
    for j in range(len(Tunings)) : 

        #plt.plot(Responses[j], color = 'black')

        plt.plot(Tunings[j]['Responses'].T, color = ['orange','lightgreen'][i])
        m_resp.append(Tunings[j]['Responses'])
    plt.plot(np.mean(np.concatenate(m_resp),axis=0), color = 'black')
    plt.show()


#%%
