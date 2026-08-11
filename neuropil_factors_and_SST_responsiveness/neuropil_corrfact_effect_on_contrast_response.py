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

#------------Build execution----------------#
from physion.assembling.dataset import read_spreadsheet

def process_file(filename, i, c, neuropil_inclusion_factor, neuropil_correction_factor):
    """
    Modified to add : 
    -arousal condition : nb of summary created = nb arousal conds * nb keys.
    -2 new keys for stats purposes compared to the original summaries : see compute sensitivity per cell.
    -correcting missing contrast responses and stds with nans : needed for sessions having no episode that match the angle*contrast*arousal cond.
        in every summary, every session, we always keep the original structure nb contrast * nb ROIs, and add nans in case of no episode matching the conditions.
    -group by contrast option : replace the 8 contrast values by 2, hence artificially grouping contrasts. not thought to be used a lot hence group values not passed as process_file arguments
    -check for locomotion issue : some nwbs have locomotion metadata set to True but no values. 
    """
    # to be a valid datafile:
    nMIN_ROIs = 4
    # calcium pre-processing params
    dFoF_parameters = dict(\
            roi_to_neuropil_fluo_inclusion_factor=neuropil_inclusion_factor,
            neuropil_correction_factor = neuropil_correction_factor,
            method_for_F0 = 'sliding_percentile',
            percentile=5., # percent
            sliding_window = 5*60, # seconds
    )
    # statistical test for visually-evoked-responses
    stat_test_props=dict(interval_pre=[-1.,0],
                         interval_post=[1.,2.],                                   
                         test='ttest',                                            
                         sign='positive')

    response_significance_threshold=5e-2

    print('%i) ' % (i+1), 'analyzing file: %s  [...] ' % filename)
    data = Data(filename, verbose=False)
    protocol_name=[p for p in data.protocols if '8contrast' in p][0]
    data.build_dFoF(**dFoF_parameters, verbose=False)


    if data.nROIs>=nMIN_ROIs:

        try:
            Episodes = EpisodeData(data, 
                                    quantities=['dFoF','running'], 
                                    protocol_name=protocol_name, 
                                    verbose=False)

            if np.sum(data.running) == 0 : 
                print('File discarded because it has no locomotion values')
                return None 
            
            Sensitivity = compute_sensitivity_per_cells(data, Episodes, 
                                                        stat_test_props=stat_test_props, 
                                                        response_significance_threshold = response_significance_threshold, 
                                                        filtering_cond = None,
                                                        quantity='dFoF', 
                                                        angle = float(c.split('angle-')[1][:3]),
                                                        nMin_episodes = nMin_episodes)
            
            Sensitivity = correct_missing_responses_and_stds_with_nans(Sensitivity, Episodes)
            Sensitivity['datafile'] = filename
            Sensitivity['nROIs_original'] = data.original_nROIs
            Sensitivity['nROIs_final'] = data.nROIs
            Sensitivity['subject'] = data.nwbfile.subject.subject_id

            print('      [v] --> included, n=%i ROIs ' % data.nROIs)
        except BaseException as be:
            print('                        [-------------------------------]')
            print(be)
            print()
            print('      [X] --> discarded, problem in datafile, CHECK [!!]')
            print('                        [-------------------------------]')
 
    else:
        print('      [X] --> discarded, n=%i ROIs ' % data.nROIs)
        Sensitivity = None
    
    return Sensitivity

for neuropil_correction_factor in neuropil_correction_factors : 
    print('COMPUTING CORRECTION FACTOR ' + str(neuropil_correction_factor))
    for neuropil_inclusion_factor in neuropil_inclusion_factors : 
        print('COMPUTING INCLUSION FACTOR ' + str(neuropil_inclusion_factor))
        if __name__=='__main__':

            import physion

            datasets = {}
            for c in folders:

                for angle in [0., 90.]:

                    datasets[c+'_angle-%.1f' % angle] =\
                        {'datafolder':os.path.join(base_path, c, 'NWBs'), 
                            'age_interval':None}
                    
                    # we split young animals into age groups
                    if 'Young' in c:
                        for interval in AGE_INTERVALS:
                            datasets[c.replace('Young', 'P%i-P%i' % interval)+'_angle-%.1f' % angle] =\
                                {'datafolder':os.path.join(base_path, c, 'NWBs'), 
                                    'age_interval':interval}

            Nstart = 0
            Nend = len(datasets) 
            
            for n in range(Nstart, Nend):

                c = list(datasets.keys())[n] #select population (ie folder + angle + age conditions)

                table = datasets[c]['datafolder'].replace('NWBs', 'DataTable.xlsx')

                dataset_table, subjects_table, analysis =\
                        read_spreadsheet(table,get_metadata_from='table')
                print()
                print()
                print('=================================================================')
                print('-----------------------------------------------------------------')
                print('------- %i) computing : %s %a' % (n, c, ' '))
                print('-----------------------------------------------------------------')
                print()

                DATASET = scan_folder_for_NWBfiles(datasets[c]['datafolder'])
                
                # FILTER

                # 1) protocol type: contrast sensitivity
                cond = np.array([np.sum(['8contrast' in p for p in protocols])\
                                for protocols in DATASET['protocols']], dtype=bool)
                
                # 2) age condition
                if datasets[c]['age_interval'] is not None: #in all the nwb files found only keep those within the age interval of this c population
                    cond = cond &\
                        (DATASET['ages']>=datasets[c]['age_interval'][0]) &\
                        (DATASET['ages']<=datasets[c]['age_interval'][1])


                if len(DATASET['files'][cond])>nMIN_DATAFILES:

                    Sensitivities = []
                    for i, f in enumerate(DATASET['files'][cond]):
                        print(DATASET['files'][cond][i])
                        Sensitivity = process_file(DATASET['files'][cond][i], i, c, neuropil_inclusion_factor, neuropil_correction_factor)
                        if Sensitivity != None :  
                            Sensitivities.append(Sensitivity)

                else:
                    print()
                    print('   [!!]   DATASET NOT LARGE ENOUGH   [!!] ')
                    print('               only N=%i sessions available' %\
                                                len(DATASET['files'][cond]))
                    print('   [!!]   DATASET not analyzed       [!!] ')
                    print()

                print('-----------------------------------------------------------------')
                print('=================================================================') 

                np.save(os.path.join(summary_path, str('corrfact_') + str(neuropil_correction_factor) + 'inclufact_' + str(neuropil_inclusion_factor) + 'Sensitivities_%s.npy' % c ), 
                        Sensitivities)
                print(c + ' summary saved!')



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

