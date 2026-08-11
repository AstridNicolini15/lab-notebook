import os, sys , shutil 
import multiprocessing

sys.path += ['/home/user/lab-notebook/astrid/physion/src']
sys.path += ['/home/user/lab-notebook/astrid']
from run_rest_responses.contrast_arousal_summary_functions import * 
#%%
#-------------------plot execution----------------------#
folders = [#"PV-cells_WT_Adult_V1", 
    #"SST-cells_WT_Adult_V1",
    "SST-cells_cond-GluN1-KO_Adult_V1"
    ]

ylims = [(-0.05,1),(-0.05,0.45),(-0.15,0.85)] #to automatize 

arousal_keys = ['', 'Run_', 'Rest_']
uncertainties = ['session sem with propagation and independance hypothesis']
group_ROIs = False
group_contrast = False 
summary_path = '/home/user/DATA/Astrid/run_rest_summary'
base_path = os.path.expanduser('~/DATA/Astrid/Cibele_data')

# age intervals in Young
AGE_INTERVALS = [\
    (15,19), (20,23), (24,27), (16,21), (22,27)]

for i, folder in enumerate(folders):

    i=2
    colors = [pt.tab10(i), 'lightgrey']

    if 'Young' not in folder  : 
        keys =  ['%s_angle-0.0' % folder, 
                                '%s_angle-90.0' % folder]
        
        fig, ax =  plot_contrast_sensitivity_with_uncertainty(keys,
                                arousal_keys,
                                summary_path=summary_path,
                                average_by='sessions',
                                uncertainties = uncertainties,
                                colors = colors,
                                plot_perc_run = True,
                                group_ROIs = False,
                                base_path = base_path, 
                                ylims = ylims[i]) 

    else :

        for interval in AGE_INTERVALS:
            keys = []
            for angle in [0.,90.] : 
                keys.append((folder[:-8] + 'P%i-P%i' % interval)+'_V1_angle-%.1f' % angle)

            fig, ax =  plot_contrast_sensitivity_with_uncertainty(keys,
                                    arousal_keys,
                                    summary_path=summary_path,
                                    average_by='sessions',
                                    uncertainties = uncertainties,
                                    colors = colors,
                                    plot_perc_run = True,
                                    group_ROIs = False,
                                    base_path = base_path) 

            i+=1

#%%

fig, ax =  plot_contrast_sensitivity_with_uncertainty_all_pop_in_one_graph(folders,
                        arousal_keys,
                        summary_path=summary_path,
                        average_by='sessions',
                        uncertainties = uncertainties,
                        COLORS = colors,
                        plot_perc_run = True,
                        group_ROIs = False,
                        base_path = base_path, 
                        YLIMS = ylims) 


#%% Variables to declare before executing

arousal_keys = ['', 'Run_', 'Rest_']

arousal_keys = ['']
folders = [#"PV-cells_WT_Adult_V1", 
    #"SST-cells_WT_Adult_V1",
    "SST-cells_cond-GluN1-KO_Adult_V1"
    ]

summary_path = '/home/user/DATA/Astrid/run_rest_summary'
base_path = os.path.expanduser('~/DATA/Astrid/Cibele_data')

nMin_episodes = 2
nMIN_DATAFILES = 2


group_contrast = False 
# age intervals in Young
AGE_INTERVALS = [\
    (15,19), (20,23), (24,27), (16,21), (22,27)]
save_summary = True

#%% 
#------------Build execution----------------#
from physion.assembling.dataset import read_spreadsheet

def process_file(filename, i, c, arousal_cond = 'Run', group_contrast = False):
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
            roi_to_neuropil_fluo_inclusion_factor=1.15,
            neuropil_correction_factor = 0.7,
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
        
            if group_contrast == True : 
                compute_group_contrast(Episodes, grouped_contrast_values = [0.05,1]) 

            if np.sum(data.running) == 0 : 
                print('File discarded because it has no locomotion values')
                return None 
            
            if arousal_cond == 'Run_' : 
                filtering_arousal_cond = compute_arousal_mask(Episodes)[0]

            elif arousal_cond == 'Rest_' : 
                filtering_arousal_cond = compute_arousal_mask(Episodes)[1]

            elif arousal_cond == '' : 
                filtering_arousal_cond = None

            Sensitivity = compute_sensitivity_per_cells(data, Episodes, 
                                                        stat_test_props=stat_test_props, 
                                                        response_significance_threshold = response_significance_threshold, 
                                                        filtering_cond = filtering_arousal_cond,
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
    
    for m in range(len(arousal_keys) ) : 
        arousal_cond = arousal_keys[m]

        for n in range(Nstart, Nend):

            c = list(datasets.keys())[n] #select population (ie folder + angle + age conditions)

            table = datasets[c]['datafolder'].replace('NWBs', 'DataTable.xlsx')

            dataset_table, subjects_table, analysis =\
                    read_spreadsheet(table,get_metadata_from='table')
            print()
            print()
            print('=================================================================')
            print('-----------------------------------------------------------------')
            print('------- %i) computing : %s %a' % (n, c, arousal_cond))
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
                    Sensitivity = process_file(DATASET['files'][cond][i], i, c, arousal_cond, group_contrast)
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


            if save_summary == True : 

                # # saving data
                if group_contrast == True : 
                    np.save(os.path.join(summary_path, arousal_cond + 'grouped_Sensitivities_%s.npy' % c ), 
                        Sensitivities)
                else : 
                    np.save(os.path.join(summary_path, arousal_cond + 'Sensitivities_%s.npy' % c ), 
                            Sensitivities)
                print(c + ' summary saved!')