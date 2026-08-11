 # Build All/Run/Rest Tuning summaries across conditions
#%%
import os, sys , shutil 
import multiprocessing

sys.path += ['/home/user/lab-notebook/astrid/physion/src']
sys.path += ['/home/user/lab-notebook/astrid']
from run_rest_responses.tuning_arousal_summary_functions import *

#%% plot


folders = [#"PV-cells_WT_Adult_V1", 
    #"SST-cells_WT_Adult_V1",
    "SST-cells_cond-GluN1-KO_Adult_V1"
    ]

arousal_keys = ['', 'Run_', 'Rest_']
uncertainties = ['session sem with propagation and independance hypothesis']
summary_path = '/home/user/DATA/Astrid/run_rest_summary'
base_path = os.path.expanduser('~/DATA/Astrid/Cibele_data')
ylims = [[-0.05,1],[0.4, 1.05],[-0.1,1.2]]

for i, folder in enumerate(folders):
    colors = [pt.tab10(i), 'lightgrey']

    if 'Young' not in folder  : 
        keys =  ['%s_contrast-1.0' % folder, 
            '%s_contrast-0.5' % folder]
        
        fig,ax = plot_orientation_tuning_curve_with_uncertainty(keys,
                                                                arousal_keys,
                                                                summary_path=summary_path,
                                                                average_by='sessions',
                                                                uncertainties = uncertainties,
                                                                colors = colors, 
                                                                plot_perc_run = False,
                                                                group_ROIs = False,
                                                                gaussian_fit= True,
                                                                base_path = base_path,
                                                                ylims=ylims[i]) 
#%%

plot_orientation_tuning_curve_with_uncertainty_all_pop_in_one_graph(folders,
                              arousal_keys,
                              summary_path=summary_path,
                              average_by='sessions',
                              uncertainties = uncertainties,
                              COLORS = colors, 
                              plot_perc_run = False,
                              group_ROIs = False,
                              gaussian_fit = True,
                              base_path = base_path,
                              YLIMS=[(-0.15,1.25),(-0.15,1.25)])
 
   
#%% Create summaries : Variables to declare before executing

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
"""
Notes on tuning run/rest summary building choices : 

    Overall : for every session, the complete response is always first calculated. Then the overall arousal dependant responses and std values are recalculated using only the arousal dependant responses. 

    - To handle sessions with nans for many variations of orientations, and to ensure ROIs have the same prefered angles within complete and arousal dependant response : 
        -ROIs response values within arousal summaries are rotated according to the ROI's prefered angle in the complete summary and missing values are filled with nans. 
        -For this reason, the ntrials key is no longer a list of dim 8 orientations but an array of shape nROIs * 8 orientations : 
            indeed, a session may have run episodes for let's say 3 orientations. all ROIs response to run will be for these 3 orientations. However, because the response are shifted 
            due to the prefered angles, the responses array may have values for all 8 orientations. hence the ntrials list of every ROIs always contains the same numbers, which are the original 
            number of run episodes. But this ntrials list is also shifted by the prefered angles. 


    -the significance of the ROIs are not recalculated based on their arousal-dependent response. 
        The significances of the complete response is just copy/paste into the arousal summary


Notes on propagated_sem for tuning arousal summaries : 
    the propagation of the variance at the ROI level is done using the summaries ROIs responses hence not normalized. 
    while the propagation of the session variance at the session level is done using the responses from get_tuning_responses which now include normalization and significance condition. 
    Also session for which, an orientations only contains 1 ROI response will have an std of 0 for this orientation. The episodes uncertainty will hence be bigger and passed on. 
    

"""
from physion.assembling.dataset import read_spreadsheet

def process_file(filename, i, c,  arousal_cond = 'Run'):

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
    protocol_name=[p for p in data.protocols if '8orientation' in p][0]
    data.build_dFoF(**dFoF_parameters, verbose=False)

    quantities = ['dFoF']
    if 'Running-Speed' in data.nwbfile.acquisition:
        quantities += ['running']

    if data.nROIs>=nMIN_ROIs:

        try:
            Episodes = EpisodeData(data, 
                                    quantities=['dFoF','running'],
                                    protocol_name=protocol_name, 
                                    verbose=False)
            
            if np.sum(data.running) == 0 : 
                print('File discarded because it has no locomotion values')
                return None 
            
            Tuning = compute_tuning_response_per_cells_with_arousal_cond(data, Episodes,
                                                        arousal_cond = arousal_cond,
                                                        quantity='dFoF', 
                                                        stat_test_props = stat_test_props, 
                                                        response_significance_threshold = response_significance_threshold, 
                                                        contrast =float(c.split('contrast-')[1][:3]),
                                                        nMin_episodes = nMin_episodes,
                                                        start_angle=-22.5, 
                                                        angle_range=180,
                                                        verbose=False)
            
            Tuning['datafile'] = filename
            Tuning['nROIs_original'] = data.original_nROIs
            Tuning['nROIs_final'] = data.nROIs
            
            Tuning['subject'] = data.nwbfile.subject.subject_id
            
            print('      [v] --> included, n=%i ROIs ' % data.nROIs)

            return Tuning

        except BaseException as be:
            print('                        [-------------------------------]')
            print(be)
            print()
            print(filename)
            print('nROIs=%i' % data.nROIs, ', protocols=%s' % data.protocols) 
            print(Episodes.varied_parameters)
            print('      [X] --> discarded, problem in datafile, CHECK [!!]')
            print('                        [-------------------------------]')

    else:
        print('      [X] --> discarded, n=%i ROIs ' % data.nROIs)
    
    
if __name__=='__main__':

    import physion
    
    datasets = {}
    for c in folders:

        for contrast in [0.5, 1.0]:

            datasets[c+'_contrast-%.1f' % contrast] =\
                {'datafolder':os.path.join(base_path, c, 'NWBs'), 
                    'age_interval':None}
            
            # we split young animals into age groups
            if 'Young' in c:
                for interval in AGE_INTERVALS:
                    datasets[c.replace('Young', 'P%i-P%i' % interval)+'_contrast-%.1f' % contrast] =\
                        {'datafolder':os.path.join(base_path, c, 'NWBs'), 
                            'age_interval':interval}
                

    Nstart = 0
    Nend = len(datasets)

    for arousal_cond in arousal_keys : 

        for n in range(Nstart, Nend):

            c = list(datasets.keys())[n]

            table = datasets[c]['datafolder'].replace('NWBs', 'DataTable.xlsx')

            dataset_table, subjects_table, analysis =\
                    read_spreadsheet(table,
                        get_metadata_from='table')
            print()
            print()
            print('=================================================================')
            print('-----------------------------------------------------------------')
            print('------- %i) computing : %s ' % (n, c))
            print('-----------------------------------------------------------------')
            print()

            DATASET = scan_folder_for_NWBfiles(datasets[c]['datafolder'])
            
            # FILTER
            # 1) protocol type: orientation tuning
            cond = np.array([np.sum(['8orientation' in p for p in protocols])\
                            for protocols in DATASET['protocols']], dtype=bool)
            # 2) age condition
            if datasets[c]['age_interval'] is not None: 
                cond = cond &\
                    (DATASET['ages']>=datasets[c]['age_interval'][0]) &\
                    (DATASET['ages']<=datasets[c]['age_interval'][1])


            if len(DATASET['files'][cond])>nMIN_DATAFILES:

                Tunings = []
                for i, f in enumerate(DATASET['files'][cond]):

                    Tuning = process_file(DATASET['files'][cond][i], i, c, arousal_cond)
                    if Tuning != None : 
                        Tunings.append(Tuning)

                # # saving data
                np.save(os.path.join(summary_path, arousal_cond + 'Tunings_%s.npy' % c), Tunings)

            else:
                print()
                print('   [!!]   DATASET NOT LARGE ENOUGH   [!!] ')
                print('               only N=%i sessions available' %\
                                            len(DATASET['files'][cond]))
                print('   [!!]   DATASET not analyzed       [!!] ')
                print()

            print('-----------------------------------------------------------------')
            print('=================================================================')

