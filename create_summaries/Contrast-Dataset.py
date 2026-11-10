# %% [markdown]
# # Build Sensitivity Dataset across Conditions

# %%
import os, sys , shutil 
import multiprocessing
import numpy as np

sys.path += ["/home/user/lab-notebook/astrid/physion/src"]

from physion.analysis.read_NWB\
                         import scan_folder_for_NWBfiles, Data
from physion.analysis.episodes.build import EpisodeData
from physion.analysis.protocols.contrast_sensitivity\
                        import compute_sensitivity_per_cells

from physion.analysis.episodes.build import EpisodeData
from dataset_tools import (get_summary_prefix_name, 
                            get_filtering_cond, 
                            build_filtering_cond_quantities,
                            build_quantity)
from contrast_dataset_tools import *

parallelized, debug = False, False 

# load the dataset locations:
from Dataset_Organization import datasets_func, quantity, summary_folder, filtering_cond_name
datasets = datasets_func('angle', [0.0, 90.])

from Preprocessing_Settings import get_dFoF_params, get_stat_test_props

if sys.argv == ['Contrast-Dataset.py']:
    special_dFoF_params = None
else : 
    special_dFoF_params = sys.argv 
# %%
# to be a valid dataset:
nMIN_DATAFILES = 2

def process_file(filename, i, c, quantity, filtering_cond_name):

    # to be a valid datafile:
    nMIN_ROIs = 4

    # CELL-dependent calcium pre-processing params 
    dFoF_parameters = get_dFoF_params(c, special_dFoF_params)

    # statistical test for visually-evoked-responses
    stat_test_props= get_stat_test_props(c)
    response_significance_threshold=5e-2

    print('%i) ' % (i+1), 'analyzing file: %s  [...] ' % filename)
    data = Data(filename, verbose=False)
    protocol_name=[p for p in data.protocols if '8contrast' in p][0]

    data, quantity, dFoF_parameters, stat_test_props = build_quantity(data, quantity, dFoF_parameters, stat_test_props)

    quantities = [quantity]

    if data.nROIs>=nMIN_ROIs:

        try:
            if filtering_cond_name :
                data, quantities = build_filtering_cond_quantities(filtering_cond_name, data, quantities)

            Episodes = EpisodeData(data, 
                                    quantities=quantities,
                                    protocol_name=protocol_name, 
                                    verbose=False)
            
            filtering_cond = get_filtering_cond(filtering_cond_name, Episodes)

            Sensitivity = compute_sensitivity_per_cells(data, Episodes, 
                                                        quantity=quantity,
                                                        stat_test_props=stat_test_props, 
                                                        response_significance_threshold = response_significance_threshold, 
                                                        filtering_cond = filtering_cond, 
                                                        angle = float(c.split('angle-')[1][:3]))
            if filtering_cond_name :
                Sensitivity = fill_missing_contrast_values_with_nans(Sensitivity, Episodes)

            Sensitivity['datafile'] = filename
            Sensitivity['nROIs_original'] = data.original_nROIs
            Sensitivity['nROIs_final'] = data.nROIs
            Sensitivity['subject'] = data.nwbfile.subject.subject_id

            if  dFoF_parameters['computed_neuropil_fact_method'] : #dFoF_parameters['roi_specific_neuropil_correction_factor'] or dFoF_parameters['with_computed_neuropil_fact'] : 
                Sensitivity['neuropil_correction_factors'] = data.neuropil_correction_factor

            np.save(os.path.join(summary_folder, 'temp', 
                                 'Sensitivity-%s-%s-%i.npy' % (quantity, c, i)),
                    Sensitivity)
            print('      [v] --> included, n=%i ROIs ' % data.nROIs)

        except ValueError as ve: #value error from no locomotion value or no prefered angles for filtered summaries
            print(f"Error: {ve}")
            print('File: %s' % filename, ' discarded')
        
        # except BaseException as be:
        #     print('                        [-------------------------------]')
        #     print(be)
        #     print()
        #     print('      [X] --> discarded, problem in datafile, CHECK [!!]')
        #     print('                        [-------------------------------]')

    else:
        print('      [X] --> discarded, n=%i ROIs ' % data.nROIs)

if __name__=='__main__':

    from physion.assembling.dataset import read_spreadsheet
    cpus = multiprocessing.cpu_count()-1 # leaving 1 cpu for the rest

    # temporary folder for parallelization
    os.makedirs(os.path.join(summary_folder, 'temp'), exist_ok=True)

    Nstart = 0
    Nend = len(datasets) 

    for n in range(Nstart, Nend):

        c = list(datasets.keys())[n]

        #table = datasets[c]['datafolder'].replace('NWBs', 'DataTable.xlsx')

        #dataset_table, subjects_table, analysis =\
        #        read_spreadsheet(table, get_metadata_from='table')
        print()
        print()
        print('=================================================================')
        print('-----------------------------------------------------------------')
        print('------- %i) computing : %s ' % (n+1, c))
        print('-----------------------------------------------------------------')
        print()

        DATASET = scan_folder_for_NWBfiles(datasets[c]['datafolder'])
        # filenames = [\
        #     '2024_08_23-14-01-41',
        #     '2024_08_23-14-26-50',
        #     '2024_08_23-15-19-10',
        #     '2024_08_23-15-44-18',
        #     '2024_08_27-10-27-32',
        #     '2024_08_27-10-52-11',
        #     '2024_08_27-11-16-53',
        #     '2024_08_27-11-45-06',
        #     '2024_08_27-12-21-14',
        #     '2024_08_27-12-46-41',
        #     '2024_08_27-14-30-18',
        #     '2024_08_27-15-19-30',
        #     '2024_08_27-15-44-31',
        #     '2024_08_27-16-11-31',
        #     '2024_08_28-13-59-30',
        #     '2024_08_28-14-26-46',
        #     '2024_08_28-15-01-17',
        #     '2024_08_28-15-26-20',
        #     '2024_08_28-15-51-23',
        #     '2024_08_28-16-16-08',
        #     '2024_08_29-14-06-46',
        #     '2024_08_29-14-32-16',
        #     '2024_08_29-14-57-04',
        #     '2024_08_29-15-24-06',
        #     '2024_09_11-15-33-02',
        #     '2024_09_11-16-06-39',
        #     '2024_09_11-16-33-38',
        #     '2024_09_11-17-04-46',
        #     '2024_12_26-13-53-20',
        #     '2024_12_26-14-19-20',
        #     '2024_12_26-14-44-22',
        #     '2024_12_26-15-19-14',
        #     '2024_12_26-15-44-29',
        #     '2025_01_06-14-03-07',
        #     '2025_01_06-14-28-07',
        #     '2025_01_06-15-03-22',
        #     '2025_01_06-15-31-32',
        #     '2025_01_06-15-58-15',
        #     '2025_01_06-16-23-43',
        #     '2025_01_17-10-36-55',
        #     '2025_01_17-11-06-52',
        #     '2025_01_17-11-32-35',
        #     '2025_01_23-13-19-27',
        #     '2025_01_23-14-07-05',
        #     '2025_02_03-14-40-35',
        #     '2025_02_03-15-01-59',
        #     '2025_02_03-15-23-15',
        #     ]
        # DATASET = scan_folder_for_NWBfiles(datasets[c]['datafolder'], from_filenames=filenames)
        
        # FILTER
        # 1) protocol type: contrast sensitivity
        cond = np.array([np.sum(['ff-gratings-2orientations-8contrasts-15repeats' in p for p in protocols])\
                        for protocols in DATASET['protocols']], dtype=bool)
        # 2) age condition
        if datasets[c]['age_interval'] is not None:
            cond = cond &\
                (DATASET['ages']>=datasets[c]['age_interval'][0]) &\
                (DATASET['ages']<=datasets[c]['age_interval'][1])

        if len(DATASET['files'][cond])>nMIN_DATAFILES:

            if parallelized:
                ################################################
                ###    parallelization here !   #################
                ################################################
                nruns = int(len(DATASET['files'][cond])/cpus)+1

                for r in range(nruns):
                    i0 = r*cpus
                    imax = np.min([i0+cpus, len(DATASET['files'][cond])]) 
                    print(' - running set of files %i:%i' % (i0, imax))

                    # start the processes
                    procs = []
                    for i in range(i0,imax):
                        proc = multiprocessing.Process(\
                                            target=process_file, 
                                            args=(DATASET['files'][cond][i], 
                                                  i, c, quantity, filtering_cond_name))
                        procs.append(proc)
                        proc.start()

                    # complete the processes
                    for proc in procs:
                        proc.join()
            else:
                #####################################
                ###### UN-PARALLELIZED VERSION ######
                for i, f in enumerate(DATASET['files'][cond]):
                    process_file(f, i, c, quantity, filtering_cond_name)
                #####################################

            # now that we have stored all datafile outputs
            Sensitivities = []
            for i, f in enumerate(DATASET['files'][cond]):

                if os.path.isfile(os.path.join(summary_folder, 'temp', 
                                              'Sensitivity-%s-%s-%i.npy' % (quantity, c, i))):
                    Sensitivity = np.load(os.path.join(summary_folder, 'temp', 
                                                'Sensitivity-%s-%s-%i.npy' % (quantity, c, i)),
                                        allow_pickle=True).item()
                    Sensitivities.append(Sensitivity)

            # # saving data
            summary_prefix_name = get_summary_prefix_name(quantity, filtering_cond_name, special_dFoF_params)
            print(summary_prefix_name)

            np.save(os.path.join(summary_folder, summary_prefix_name + 'Sensitivities_%s.npy' % c), 
                    Sensitivities)

        else:
            print()
            print('   [!!]   DATASET NOT LARGE ENOUGH   [!!] ')
            print('               only N=%i sessions available' %\
                                        len(DATASET['files'][cond]))
            print('   [!!]   DATASET not analyzed       [!!] ')
            print()

        print('-----------------------------------------------------------------')
        print('=================================================================')
    #shutil.rmtree(os.path.join(summary_folder, 'temp'))

# %%
#f = '/home/user/CURATED/Cibele/SST-cells_cond-GluN1-KO_Adult_V1/NWBs/2026_02_18-14-42-43.nwb'
#d = Data(f, verbose=True)

# %%
if False:
    from Dataset_Organization_cibele import summary_folder
    from physion.analysis.protocols.contrast_sensitivity\
            import plot_contrast_sensitivity, plot_contrast_responsiveness
    fig, ax = plot_contrast_sensitivity(\
                            ['PV-cells_WT_Adult_V1_contrast-1.0', 
                             'PV-cells_WT_Adult_V1_contrast-0.5'],
                            #   average_by='ROIs',
                            path=summary_folder)


# %%
