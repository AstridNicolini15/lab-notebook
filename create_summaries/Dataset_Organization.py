import os
import sys 

folders = [
    "PV-cells_WT_Adult_V1", 
    #"PV-cells_WT_Young_V1", 
    #"PV-cells_cond-GluN1-KO_Adult_V1", 
    #"PYR-SynGCaMP_WT_V1",
    #"PYR-cells_WT_Adult_V1",
    #"SST-cells_cond-GluN1-KO_Young_V1",
    #"SST-cells_cond-GluN1-KO_Adult_V1",
    #"SST-cells_WT_Adult_V1",
    #"SST-cells_WT_Young_V1",
    #"Wild-Type",
    #"GluN1-KO",
    #"SST_WT",
    #"SST_GluN1",
    #"SST_GluN3"
]
#base_path =  "/home/user/DATA/Astrid/Taddy_data/OneDrive_1_9-3-2026" #"/home/user/DATA/Astrid/Taddy_data/photodiode_aligned" # "/home/user/DATA/Astrid/Taddy_data/OneDrive_1_9-3-2026" #"/home/user/DATA/Astrid/Taddy_data"
base_path =  '/home/user/DATA/Astrid/Cibele_data'
summary_folder = base_path + '/summary'

quantity = 'dFoF' #'Deconvolved_correctedFluo' or 'Deconvolved', Zscore_correctedFluo'
filtering_cond_name = None
keep_all_orientations_significances = False
#sys.argv = {#'roi_to_neuropil_fluo_inclusion_factor' : neuropil_inclusion_factor,
            #'neuropil_correction_factor' : neuropil_correction_factor,
            #'with_computed_neuropil_fact' : True,
            #'computed_neuropil_fact_method' : 'slope'}

# age intervals in Young
AGE_INTERVALS = [\
    (15,19), (20,23), (24,27) #, (16,21), (22,27)
    ]

# to be a valid dataset:
nMIN_DATAFILES = 2

parallelized = False
debug = False

def datasets_func(key='contrast',
        values=[0.5,1.0]):

    datasets = {}
    for c in folders:

        for val in values:

            # we split young animals into age groups
            if 'Young' in c:
                for interval in AGE_INTERVALS:
                    datasets[c.replace('Young', 'P%i-P%i' % interval)+'_%s-%.1f' % (key,val)] =\
                        {'datafolder':os.path.join(base_path, c, 'NWBs'), 
                            'age_interval':interval}
            else : 
                datasets[c+'_%s-%.1f' % (key, val)] =\
                    {'datafolder':os.path.join(base_path, c, 'NWBs'), 
                        'age_interval':None}
                
    return datasets


if __name__=='__main__':
    from pprint import pprint
    pprint(datasets)

