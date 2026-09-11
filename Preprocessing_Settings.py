    
def get_dFoF_params(dataset, special_dFoF_params):

    dFoF_parameters = dict(\
                roi_to_neuropil_fluo_inclusion_factor=1.15, # no factor here
                neuropil_correction_factor = 0.7,
                method_for_F0 = 'sliding_percentile',
                percentile=5., # percent
                sliding_window = 5*60, # seconds
        )
    
    if dataset[:4]=='PYR-':
        # means pyramidal cells
        dFoF_parameters['roi_to_neuropil_fluo_inclusion_factor'] = 0 

    if special_dFoF_params is not None :
        dFoF_parameters['roi_to_neuropil_fluo_inclusion_factor'] = special_dFoF_params['neuropil_inclusion_factor']
        dFoF_parameters['neuropil_correction_factor'] = special_dFoF_params['neuropil_correction_factor']

    return dFoF_parameters

def get_stat_test_props(dataset): 

    stat_test_props=dict(interval_pre=[-1.,0],
                        interval_post= [1.,2.],                                   
                        test='ttest',                                            
                        sign='positive')

    if  dataset[:4]=='PYR-' or  dataset[:3]=='PV-' : 
        stat_test_props['interval_post'] = [0,1.]

    return stat_test_props

