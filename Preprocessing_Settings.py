    
def get_dFoF_params(dataset, special_dFoF_params):

    if dataset[:4]=='PYR-':
        # means pyramidal cells
        dFoF_parameters = dict(\
                roi_to_neuropil_fluo_inclusion_factor=0., # no factor here
                neuropil_correction_factor = 0.7,
                method_for_F0 = 'sliding_percentile',
                percentile=5., # percent
                sliding_window = 5*60, # seconds
        )
    else:
        # means interneurons
        dFoF_parameters = dict(\
                roi_to_neuropil_fluo_inclusion_factor=1.15,
                neuropil_correction_factor = 0.7,
                method_for_F0 = 'sliding_percentile',
                percentile=5., # percent
                sliding_window = 5*60, # seconds
                with_correctedFluo_and_F0=True,
        )

    if special_dFoF_params is not None :
        dFoF_parameters['roi_to_neuropil_fluo_inclusion_factor'] = special_dFoF_params['neuropil_inclusion_factor']
        dFoF_parameters['neuropil_correction_factor'] = special_dFoF_params['neuropil_correction_factor']

    return dFoF_parameters
