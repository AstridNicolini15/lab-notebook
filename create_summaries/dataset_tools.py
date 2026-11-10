from notebooks.arousal import compute_arousal_masks, check_presence_locomotion_values

def get_summary_prefix_name(quantity, filtering_cond_name, special_dFoF_params) :

    if quantity == 'dFoF' : 
        summary_prefix_name = ''
    else : 
        summary_prefix_name = quantity + '_'

    if filtering_cond_name is not None : 
        summary_prefix_name += filtering_cond_name  + '_' #Run or Rest for exemple

    if special_dFoF_params is not None: 
        for key in special_dFoF_params.keys(): 
            summary_prefix_name += key + '_' + str(special_dFoF_params[key]) + '_'

    return summary_prefix_name


def build_filtering_cond_quantities(filtering_cond_name, data, quantities):

    if filtering_cond_name == 'Run' or filtering_cond_name == 'Rest' :
        if data.has_running():
            quantities += ['running']
            data.build_running(verbose=False)
            check_presence_locomotion_values(data)

    return data, quantities


def get_filtering_cond(filtering_cond_name, Episodes):

    if filtering_cond_name == None : 
        filtering_cond = None
    
    elif filtering_cond_name == 'Run' : 
        filtering_cond =  compute_arousal_masks(Episodes)[0]

    elif filtering_cond_name == 'Rest' : 
        filtering_cond =  compute_arousal_masks(Episodes)[1]

    return filtering_cond



def build_quantity(data, quantity, dFoF_parameters, stat_test_props):

    if '_' in quantity : #ie calculation on top of a quantity
        quantity_type = quantity.split('_')[0]
        quantity = quantity.split('_')[1]
    else : 
        quantity_type = None

    if quantity in ['correctedFluo', 'F0'] : 
        dFoF_parameters['with_correctedFluo_and_F0'] = True

    if quantity =='dFoF' :
        data.build_dFoF(**dFoF_parameters, verbose=False)
    elif quantity == 'correctedFluo' :
        data.build_correctedFluo(specific_time_sampling=None, interpolation='linear',
                            neuropil_correction_factor=dFoF_parameters['neuropil_correction_factor'],
                            verbose=False)
    #or build only the quantity of interests
    #getattr(data, 'build_%s' % quantity)(quantity, **dFoF_parameters)

    #build the calculation on top of the quantity if needed
    if quantity_type is not None :
        if quantity_type == 'Zscore' :
            data.build_Zscore(quantity)

        if quantity_type == 'Deconvolved':
            data.build_Deconvolved(Tau = 1.5, quantity = quantity)
            stat_test_props['interval_post']=[0.0, 1.0] 

    return data, quantity, dFoF_parameters, stat_test_props

