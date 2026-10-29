#%%
import os, sys 
os.chdir('/home/user/lab-notebook/astrid')
import numpy as np
import itertools
from create_summaries.Preprocessing_Settings import get_dFoF_params

def extract_special_dict(special_dict): 

    corrected_special_dict = { 'title' : '', 'ylims' : None}

    if special_dict is not None :
        for key in special_dict.keys() :
            corrected_special_dict[key] = special_dict[key]

    title = corrected_special_dict['title']
    ylims = corrected_special_dict['ylims']
    
    return title, ylims


def get_protocol_keys(summary_protocol = 'Tunings') : 

    if summary_protocol == 'Tunings' : 
        keys = ['_contrast-1.0', 
            '_contrast-0.5' ]

    elif summary_protocol == 'Sensitivities' :
        keys =  ['_angle-90.0' , 
            '_angle-0.0' ] 

    return keys

         
def get_nwb_filenames_from_summary(folder, summary_path, summary_protocol):

    key = get_protocol_keys(summary_protocol = summary_protocol)[0]

    summary= np.load(summary_path + '/' + summary_protocol + '_' + folder + key + '.npy', allow_pickle=True)  #only one control cond is sufficient
    nwb_filenames = [s['datafile'] for s in summary]

    return nwb_filenames


def get_summary_names(folders, summary_paths, summary_protocols, summary_prefixe_names = [''], summary_sufixe_names = [''], template = [], add_keys = False): 

    def input_form_check(input, ref = [1]): 
        if type(input) == str:
            input = [input]*len(ref)
        if type(input) == list and len(input) != len(ref):
            input = input *len(ref)
        return input
    
    folders = input_form_check(folders)
    summary_protocols = input_form_check(summary_protocols, folders)
    summary_paths = input_form_check(summary_paths, folders)

    summary_filenames = []

    for i,folder in enumerate(folders) : 
        summary_prefixe_names = ['']
        summary_sufixe_names = ['']

        for item in template : 

            if item == 'multiple_corfact': 
                for neuropil_correction_factor in [0,0.15,0.3,0.55,0.7,0.85,1] : 
                    inclufact = get_dFoF_params(folder, None)['roi_to_neuropil_fluo_inclusion_factor']
                    summary_prefixe_names.append('roi_to_neuropil_fluo_inclusion_factor_' + str(inclufact) + '_neuropil_correction_factor_' + str(neuropil_correction_factor) + '_') 
            
            if item == 'inclufact_2_3' : 
                #default_inclufact = get_dFoF_params(folder, None)['roi_to_neuropil_fluo_inclusion_factor']
                for inclufact in [2, 3]:
                    summary_prefixe_names.append('roi_to_neuropil_fluo_inclusion_factor_' + str(inclufact) + '_neuropil_correction_factor_0.7_') 

            if item == '' and type(template) == list: 
                summary_prefixe_names.append('')

            if item == 'young' : 
                for interval in (15,19), (20,23), (24,27) : 
                        summary_sufixe_names.append('_P%i-P%i_' % interval)

        if summary_prefixe_names[0] == '' and len(summary_prefixe_names) > 1 : 
            del summary_prefixe_names[0]
        if summary_sufixe_names[0] == '' and len(summary_sufixe_names) > 1 : 
            del summary_sufixe_names[0]

        for summary_prefixe_name, summary_sufixe_name in list(itertools.product(summary_prefixe_names,summary_sufixe_names)) :

            if add_keys : 
                keys = [key + '.npy' for key in get_protocol_keys(summary_protocols[i])]
            else : 
                keys = ['']
                
            for key in keys : 
                    summary_name = summary_paths[i] + '/' + summary_prefixe_name + summary_protocols[i] + '_' + folder + summary_sufixe_name + key 
                    summary_filenames.append(summary_name)

    return summary_filenames
