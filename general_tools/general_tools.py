#%%
import os, sys 
os.chdir('/home/user/lab-notebook/astrid')
import numpy as np

def extract_special_dict(special_dict): #to generalize with keys (allow new key and to not have to put all the keys ex : no title)

    corrected_special_dict = {'name' : '', 'title' : '', 'ylims' : None}

    if special_dict is not None :
        for key in special_dict.keys() :
            corrected_special_dict[key] = special_dict[key]

    summary_prefix_name =  corrected_special_dict['name']
    title = corrected_special_dict['title']
    ylims = corrected_special_dict['ylims']
    
    return summary_prefix_name, title, ylims


def get_folder_keys(folder, summary_protocol = 'Tunings') : 

    if summary_protocol == 'Tunings' : 
        keys = ['%s_contrast-0.5' % folder, 
            '%s_contrast-1.0' % folder]

    elif summary_protocol == 'Sensitivities' :
        keys =  ['%s_angle-90.0' % folder, 
            '%s_angle-0.0' % folder] 

    return keys

         
def get_summary_protocol(summary):

    keys = summary[0].keys()
    if 'shifted_angle' in keys : 
        summary_protocol = 'Tunings'

    elif 'contrast' in keys :
        summary_protocol = 'Sensitivities'

    return summary_protocol

def correct_summary_form(summary_path, folders) :
    if type(summary_path) == str :
        summary_path = [summary_path] * len(folders)
        print('summary_path is not a list, assuming same path for all folders')
    return summary_path


def load_summaries_one_folder(folder, special_dict = None, summary_protocol = 'Tunings', summary_path =  '/home/user/DATA/Astrid/summary') : 

    summary_prefix_name, title, ylims = extract_special_dict(special_dict)
    keys = get_folder_keys(folder, summary_protocol = 'Tunings')
    summary_cond1 = np.load(summary_path + '/' + summary_prefix_name + summary_protocol + '_' + keys[0] + '.npy', allow_pickle=True)  #only one control cond is sufficient
    summary_cond2 = np.load(summary_path + '/' + summary_prefix_name +  summary_protocol + '_' + keys[1] + '.npy', allow_pickle=True)  #only one control cond is sufficient


    return summary_cond1, summary_cond2