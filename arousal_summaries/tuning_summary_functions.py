#%%
import os, sys , shutil 
import numpy as np
from scipy import stats
import pandas as pd 

os.chdir('/home/user/lab-notebook/astrid')
sys.path += ['./physion/src']
import physion.utils.plot_tools as pt
from physion.analysis.read_NWB\
                         import scan_folder_for_NWBfiles, Data
from physion.analysis.episodes.build import EpisodeData
from physion.analysis.protocols.orientation_tuning import *
from arousal_summaries.arousal_common_fcts import * 

#%%


def get_prefered_angles(prefered_angles_dataset, filename):

    short_filename = os.sep.join(filename.split('/')[-3:])
    prefered_angles = prefered_angles_dataset['prefered_angles'][(prefered_angles_dataset['short_filenames'] == short_filename)]

    if not prefered_angles.empty : 
        return prefered_angles.values[0]
    else : 
        raise ValueError('file %s not present in the prefered_angles_dataset' % short_filename)


def get_prefered_angles_dataset(data_folder, summary_folder, contrast) :

    folder = data_folder.split('/')[-2]
    prefered_angles_dataset = pd.DataFrame({})
    try : 
        Tunings = np.load(summary_folder + '/Tunings_%s_contrast-%s.npy' % (folder, contrast), allow_pickle=True)
        prefered_angles_dataset['short_filenames'] = [os.sep.join(Tuning['datafile'].split('/')[-3:]) for Tuning in Tunings] #ignore basepath
        prefered_angles_dataset['contrasts'] = [contrast]*len(prefered_angles_dataset['short_filenames'])
        prefered_angles_dataset['prefered_angles'] = [Tuning['prefered_angles'] for Tuning in Tunings]
        
    except: 
        print('No Tunings file found for folder %s and contrast %s' % (folder, contrast))
        print('You need to build the classical tunings summary first, before building one with a filtering condition')


    return prefered_angles_dataset


#%%



#--------------------------#
#a priori this one can be delete. to do after checking creation new run tuning summary are ok
def compute_tuning_response_per_cells_with_arousal_cond(data, Episodes, 
                                      arousal_cond,
                                      stat_test_props,
                                      response_significance_threshold = 0.05,
                                      filtering_cond=None,
                                      quantity='dFoF',
                                      contrast=1.0,
                                      nMin_episodes = 2,
                                      start_angle=-22.5, 
                                      angle_range=180,
                                      verbose=False) : 
    """
    /!\ Choice of significance based on summary without condition. 
    If statistical signficiance test were to be done on the actual run response => way less significant ROI ie cf commented lines 
    """
    Tuning_wo_cond = compute_tuning_response_per_cells(data, Episodes, 
                                                quantity='dFoF', 
                                                stat_test_props = stat_test_props, 
                                                response_significance_threshold = response_significance_threshold, 
                                                filtering_cond = None,
                                                contrast = contrast,
                                                nMin_episodes = nMin_episodes,
                                                start_angle=-22.5, 
                                                angle_range=180,
                                                verbose=False)


    if arousal_cond == '' : 
        Tuning = Tuning_wo_cond

    else : 
        Tuning  = {}
        for key in  Tuning_wo_cond.keys() : 
            if key == 'ntrials' : 
                temp_arr = np.zeros(Tuning_wo_cond['Responses'].shape)
            else : 
                temp_arr = np.zeros(Tuning_wo_cond[key].shape)
            temp_arr[:] = np.nan
            Tuning[key] = temp_arr
        
        Tuning['prefered_angles'] = Tuning_wo_cond['prefered_angles']
        Tuning['shifted_angle'] = Tuning_wo_cond['shifted_angle']
        Tuning['significant_ROIs'] = Tuning_wo_cond['significant_ROIs']

        filtering_arousal_cond = get_arousal_filtering_cond(arousal_cond, Episodes)

        cond = Episodes.find_episode_cond(key='contrast', 
                                            value=contrast) &\
                                            filtering_arousal_cond
        
        summary = Episodes.pre_post_statistics(stat_test_props = stat_test_props,
                                                    episode_cond=cond,
                                                    repetition_keys=['repeat', 'contrast'],
                                                    nMin_episodes = nMin_episodes,
                                                    response_args=dict(quantity=quantity),
                                                    response_significance_threshold=response_significance_threshold,
                                                    multiple_comparison_correction=False,
                                                    loop_over_cells=True,
                                                    verbose=verbose)
        #significant = np.zeros(data.nROIs)
        for i,angle in enumerate(summary['angle']) : 
            for roi in range(data.nROIs) : 
                new_angle = shift_orientation_according_to_pref(angle,
                                                            pref_angle=Tuning['prefered_angles'][roi],
                                                            start_angle=start_angle,
                                                            angle_range=angle_range)
                iangle = np.flatnonzero(Tuning['shifted_angle']==new_angle)[0]

                Tuning['Responses'][roi][iangle] = summary['value'][roi,i]
                Tuning['std-values'][roi][iangle] = summary['std-value'][roi,i]
                Tuning['ntrials'][roi][iangle] = summary['ntrials'][i]
                #significant[roi] += np.sum(summary['significant'][roi])

            #new_angle = shift_orientation_according_to_pref(summary['angle'][i],
            #                                        pref_angle=0,
            #                                        start_angle=start_angle,
            #                                        angle_range=angle_range)
            #iangle = np.flatnonzero(Tuning['shifted_angle']==new_angle)[0]

            #Tuning['ntrials'][iangle] = summary['ntrials'][i]
        
        #Tuning['significant_ROIs'] = significant > 0 
    return Tuning 

