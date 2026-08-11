#%%

from scipy import stats
import matplotlib.pyplot as plt
from physion.utils import plot_tools as pt
from physion.analysis.protocols.orientation_tuning import *

#%%
def plot_tuning_responses_of_sign_and_non_sign_cells(folders, ylims, averaged_by_sessions = True) :
    
    """
    comments : not normalised of course
    """

    for i, folder in enumerate(folders) : 

        keys =  ['%s_contrast-1.0' % folder, 
                    '%s_contrast-0.5' % folder]

        RESPONSES_sign = []
        RESPONSES_non_sign = []

        for key in keys : 

            Tunings = np.load('/home/user/DATA/Astrid/Cibele_data/summary/mccf_Tunings_%s.npy' % key, allow_pickle=True)   
            
            Responses_sign = []
            Responses_non_sign = []

            for Tuning in Tunings : 
                Responses_sign.append([Tuning['Responses'][i] for i in range(len(Tuning['significant_ROIs'])) if Tuning['significant_ROIs'][i]])
                Responses_non_sign.append([Tuning['Responses'][i] for i in range(len(Tuning['significant_ROIs'])) if ~Tuning['significant_ROIs'][i]])
            
            Responses_sign = [resp for resp in Responses_sign if len(resp)!=0] #get rid of empty lists ie when no sign rois
            Responses_non_sign = [resp for resp in Responses_non_sign if len(resp)!=0] #get rid of empty lists ie when no non sign rois

            RESPONSES_sign.append(Responses_sign)
            RESPONSES_non_sign.append(Responses_non_sign)

        if averaged_by_sessions : 
            y_sign0 = np.mean([np.mean(RESPONSES_sign[0][i], axis = 0) for i in range(len(RESPONSES_sign[0]))], axis = 0)
            y_sign1 = np.mean([np.mean(RESPONSES_sign[1][i], axis = 0) for i in range(len(RESPONSES_sign[1]))], axis = 0)
            y_non_sign0 = np.mean([np.mean(RESPONSES_non_sign[0][i], axis = 0) for i in range(len(RESPONSES_non_sign[0]))], axis = 0)
            y_non_sign1 = np.mean([np.mean(RESPONSES_non_sign[1][i], axis = 0) for i in range(len(RESPONSES_non_sign[1]))], axis = 0)

        elif not averaged_by_sessions : 
            y_sign0 = np.mean(np.concatenate(RESPONSES_sign[0]), axis = 0)
            y_sign1 = np.mean(np.concatenate(RESPONSES_sign[1]), axis = 0)
            y_non_sign0 = np.mean(np.concatenate(RESPONSES_non_sign[0]), axis = 0)
            y_non_sign1 = np.mean(np.concatenate(RESPONSES_non_sign[1]), axis = 0)


        #plot

        colors=[pt.tab10(i), 'lightgrey']
        fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(9, 3), sharey = True)
        fig.tight_layout(w_pad = 5.5)

        x = np.linspace(-30, 180-30, 100)
        for j,y in enumerate([y_sign0, y_sign1, y_non_sign0, y_non_sign1]) : 

            C, func = fit_gaussian(Tunings[0]['shifted_angle'], y)
            axes[int(j/2)].plot(x, func(x), lw=2, alpha=.5,color = colors[np.mod(j,2)], label = keys[0][-4:])
            axes[int(j/2)].scatter(Tuning['shifted_angle'], y, color = colors[np.mod(j,2)], label = keys[1][-4:])


            axes[int(j/2)].set_xticks(ticks = Tunings[0]['shifted_angle'], labels = ['%i' % a if (a in [0, 90]) else '' for a in Tunings[0]['shifted_angle'] ])
            axes[int(j/2)].set_yticks(np.arange(3)*0.5)
            axes[int(j/2)].set_xlabel('angle ($^o$) from pref.')
            axes[int(j/2)].set_ylabel('$\\delta$ $\\Delta$F/F')



        axes[0].set_title('Responses of responsive cells')
        axes[1].set_title('Responses of non responsive cells')
        
        axes[0].text(x = -50, y= -0.5, s = folder)
        axes[0].text(x = -50, y= -0.6, s = 'averaged by sessions = ' + str(averaged_by_sessions))
        axes[0].text(x = -50, y= -0.7, s = 'multiple comparaison correction = False')
        axes[0].set_ylim(ylims[i])




#effect of ROI_TO_NEUROPIL_INCLUSION_FACTOR on nb of valid rois 
def nb_valid_rois_by_inclusion_factor(folders, summary_folder, neuropil_inclusion_factors, neuropil_correction_factor = 0.7) :
    
    plt.figure(figsize = (7,7))

    for i,folder in enumerate(folders):

        nb_rois_per_factor = []

        for neuropil_inclusion_factor in neuropil_inclusion_factors :
        
            Tunings = np.load(os.path.join(summary_folder, str('corrfact_') + str(neuropil_correction_factor) + 'inclufact_' + str(neuropil_inclusion_factor) + 'Tunings_%s.npy' % key), allow_pickle=True)
            nb_rois_per_factor.append(np.sum([Tuning['nROIs_final'] for Tuning in Tunings]))

        plt.plot(neuropil_inclusion_factors, nb_rois_per_factor, label = '%s_contrast-1.0' % folders[i], color = pt.tab10(i), alpha = 0.7)
        plt.scatter(neuropil_inclusion_factors, nb_rois_per_factor, color = pt.tab10(i))

    plt.text(s='fixed correction factor at : 0.7', x = -1, y = -120)
    plt.xlabel('ROI_TO_NEUROPIL_INCLUSION_FACTOR ')
    plt.ylabel('Total nb of valid rois')
    plt.legend()



#effect of NEUROPIL_CORRECTION_FACTOR on nb of valid and significant rois 
def nb_valid_and_sign_rois_by_correction_factor(folders, summary_folder, neuropil_correction_factors, neuropil_inclusion_factor = 1.15) :

    nb_rois_per_factor = [[]]*len(folders)
    nb_responsive_rois_per_factor = [[],[]]*len(folders)
    plt.figure(figsize = (7,7))

    for i, folder in enumerate(folders):
        keys = ['%s_contrast-1.0' % folder, 
                    '%s_contrast-0.5' % folder]

        for neuropil_correction_factor in neuropil_correction_factors :
        
            Tunings = np.load(os.path.join(summary_folder, str('corrfact_') + str(neuropil_correction_factor) + 'inclufact_' + str(neuropil_inclusion_factor) + 'Tunings_%s.npy' % keys[0]), allow_pickle=True)
            nb_rois_per_factor.append(np.sum([Tuning['nROIs_final'] for Tuning in Tunings]))
            nb_responsive_rois_per_factor[0].append(np.sum([np.sum(Tuning['significant_ROIs']) for Tuning in Tunings]))

            Tunings = np.load(os.path.join(summary_folder, str('corrfact_') + str(neuropil_correction_factor) + 'inclufact_' + str(neuropil_inclusion_factor) + 'Tunings_%s.npy' % keys[1]), allow_pickle=True)
            nb_responsive_rois_per_factor[1].append(np.sum([np.sum(Tuning['significant_ROIs']) for Tuning in Tunings]))



        plt.plot(neuropil_correction_factors, nb_rois_per_factor[0], label = 'valid rois', color = pt.tab10(i), alpha = 0.7)
        plt.scatter(neuropil_correction_factors, nb_rois_per_factor[0], color = pt.tab10(i))

        plt.plot(neuropil_correction_factors, nb_responsive_rois_per_factor[0], label = 'responsive rois' , color = pt.tab10(i), alpha = 0.7, lw = 2)
        plt.plot(neuropil_correction_factors, nb_responsive_rois_per_factor[1], label = 'responsive rois' , color = 'lightgrey', alpha = 0.7, lw = 2)
        plt.scatter(neuropil_correction_factors, nb_responsive_rois_per_factor[0], color = pt.tab10(i))
        plt.scatter(neuropil_correction_factors, nb_responsive_rois_per_factor[1], color = 'lightgrey')
        plt.text(x = -0.1, y = -275 -25*i, s='%s_contrast-1.0' % folders[i], color = pt.tab10(i))



    plt.text(x = -0.1, y = -300, s='fixed inclusion factor at : 1.15')

    plt.xticks(neuropil_correction_factors)
    plt.xlabel('NEUROPIL_CORRECTION_FACTOR ')
    plt.ylabel('Total nb of rois')
    #plt.ylim((500,700))
    plt.legend(loc='center', bbox_to_anchor=(0.8, -0.3))



def plot_orientation_tuning_variations_both_parameters(keys,
                      path=os.path.expanduser('~'),
                      average_by='sessions',
                      colors=None,
                      fig_args={'right':20, 'figsize' : (7,7)},
                      neuropil_inclusion_factors = [],
                      neuropil_correction_factors = []):
    

    fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(9, 3))
    fig.tight_layout(w_pad = 5.5)
    fig.text(x = -0.02, y= 1.12, s = keys[0][:-13])

    if type(keys)==str:
        keys, colors = [keys], [colors[0]]

    x = np.linspace(-30, 180-30, 100)

    for j,neuropil_inclusion_factor in enumerate(neuropil_inclusion_factors) :
        a = 0
        for i, (key, color) in enumerate(zip(keys, colors)):
            mean_nroi = []
            for k,neuropil_correction_factor in enumerate(neuropil_correction_factors) :
                # load data
                Tunings = \
                    np.load(os.path.join(path, 'corrfact_' + str(neuropil_correction_factor) + 'inclufact_' + str(neuropil_inclusion_factor) + 'Tunings_%s.npy' % key), 
                            allow_pickle=True)

                Responses = get_tuning_responses(Tunings,
                                            average_by=average_by)

                # Gaussian Fit
                C, func = fit_gaussian(Tunings[0]['shifted_angle'],
                                    np.nanmean([r/r[1] for r in Responses], axis=0))

                axes[j].plot(x, func(x), lw=2, color=color, alpha =0.3*(k+1))

                axes[j].scatter(Tunings[0]['shifted_angle'], np.nanmean([r/r[1] for r in Responses], axis=0),  color=color, alpha = 0.3*(k+1))
                
                axes[j].text(x = -60, y= -0.6-(k/10)+a, s = 'neuropil correction factor :' + str(neuropil_correction_factor), color=color, alpha = 0.3*(k+1))
                mean_nroi.append(np.sum([Tunings[i]['nROIs_responsive'] for i in range(len(Tunings))]))
                axes[j].text(x = 120 + a*700, y= -0.6-(k/10)+a, s = str(mean_nroi[-1]), color='black', alpha = 0.3*(k+1))

            axes[j].text(x = 120+ a*700, y= -0.95, s = str(np.round(np.mean(mean_nroi),1)), color='black')
            a = 0.05

        axes[j].set_xticks(ticks = Tunings[0]['shifted_angle'], labels = ['%i' % a if (a in [0, 90]) else '' for a in Tunings[0]['shifted_angle'] ])
        axes[j].set_yticks(np.arange(3)*0.5)
        axes[j].set_xlabel('angle ($^o$) from pref.')
        axes[j].set_ylabel('norm. $\\delta$ $\\Delta$F/F')
        axes[j].set_ylim([-0.05, 1.05])
        axes[j].text(x = -60, y= -0.4, s = 'neuropil inclusion factor = ' + str(neuropil_inclusion_factor))
        axes[j].text(x = 120, y= -0.4, s = '/ nROIs resp')
    



def plot_orientation_tuning_curve_multiple_corrfact_one_graph(keys,
                      path=os.path.expanduser('~'),
                      average_by='sessions',
                      colors=None,
                      neuropil_correction_factors = [],
                      neuropil_inclusion_factor = 1.15):
        
    if colors is None:
        colors = pt.plt.rcParams['axes.prop_cycle'].by_key()['color']

    if type(keys)==str:
        keys, colors = [keys], [colors[0]]

    plt.figure(figsize = (6,6))
    x = np.linspace(-30, 180-30, 100)

    a=0 
    for i, (key, color) in enumerate(zip(keys, colors)):
        for k,neuropil_correction_factor in enumerate(neuropil_correction_factors) :
            # load data
            Tunings = \
                np.load(os.path.join(path, 'corrfact_' + str(neuropil_correction_factor) + 'inclufact_' + str(neuropil_inclusion_factor) + 'Tunings_%s.npy' % key), 
                        allow_pickle=True)

            Responses = get_tuning_responses(Tunings,
                                           average_by=average_by)

            # Gaussian Fit
            C, func = fit_gaussian(Tunings[0]['shifted_angle'],
                                np.mean([r/r[1] for r in Responses], axis=0))

            plt.scatter(Tunings[0]['shifted_angle'], np.mean([r/r[1] for r in Responses], axis=0), 
                        color=color, alpha = 0.3+k/10)


            plt.plot(x, func(x), lw=2, color=color, alpha =0.3+ k/10)

            plt.text(x = -60, y= -0.5-(k/10)+a, s = 'neuropil correction factor :' + str(neuropil_correction_factor), color=color, alpha =0.3+ k/10)


        a=0.05

    plt.xticks(ticks = Tunings[0]['shifted_angle'], labels = ['%i' % a if (a in [0, 90]) else '' for a in Tunings[0]['shifted_angle'] ])
    plt.yticks(np.arange(3)*0.5)
    plt.xlabel('angle ($^o$) from pref.')
    plt.ylabel('norm. $\\delta$ $\\Delta$F/F')
    plt.ylim([-0.05, 1.05])
    plt.title('NEUROPIL_CORRECTION_FACTOR variations')
    plt.text(x = -60, y= -0.4, s = key[:-13])
    plt.text(x = -60, y= -0.3, s = 'ROI_TO_NEUROPIL_INCLUSION_FACTOR = 1.15')
