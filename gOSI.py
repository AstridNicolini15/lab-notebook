#%% gOSI 
import sys
sys.path += ["/home/user/lab-notebook/astrid/physion/src"]

from physion.analysis.protocols.orientation_tuning import get_tuning_responses
from math import radians




#%% Cibele
summary_path = "/home/user/DATA/Astrid/Cibele_data/summary"

folders = [
    "SST-cells_WT_Adult_V1",
    "SST-cells_cond-GluN1-KO_Adult_V1",
    ]

colors = [[(to_rgb('#12522eff'), 0.8),'lightgrey'], [(to_rgb('#5d1490ff'), 0.8),'lightgrey']]

#%% Taddy
summary_path = "/home/user/DATA/Astrid/Taddy_data/OneDrive_1_9-3-2026/summary"

folders = [
    "Wild-Type",
    "GluN1-KO"
    ]
colors = [[pt.tab10(1),'lightgrey'], [pt.tab10(4),'lightgrey']]

#%% plot inclusion factor curves and corresponding gOSI (average of sessions averages)

mosaic = """
    AABBCCDDD
    AABBCCDDD
    """

fig = plt.figure(layout="constrained", figsize = (20,5))
ax_dict = fig.subplot_mosaic(mosaic, gridspec_kw={"wspace": 0.3})


for neuropil_inclusion_factor, ax in zip([1.15, 2, 3], [ax_dict["A"], ax_dict["B"], ax_dict["C"]]) : 
    special_dict = {'name' : 'inclusion_factor-' + str(neuropil_inclusion_factor) + '_correction_factor-0.7_' , 
                    'title' : 'neuropil inclusion \n factor = '+ str(neuropil_inclusion_factor),
                    'ylims' : (0,1.1)}
    
    plot_tuning_responses_many_pop(folders, colors, ax, special_dict = special_dict, summary_path = summary_path)

#colors = [pt.tab10(1)]*3 + [(pt.tab10(1)[:3], 0.3)]*3 + [pt.tab10(4)]*3 + [(pt.tab10(4)[:3], 0.3)]*3
colors = [to_rgb('#12522eff')]*3 + [(to_rgb('#12522eff'), 0.3)]*3 + [to_rgb('#5d1490ff')]*3 + [(to_rgb('#5d1490ff'), 0.3)]*3



x = [0,1,2,4,5,6,10,11,12,14,15,16]
height = []
labels = []
plot_osi_final_curve = False
for folder in folders : 

    for key in ['%s_contrast-1.0' % folder,
                '%s_contrast-0.5' % folder]:

        for inclusion_factor in [1.15, 2, 3] : 

            labels.append(key)
            prefixe = 'inclusion_factor-' + str(inclusion_factor) + '_correction_factor-0.7_'
            Tunings = np.load(summary_path + '/' + prefixe + 'Tunings_%s.npy' % key, allow_pickle=True)

            if plot_osi_final_curve : 
                Responses = get_tuning_responses(Tunings)
                Response = np.nanmean(Responses, axis=0)
                gOSI = get_cell_gOSI(Response, Radians) 

                height.append(gOSI)

            else : 
                gOSIs = get_gOSIs(Tunings)
                height.append(np.mean([np.mean(session_gOSIs) for session_gOSIs in gOSIs]))
            

ax_dict["D"].grid(True, linewidth = 1, axis = 'y')
ax_dict["D"].bar(x = x, height = height, width = 0.75, color = colors)
ax_dict["D"].set_xticks(ticks = x, labels = [1.15,2,3,1.15,2,3,1.15,2,3,1.15,2,3], rotation = 45)
ax_dict["D"].set_xlabel('Neuropil inclusion factor')
ax_dict["D"].set_ylabel('gOSI')
ax_dict["D"].set_ylim(0,0.45)

#%%


#%%
summary_path = "/home/user/DATA/Astrid/Cibele_data/summary"

folders = [
    "SST-cells_WT_Adult_V1",
    "SST-cells_cond-GluN1-KO_Adult_V1",
    ]

colors = [to_rgb('#12522eff'), (to_rgb('#5d1490ff'))]
#%%

from matplotlib.transforms import Bbox
summary_path = "/home/user/DATA/Astrid/Taddy_data/OneDrive_1_9-3-2026/summary"

folders = [
    "Wild-Type",
    "GluN1-KO"
    ]
colors = [pt.tab10(1), pt.tab10(4)]


#%%plot cell to cell gOSI
for i,folder in enumerate(folders) : 

    gOSIs = []
    concatenate_significant_cells = []

    for key in ['%s_contrast-1.0' % folder,
                '%s_contrast-0.5' % folder]:

        Tunings = np.load(summary_path + '/' + 'Tunings_%s.npy' % key, allow_pickle=True)

        gOSIs.append(np.concatenate(get_gOSIs(Tunings, only_significant = False)))
        concatenate_significant_cells.append(np.concatenate([T['significant_ROIs'] for T in Tunings]))

    significant_cells_in_both_cond = concatenate_significant_cells[0]*concatenate_significant_cells[1]

    diff_gOSI = gOSIs[0][significant_cells_in_both_cond] - gOSIs[1][significant_cells_in_both_cond]
    mean_diff_gOSI = np.mean(diff_gOSI)

    #plot
    fig,ax = plt.subplots(1,2,figsize = (11,7))
    ax[0].set_position(Bbox([[0, 0], [0.7, 1]]))
    ax[1].set_position(Bbox([[0.8, 0], [1, 0.3]]))
    for y1,y2 in zip(gOSIs[0][significant_cells_in_both_cond], gOSIs[1][significant_cells_in_both_cond]) : 
        ax[0].plot([0,1], [y1,y2], color = 'black', linewidth = 0.2)
        ax[0].scatter([0,1], [y1,y2], color = colors[i])


    ax[0].set_xticks(ticks = [0,1], labels = ['contrast = 1', 'contrast = 0.5'])
    ax[0].set_xlim(-0.15,1.15)
    ax[0].set_ylim(-0.05,0.9)
    ax[0].set_ylabel('gOSI of cells significant in both contrast')

    counts, bins = np.histogram(diff_gOSI, 20)
    ax[1].stairs(counts/len(diff_gOSI), bins, color = colors[i])
    ax[1].set_xlim(-0.8,0.8)
    ax[1].set_yticks(ticks = [0,0.125,0.25], labels = [0,0.125,0.25])
    ax[1].set_ylabel('%')
    ax[1].set_xlabel('diff gOSI values')



    ax[0].annotate(text = 'mean diff gOSI \n =' + str(np.round(mean_diff_gOSI, 3)), xy = (1.05,0.06), fontsize = 9)
    ax[0].annotate(text = folder, xy = (1.05,0.80), color = colors[i])
    ax[0].annotate(text = 'N cell = ' + str(np.sum(significant_cells_in_both_cond)), xy = (1.05,0.13), color = colors[i])

        

# %%

def get_gOSIs(Tunings, only_significant = True) : 

    gOSIs = []

    angles = Tunings[0]['shifted_angle']
    Radians = [radians(2*(angle + 22.5)) for angle in angles]

    for T in Tunings : 
        session_gOSIs = []

        if only_significant : 
            Responses = T['Responses'][T['significant_ROIs']]
        else : 
            Responses = T['Responses']

        for R in Responses : #to note : do not work with summary with 2d significances

            cell_gOSI = get_cell_gOSI(R, Radians)
            session_gOSIs.append(cell_gOSI)

        gOSIs.append(session_gOSIs)
    
    return gOSIs


def get_cell_gOSI(Responses, Radians) : 

    """
    Responses is a 1D array of dim = (1*8)
    Radians is 1D array or list of dim = (1*8)
    """

    cell_gOSI = np.sqrt(np.sum(Responses*np.cos(Radians))**2 + np.sum(Responses*np.sin(Radians))**2) / np.sum(np.abs(Responses))

    return cell_gOSI