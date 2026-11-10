#%%
from physion.analysis.read_NWB import Data
import numpy as np
import matplotlib.pyplot as plt



def plot_file_photodiode_shift(filename):

    #build
    data = Data(filename, verbose=False)
    data.build_rawFluo()
    quantities = ['rawFluo']

    if data.has_photodiode(): 
        data.build_photodiode()
        quantities.append('photodiode')

    Episodes = EpisodeData(data, 
                            quantities=quantities, 
                            protocol_name = data.protocols[0],
                            verbose=False)

    Episodes.rawFluo = np.mean(Episodes.rawFluo, axis = 1)
    N_ep = Episodes.rawFluo.shape[0]

    #compute values of interest
    period_evoked = []
    for quantity in quantities : 
        val = getattr(Episodes, quantity)

        for i in range(3):
            period_start = int((i/3) * N_ep)
            period_stop = int(((i+1)/3)*N_ep)
            period_evoked.append(np.mean(val[period_start:period_stop, :], axis = 0)) 

    #plot

    mosaic = "ABC" 
    fig = plt.figure(layout="constrained", figsize = (12,4))

    ax_dict = fig.subplot_mosaic(mosaic)
    labels = ['1st period', '2nd period', '3rd period']
    for i in range(3):
        ax_dict["A"].plot(Episodes.t, period_evoked[i]/np.max(period_evoked[i]), label = labels[i])
    ax_dict["A"].set_xlabel('time')
    ax_dict["A"].set_ylabel('norm mean rawFluo')
    ax_dict["A"].legend()

    if data.has_photodiode(): 
        for i in range(3,6):
            ax_dict["B"].plot(Episodes.t, period_evoked[i])

        ax_dict["B"].set_xlabel('time')
        ax_dict["B"].set_ylabel('mean Photodiode')

        mask_first_min = [data.t_photodiode <= 60][0]
        ax_dict["C"].plot(data.t_photodiode[mask_first_min], data.photodiode[mask_first_min], color = 'black')
        ax_dict["C"].set_xlabel('time')
        ax_dict["C"].set_ylabel('Photodiode')
    else : 
        ax_dict["B"].annotate(text = 'no photodiode', xy = (0.5,0.5), fontsize = 15)
        ax_dict["B"].axis("off")
        ax_dict["C"].axis("off")
    
    return fig 

