#%%
import os, sys
                          

neuropil_inclusion_factors = [2,3] #[1.15]
neuropil_correction_factors = [0.7] #[0,0.15,0.3,0.55,0.7,0.85,1]

for neuropil_inclusion_factor in neuropil_inclusion_factors : 

    for neuropil_correction_factor in neuropil_correction_factors : 
        print('COMPUTING : inclusion_factor = ' + str(neuropil_inclusion_factor) + ' and correction_factor = ' + str(neuropil_correction_factor))

        sys.argv = {'neuropil_inclusion_factor' : neuropil_inclusion_factor,
                    'neuropil_correction_factor' : neuropil_correction_factor}
        
        exec(open('/home/user/lab-notebook/astrid/Tuning-Dataset.py').read())
