import os
import numpy as np
import pymulti as pm


program_name='Multi-1D'
new_dir=os.path.join(program_name,'database')
n_splits=4
inp_data = pm.data1D_process_inp(program_name, n_splits ** 8)
V=np.zeros((n_splits ** 8,1))
processor = pm.DataProcessor()
for i in range(n_splits ** 8):
    processor.read_fort10_file(os.path.join(new_dir,'fort_'+str(i)+'.10'))
    V[i,0]= processor.extract('V')/1e5
np.save('1D_database_final.npy',np.concatenate((inp_data, V), axis=1))

