import time
start_time=time.time()
import sys
import os
import moses
import pandas as pd
import numpy as np
print(f'loaded libraries in : {time.time()-start_time:.2f}s')

### changing the working directory 
script_dir = os.path.dirname(os.path.realpath(__file__))
if __name__ == '__main__':
    sys.path.append(os.path.join(script_dir, '../'))
print(f'working dir is {os.getcwd()}')


def download_moses():
    ### load and concatenate SMILES from MOSES datasets 
    print('>>> Loading data from moses')
    smiles_list = moses.get_dataset('train')
    smiles_list= np.hstack([smiles_list,moses.get_dataset('test')]) 
    smiles_list= np.hstack([smiles_list,moses.get_dataset('test_scaffolds')]) 
    
    ### format to Panda dataframe
    dataset = pd.DataFrame(smiles_list).rename(columns={0: 'smiles'})
    
    ### save as csv file 
    print('>>> Saving data to csv file in data dir')
    dataset.to_csv('data/whole_MOSES_smiles.csv',index=False)

download_moses()
print(f'it took {time.time()-start_time:.2f}s to execute the whole script')