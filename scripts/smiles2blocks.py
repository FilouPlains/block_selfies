
import time 
global_start = time.perf_counter()
import pandas as pd
import selfies as sf 
import os
import sys
from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np
import shutil
import multiprocessing as mp

from io import StringIO
from functools import partial

script_dir = os.path.dirname(os.path.realpath(__file__))
if __name__ == '__main__':
    sys.path.append(os.path.join(script_dir, '..'))

from fragment_SELFIES import process_selfies,selfies_split,token_classify
duration=time.perf_counter()-global_start

def block_quality_control(block,max_length):

    ### assess if nb of tokens below or equal  to threshold 
    if len(selfies_split(block)) > max_length:
        return False
    else:
        
        ### assess insertability 
        new_block="[C]"+block+"[C]"
        smi=sf.decoder(new_block)
        temp_block=sf.encoder(smi)
        if new_block==temp_block:
            
            ### assess generability 
            mol=Chem.MolFromSmiles(smi)
            sio = sys.stderr = StringIO()
            AllChem.EmbedMolecule(mol)
            error_log =sio.getvalue()
            if error_log:
                return False
            else:
                
                return True
        else:
            return False



def assess_blocks(smiles,max_length):
    
    ### try to get blocks
    selfies=sf.encoder(smiles)
    try :
        blocks=process_selfies(selfies)
    except:
        return [0,None]
    
    ### test if blocks are good 
    for block in blocks:
        if block_quality_control(block,max_length):
            pass
        else:
            return [0,None]
  
    return [len(blocks),blocks]
 
    

def smiles2blocks(smiles,params):
    
    ### instantiate start time + metrics dictionnary
    start= time.perf_counter()
    metrics_dict= {
        'smiles': smiles,
        'nb_blocks': 0,
        'combination_search_time' : float(),
        'blocks_assement_time' : float(),
        'total_time' : float()
    }

    ### find combinaison 
    temp_start = time.perf_counter() 
    mol= Chem.MolFromSmiles(smiles)
    Chem.Kekulize(mol)
    unique_SMILES_list=list(set(Chem.MolToRandomSmilesVect(mol,params["nb_random_SMILES"])))
    metrics_dict['combination_search_time']= time.perf_counter() - temp_start
    
    ### assess block quality 
    temp_start = time.perf_counter()
    results= [assess_blocks(smi,params['nb_max_tokens']) for smi in unique_SMILES_list]
    metrics_dict['blocks_assement_time']= time.perf_counter() - temp_start

    ### unpack and select best block 
    nb_blocks_list,blocks_list= zip(*results)
    nb_blocks_list= np.array(nb_blocks_list)
    good_block= blocks_list[nb_blocks_list.argmax()]
    if good_block:
        metrics_dict['nb_blocks']= len(good_block)
        good_block= np.array(good_block)
    else:
        good_block=  np.array(["[nop]"])
    metrics_dict['total_time']= time.perf_counter() - start    
    
    return [good_block,metrics_dict]
    
def main (params):

    ### create time dictionnary
    filename=params['input_path'].split("/")[-1]
    times_dict={
        'filename' : filename,
        'libraries_loading': params["import_time"]
        }

    ### load data 
    temp_start = time.perf_counter()
    subset=pd.read_csv(params['input_path'])
    times_dict['data_loading']= time.perf_counter()-temp_start

    ### generate dir for results and input relocation 
    temp_start = time.perf_counter()
    params['dir_path']=params['input_path'][:-4]+"/"
    os.makedirs(params['dir_path'],exist_ok=True)
    params['new_location']=os.path.join(params['dir_path'],filename) 
    shutil.move(params['input_path'],params['new_location'])
    times_dict['files_management']= time.perf_counter()-temp_start

    ### generate block
    temp_start = time.perf_counter()
    with mp.Pool(params['processes']) as pool:
        results=pool.map(partial(smiles2blocks,params=params),subset['sanitized_smiles'])
        pool.close()
        pool.join()
    times_dict['multiprocessing']= time.perf_counter()-temp_start

    ### repacking results  and timelog
    temp_start = time.perf_counter()
    blocks_list,metric_dict_list= zip(*results)
    metric_df= pd.DataFrame.from_dict(metric_dict_list)
    df_total= pd.concat([subset,metric_df],axis=1)

    blocks_list= np.concatenate(blocks_list)
    blocks_set= set(blocks_list)
    if "[nop]" in blocks_set:
        blocks_set.remove("[nop]")
    else: 
        pass
    blocks_df= pd.DataFrame(list(blocks_set),columns=["blocks"])
    
    times_dict['results_repacking']= time.perf_counter()-temp_start
    times_dict['total_time']= time.perf_counter()- params['start_time']
    timelog_df= pd.DataFrame(times_dict,index=[0])

    ### write results and timelog 
    timelog_df.to_csv(os.path.join(params['dir_path'],f"timelog_{filename}"),index=False)
    df_total.to_csv(os.path.join(params['dir_path'],f"metrics_{filename}"),index=False)
    blocks_df.to_csv(os.path.join(params['dir_path'],f"blocks_{filename}"),index=False)

if __name__ == '__main__':
    params={
      'input_path' : str(sys.argv[1]),
      'processes' : int(sys.argv[2]),
      'nb_random_SMILES': int(sys.argv[3]),
      'nb_max_tokens': int(sys.argv[4]),
      'import_time' : duration,
      'start_time' : global_start
    }
    main(params)
    

    
