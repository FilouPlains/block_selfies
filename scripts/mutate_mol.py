
import time 
start = time.perf_counter()
import json
import pandas as pd
import numpy as np
import selfies as sf 
import os
import sys
import multiprocessing as mp
from functools import partial
from numba import njit,prange,int64,vectorize
from rdkit import Chem
from rdkit.Chem import RDConfig
sys.path.append(os.path.join(RDConfig.RDContribDir, 'SA_Score'))
# now you can import sascore!
import sascorer
from rdkit.Chem import QED,AllChem
from rdkit.Chem import MolToInchi
from rdkit.DataStructs.cDataStructs import TanimotoSimilarity


### change the working directory 
script_dir = os.path.dirname(os.path.realpath(__file__))
if __name__ == '__main__':
    sys.path.append(os.path.join(script_dir, '..'))
    print(os.getcwd())

from fragment_SELFIES import process_selfies,selfies_split
from parsers import mutate_mol_parser
duration=time.perf_counter()-start

def blocks2intarray(blocks, vocab_btoi):
    return np.array([vocab_btoi[block] for block in blocks])

def intarray2selfies(int_array,vocab_itob):
    return "".join([vocab_itob[int_block] for int_block in int_array])

def make_selfies_smi_pair(int_array,vocab_itob):
    selfies= intarray2selfies(int_array,vocab_itob)
    return [selfies,sf.decoder(selfies)]

def gen_starting_blocks(blocks,vocab_btoi,int_pad):
    int_array=  blocks2intarray(blocks,vocab_btoi)
    starting_blocks_list=np.tile(int_array,(len(int_array)+1,1))
    for i in range(len(int_array)):
        starting_blocks_list[i,i]=int_pad
    return starting_blocks_list

@vectorize([int64(int64,int64)])
def len_resize(int_len,choice):
    if choice==0:
        int_len+=1
    else:
        pass
    return int_len

@njit
def mutate_selfies(starting_blocks,choice,position,int_block,padding_int):
    ### insert block
    if choice==0:
        return np.concatenate((starting_blocks[:position],np.array([int_block]),starting_blocks[position:]),axis=0)
    ### replace block
    elif choice ==1:
        duplicated_blocks=np.copy(starting_blocks)
        duplicated_blocks[position]=int_block
        return np.concatenate((duplicated_blocks,np.array([padding_int])),axis=0)
    else:
        return None

@njit(parallel=True)
def multiprocess_mutate(nb_mutation,starting_blocks,position_list,triplet_list,padding_int):
    results=np.zeros((nb_mutation,starting_blocks.shape[1]+1),dtype=np.int64)
    for i in prange(nb_mutation):
        results[i]=mutate_selfies(starting_blocks[triplet_list[i,0],:],triplet_list[i][1],position_list[i],triplet_list[i][2],padding_int)
    return results

def asses_stability(int_array,vocab_itob):
    original_blocks= [vocab_itob[int_block] for int_block in int_array if vocab_itob[int_block]!= "[nop]"]
    selfies= "".join(original_blocks)
    smiles=  sf.decoder(selfies)
    if smiles :
        stable_selfies= sf.encoder(smiles)
        if stable_selfies == selfies:
            return [selfies,stable_selfies,smiles,"stable",0.]
        else:
            stable_tokens= selfies_split(stable_selfies)
            original_tokens= selfies_split(selfies)
            length_original=len(original_tokens)
            length_stable=len(stable_tokens)
            if length_stable == length_original:
                return [selfies,stable_selfies,smiles,"rescued",0.]
            else:
                return [selfies,stable_selfies,smiles,"truncated",np.absolute(length_stable-length_original)/length_original*100]
    else:
        return [selfies,None,None,"failure",100.]
    


def mol_quality_control(smiles,target_fps,max_chiral=3,max_ringsize=6,SA_threshold=4.5,QED_threshold=0.5,tanimoto_lower=0.2,tanimoto_upper=0.8):
    
    ### compte 2D strucutre
    mol=Chem.MolFromSmiles(smiles)

    ### compute metrics 
    sa_score = sascorer.calculateScore(mol)
    qed_score=QED.qed(mol)
    nb_chiral=Chem.rdMolDescriptors.CalcNumAtomStereoCenters(mol)
    ssr=Chem.GetSymmSSSR(mol)
    fps=AllChem.GetMorganFingerprint(mol, 2)
    tanimoto_coef=TanimotoSimilarity(fps,target_fps)

    ### check ring size 
    number_rings=0
    mol_max_ringsize=0
    ringsize_list=[]
    if ssr:
        for ring in ssr:
            ringsize_list.append(len(list(ring)))
            number_rings+=1
    else:
        pass
    if  ringsize_list:
        mol_max_ringsize= max(ringsize_list)
    else:
        pass

    ### check if all conditions are met 
    if (sa_score <= SA_threshold 
        and qed_score >= QED_threshold 
        and nb_chiral <= max_chiral 
        and mol_max_ringsize <= max_ringsize
        and tanimoto_coef >= tanimoto_lower
        and tanimoto_coef <= tanimoto_upper):
        return [MolToInchi(mol),nb_chiral,number_rings,mol_max_ringsize,sa_score,qed_score,tanimoto_coef,True]
    else:
        return [MolToInchi(mol),nb_chiral,number_rings,mol_max_ringsize,sa_score,qed_score,tanimoto_coef,False]


def assess_uniqueness(mol_pair_list):

    ### initate molecule dictionnary : InChI are keys, smiles are values 
    mol_dict={}

    ### iterate over pair list (InChI, SMILES)
    for inchi,smi in mol_pair_list:

        ### if inchi already in dict then choose SMILES that yields the most blocks 
        if inchi in mol_dict:
            recorded_selfies=sf.encoder(mol_dict[inchi])
            challenger_selfies= sf.encoder(smi)
            try:
                blocks_len_recorded=len(process_selfies(recorded_selfies))
                blocks_len_challenger=len(process_selfies(challenger_selfies))
            except:
                continue
            if blocks_len_challenger> blocks_len_recorded:
                mol_dict[inchi]=smi
            else:
                pass

        ### just add the pair if not already present in dict 
        else:
            mol_dict[inchi]=smi

    return mol_dict

def make_subset(subset_id,df,params):

    ### comprehensive slicing of dataframe 
    chunk_min, chunk_max = subset_id * params['chunck_size'], min((subset_id + 1) * params['chunck_size']-1, params['n']-1)
    subset = df.loc[chunk_min:chunk_max]
    if subset_id < params['nb_unatributed_smiles']:
        index2add= params['n'] -(subset_id+1)
        subset.loc[index2add] = df.loc[index2add]
    
    ### write subset
    susbset_file_path=os.path.join(params["output_structpath"],f"subset_{subset_id}.csv")
    subset.to_csv(susbset_file_path,index=False)

def main(params):

    ### initiate metrics dictionnary and create output dir if it doesn't already exist
    metrics_dict= {"import_time": params["import_time"]}
    os.makedirs(params["output_dirpath"],exist_ok=True)

    ### load fragment library & starting blocks 
    start = time.perf_counter()
    df=pd.read_csv(params["input_libpath"])
    blocks_set= set(df["blocks"])
    with open(params["input_blokcspath"],"r") as json_file:
        best_blocks= json.load(json_file)
        json_file.close()
    start_mol= Chem.MolFromSmiles(sf.decoder("".join(best_blocks)))
    metrics_dict["loading_time"]=time.perf_counter()-start

    ### update set and generate dicts
    start= time.perf_counter() 
    for block in best_blocks :
        blocks_set.add(block)
    vocab_btoi={unique_block:i for i,unique_block in enumerate(blocks_set)}
    vocab_itob= {v: k for k, v in vocab_btoi.items()}
    vocab_btoi["[nop]"]=len(vocab_btoi)
    vocab_itob[len(vocab_itob)]= "[nop]"
    metrics_dict["gen_vocabs_time"]=time.perf_counter()-start

    ### generate starting_blocks
    start = time.perf_counter()
    starting_blocks=gen_starting_blocks(best_blocks,vocab_btoi,len(vocab_btoi)-1)
    with mp.Pool(mp.cpu_count()) as pool:
        selfies_smi_pair_list=pool.map(partial(make_selfies_smi_pair,vocab_itob=vocab_itob),starting_blocks)
        pool.close()
        pool.join()
    df=pd.DataFrame(selfies_smi_pair_list,columns=["selfies","smiles"])
    output_path=os.path.join(params["output_dirpath"],f'{params["name"]}_starting_blocks.csv')
    df.to_csv(output_path,index=False)
    metrics_dict["starting_block_gen_time"]=time.perf_counter()-start

    ### randomly select starting blocks,action(insert or replace),block from block library
    start = time.perf_counter()
    triplet_list=np.random.randint(0,[starting_blocks.shape[0],2,len(vocab_btoi)-1],size=(params["nb_mutations"],3))
    position_max_value_list= len_resize(np.repeat(starting_blocks.shape[1],params["nb_mutations"]),triplet_list[:,1])
    position_list=np.random.randint(0,position_max_value_list)
    metrics_dict["random_sel_time"]=time.perf_counter()-start

    ### mutation process
    start = time.perf_counter()
    mutated_blocks_list= multiprocess_mutate(params["nb_mutations"],starting_blocks,position_list,triplet_list,len(vocab_btoi)-1)
    metrics_dict["mutation_time"]=time.perf_counter()-start
    
    ### filter out truncated selfies
    start = time.perf_counter()
    with mp.Pool(mp.cpu_count()) as pool:
        results=pool.map(partial(asses_stability,vocab_itob=vocab_itob),mutated_blocks_list)
        pool.close()
        pool.join()   
    df=pd.DataFrame(results,columns=["selfies","stable_selfies","smiles","stability_status","truncation_rate"])
    output_path=os.path.join(params["output_dirpath"],f'{params["name"]}_stability_metrics.csv')
    df.to_csv(output_path,index=False)
    stable_subset=df[df["stability_status"]=="stable"]
    stable_subset= stable_subset.reset_index()
    metrics_dict["mean_trunc_rate"]= np.mean(df["truncation_rate"])
    metrics_dict["truncation_filter_time"]=time.perf_counter()-start

    ### filter out non drug like molecules 
    start = time.perf_counter()
    target_fps=AllChem.GetMorganFingerprint(start_mol, 2)
    with mp.Pool(mp.cpu_count()) as pool :
        results=pool.map(partial(mol_quality_control,
        max_chiral= params["max_nb_chiral"],
        max_ringsize=params["max_ringsize"],
        SA_threshold=params["SA_threshold"],
        QED_threshold=params["QED_threshold"],
        tanimoto_lower=params["tanimoto_lower_bound"],
        tanimoto_upper=params["tanimoto_upper_bound"],
        target_fps=target_fps),
        stable_subset["smiles"])
        pool.close()
        pool.join()
    
    df= pd.DataFrame(results,columns=["InChI","nb_chiral","nb_rings","max_ring_size","SA_score","QED_score","tanimoto_coef","status"])
    df= pd.concat([stable_subset,df],axis=1)
    output_path=os.path.join(params["output_dirpath"],f'{params["name"]}_quality_metrics.csv')
    df.to_csv(output_path,index=False)
    druglike_subset= df[df["status"]==True]
    metrics_dict["perct_druglike_mol"]= len(druglike_subset)/params["nb_mutations"]*100
    metrics_dict["quality_control_time"]=time.perf_counter()-start

    ### assess mol uniqueness 
    start = time.perf_counter()
    unique_mols_dict=assess_uniqueness(list(zip(druglike_subset["InChI"],druglike_subset["smiles"])))
    df=pd.DataFrame(unique_mols_dict.values(),columns=["unique_smiles"])
    df= df.reset_index()
    output_path= os.path.join(params["output_dirpath"],f'{params["name"]}_unique_smiles.csv')
    df.to_csv(output_path,index=False)
    metrics_dict["perct_unique_mol"]= len(df)/len(druglike_subset)*100
    metrics_dict["uniqeness_control_time"]=time.perf_counter()-start

    ### dataset split for structure generation 
    params['n']=len(df)
    params['chunck_size']= params['n'] // params['array_task']
    params['nb_unatributed_smiles']= params['n'] % params['chunck_size']
    params["output_structpath"]= os.path.join(params["output_dirpath"],f"{params['name']}_3D_structure")
    os.makedirs(params["output_structpath"],exist_ok=True)
    subset_id_list=range(params['array_task'])
    with mp.Pool(mp.cpu_count()) as pool:
        pool.map(partial(make_subset,df=df,params=params), subset_id_list)
        pool.close()
        pool.join()
    
    
    ### compute time allocation + add 20% safety margin and convert to minutes  
    mem2allocate_per_task= params['cores']*256 
    duration= params['chunck_size'] * params['dt_smiles'] + params['dt_lib']
    duration += 0.5*duration
    duration= int(np.ceil(duration/60))

    ### create logs dir 
    os.makedirs("logs/gen_3D_struct/output_logs/",exist_ok=True)
    os.makedirs("logs/gen_3D_struct/error_logs/",exist_ok=True)

    ### write slurm file 
    slurm_file_path=os.path.join(script_dir,"../slurm_gen_3D_struct.sh")
    with open(slurm_file_path,"w") as slurm_file:
        slurm_file.write("#!/bin/sh\n")
        slurm_file.write("#SBATCH --account=def-jeromew\n")
        slurm_file.write(f"#SBATCH --time=00:{duration:02d}:00\n")
        slurm_file.write("#SBATCH --job-name=blocks_generation\n")
        slurm_file.write("#SBATCH --output=/home/retienne/projects/def-jeromew/retienne/blocks_SELFIES/logs/gen_3D_struct/output_logs/gen_3D_struct_%A_%a.out\n")
        slurm_file.write("#SBATCH --error=/home/retienne/projects/def-jeromew/retienne/blocks_SELFIES/logs/gen_3D_struct/error_logs/gen_3D_struct_%A_%a.err\n")
        slurm_file.write(f"#SBATCH --cpus-per-task={params['cores']}\n")
        slurm_file.write(f"#SBATCH --mem={mem2allocate_per_task}M\n")
        slurm_file.write(f"#SBATCH --array=0-{params['array_task']-1}\n")
        slurm_file.write("\n")
        slurm_file.write("### make sure good module are loaded\n")
        slurm_file.write("module purge\n")
        slurm_file.write("module restore Blocks_modules\n")
        slurm_file.write("\n")
        slurm_file.write("### keep track of modules and python packages\n")
        slurm_file.write("module list\n")
        slurm_file.write("pip list\n")
        slurm_file.write("\n")
        slurm_file.write("### skip line for log cleanliness\n")
        slurm_file.write("echo\n")
        slurm_file.write("\n")
        slurm_file.write("### launch script\n")
        instruction="python scripts/gen_3D_struct.py " + f"{params['output_dirpath']}" + "subset_${SLURM_ARRAY_TASK_ID}.csv \n"
        slurm_file.write(instruction)
        slurm_file.close()
    
    ### write metrics
    metrics_dict["total_time"]=   time.perf_counter() - params["start_time"]
    df=pd.DataFrame(metrics_dict,index=[0])
    output_path= os.path.join(params["output_dirpath"],f'{params["name"]}_exp_metrics.csv')
    df.to_csv(output_path,index=False)

if __name__ == '__main__':
    params=mutate_mol_parser()
    params["start_time"]= start
    params["import_time"]=duration
    main(params)