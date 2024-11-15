import time 
time1 = time.perf_counter()
import pandas as pd
import os
import sys
import math
from multiprocessing import Pool
from functools import partial
duration=time.perf_counter()-time1
print(f'It took {duration:.2f}s to load python libraries')

### set working directory
script_dir = os.path.dirname(os.path.realpath(__file__))
if __name__ == '__main__':
    sys.path.append(os.path.join(script_dir, '..'))
from parsers import prep_blocks_generation_parser


def make_subset(subset_id,smiles_list,params):
     ### comprehensive slicing of SMILES for subset attribution 
        chunk_min, chunk_max = subset_id * params['chunck_size'], min((subset_id + 1) * params['chunck_size'], params['n'])
        list_data = smiles_list[chunk_min:chunk_max]
        if subset_id < params['nb_unatributed_smiles']:
            list_data.append(smiles_list[-(subset_id+1)])
        
        ### subset generation 
        subset=pd.DataFrame(list_data).rename(columns={0: 'sanitized_smiles'})
        susbset_file_path=os.path.join(params['output_dirpath'],f"subset_{subset_id}.csv")
        subset.to_csv(susbset_file_path,index=False)



def main(params):

    ###loading sanitized smiles from csv  
    dataset = pd.read_csv(params['input_path'])
    smiles_list= list(dataset['sanitized_smiles'])
    print("Data Loading Successful")

    ### compute slurm parameters
    params['n']=len(dataset)
    mem2allocate_per_task= params['cores']*256 
    params['chunck_size']= params['n'] // params['array_task']
    params['nb_unatributed_smiles']= params['n'] % params['chunck_size']
    
    ### compute time allocation + add 20% safety margin and convert to minutes  
    duration= params['chunck_size'] * params['dt_smiles'] + params['dt_lib']
    duration += 0.5*duration
    duration= math.ceil(duration/60)

    ### write slurm file 
    slurm_file_path=os.path.join(script_dir,"../slurm_generate_blocks.sh")
    with open(slurm_file_path,"w") as slurm_file:
        slurm_file.write("#!/bin/sh\n")
        slurm_file.write("#SBATCH --account=def-jeromew\n")
        slurm_file.write(f"#SBATCH --time=00:{duration:02d}:00\n")
        slurm_file.write("#SBATCH --job-name=blocks_generation\n")
        slurm_file.write("#SBATCH --output=/home/retienne/projects/def-jeromew/retienne/blocks_SELFIES/logs/generate_blocks/output_logs/generate_blocks_%A_%a.out\n")
        slurm_file.write("#SBATCH --error=/home/retienne/projects/def-jeromew/retienne/blocks_SELFIES/logs/generate_blocks/error_logs/generate_blocks_%A_%a.err\n")
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
        instruction="python scripts/smiles2blocks.py " + f"{params['output_dirpath']}" + "subset_${SLURM_ARRAY_TASK_ID}.csv "+f"{params['cores']} {params['nb_random_SMILES']} {params['nb_max_tokens']}\n"
        slurm_file.write(instruction)
        slurm_file.close()
    
    ### write subset csv file 
    os.makedirs(params['output_dirpath'],exist_ok=True)
    subset_id_list=range(params['array_task'])
    time1 = time.perf_counter()
    pool = Pool(params['processes'])
    pool.map(partial(make_subset,smiles_list=smiles_list,params=params), subset_id_list)
    duration=time.perf_counter()-time1
    pool.close()
    pool.join()
    print(f'It took {duration:.2f}s to generate {params["array_task"] } subsets for {len(smiles_list)} SMILES')    

if __name__ == '__main__':
    params=prep_blocks_generation_parser()
    main(params)
