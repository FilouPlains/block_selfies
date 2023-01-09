import time 
start = time.perf_counter()
import pandas as pd 
import sys
from io import StringIO
import os
import time
from pathlib import Path
import numpy as np
import multiprocessing as mp
import subprocess
from functools import partial
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import RDConfig
from rdkit.Chem.Pharm3D import Pharmacophore, EmbedLib
from rdkit.Chem.rdFMCS import FindMCS
### change the working directory 
script_dir = os.path.dirname(os.path.realpath(__file__))
if __name__ == '__main__':
    sys.path.append(os.path.join(script_dir, '..'))
    print(os.getcwd())

from parsers import gen_3D_struct_parser

duration=time.perf_counter()-start

def write_mol(path,mol):
    with open(path,'w+') as file:
        file.write(Chem.MolToMolBlock(mol))
        file.close()

def molfile2pdbqtfile(path,pH=7):
    new_path= path.replace(".mol",".pdbqt")
    subprocess.run(f"obabel {path} -O {new_path} -p {pH} --partialcharge gasteiger".split())
    os.remove(path)

def assess_pharmacophore_match(mol,smiles_dict,template_3D,feature_factory,pcophore):

    ### initiate a new metrics dict, prevent modifying original dict
    metrics_dict= {
        "original_smiles": smiles_dict["original_smiles"],
        "selected_smiles": smiles_dict["selected_smiles"],
        "perct_core_template_atom_match" : 0.,
        "perct_mol_template_atom_match" : 0.,
        "core_mol_RMSD" : np.nan,
        "can_match_pharmacophore" : False,
        "number_matches" : 0,
        "selected" : False
    }

    ### assess if 2D structure exists 
    if not mol :
        return None,metrics_dict
    else:
        pass

    ### generate the core of mol to be constructed 
    Chem.Kekulize(mol)
    mol_pair=[mol,template_3D]
    res= FindMCS(mol_pair,threshold=0.9,completeRingsOnly=True)
    p = Chem.MolFromSmarts(res.smartsString)
    core=AllChem.DeleteSubstructs(AllChem.ReplaceSidechains(template_3D,p),Chem.MolFromSmiles('*'))

    ### assess if a core was succesfully created
    if not core :
        return None,metrics_dict
    else:
        pass

    ### sanitize core for construction 
    try:
        Chem.Kekulize(core)
    except:
        return None,metrics_dict
    
    ### get metric from the core 
    metrics_dict["perct_core_template_atom_match"]=core.GetNumAtoms()/ template_3D.GetNumAtoms() * 100
    metrics_dict["perct_mol_template_atom_match"]=core.GetNumAtoms()/ mol.GetNumAtoms() * 100
    
    
    ### generate 3D structure and then constrained it with core 
    try :
        ### check for warning in structure generation prevent use of aberrant structure 
        sio = sys.stderr = StringIO()
        AllChem.EmbedMolecule(mol)
        error_log =sio.getvalue()
        if error_log: 
            raise Exception
        else:
            AllChem.ConstrainedEmbed(mol,core)
        metrics_dict["core_mol_RMSD"]= float(mol.GetProp('EmbedRMS'))  
    except:
        return None,metrics_dict
    
    ### check if it can be matched to pharmacophore 
    canMatch,allMatches = EmbedLib.MatchPharmacophoreToMol(mol,feature_factory,pcophore)
    metrics_dict["can_match_pharmacophore"]= canMatch
    if canMatch:
        metrics_dict["number_matches"]= len(allMatches)
        return mol,metrics_dict
    else:
        return None,metrics_dict
    
def structure_3D_pipeline(pair,template_3D,pcophore):

    ### initiate variables 
    feature_factory = AllChem.BuildFeatureFactory(str(Path(RDConfig.RDDataDir) / "BaseFeatures.fdef")) # not optimal to compute each time but doesn't work with partial
    smiles= pair[0]
    path = pair[1]
    smiles_dict={
        "original_smiles": smiles,
        "selected_smiles" : smiles
        }
    
    ### compute number of stereocenter(s)
    mol=Chem.MolFromSmiles(smiles)
    nb_chiral= Chem.rdMolDescriptors.CalcNumAtomStereoCenters(mol)

    ### check if multiple isomers may exists 
    if nb_chiral >0:
        ### initiate lists 
        mols_list=[]
        metrics_dicts_list=[]

        ### compute isomers 2D structure 
        opts = Chem.EnumerateStereoisomers.StereoEnumerationOptions(unique=True)
        isomers = list(Chem.EnumerateStereoisomers.EnumerateStereoisomers(mol, options=opts))
        
        ### assess  for each isomers if a match to pharmacopore is possible 
        for isomer in isomers:
            smiles_dict["selected_smiles"]= Chem.MolToSmiles(isomer,canonical=False) # keep smiles as close to origin smiles as it is possible
            temp_mol,metrics_dict= assess_pharmacophore_match(isomer,smiles_dict,template_3D,feature_factory,pcophore)
            mols_list.append(temp_mol)
            metrics_dicts_list.append(metrics_dict)

        ### create dataframe from lists of dicts
        df= pd.DataFrame.from_dict(metrics_dicts_list)

        ### find mol with most matches if at least one mol match 
        if df["can_match_pharmacophore"].any() :
            sel_index= df["number_matches"].idxmax()
            df["selected"][sel_index]= True
            write_mol(path,mols_list[sel_index])
        else:
            pass
        
    else:
        ### assess if a match to pharmacopore is possible
        temp_mol,metrics_dict= assess_pharmacophore_match(mol,smiles_dict,template_3D,feature_factory,pcophore)

        ### create dataframe from dict 
        df= pd.DataFrame(metrics_dict,index=[0])
        
        ### check if there is a match 
        if df["can_match_pharmacophore"][0] :
            df["selected"][0]= True
            write_mol(path,temp_mol)
        else:
            pass
    
    ### write metrics dataframe 
    path= path.replace(".mol","_metrics.csv").replace("structures_files","metrics")
    df.to_csv(path,index=False)



def main(params):

    ### initiate metrics dictionnary 
    metrics_dict= {"import_time": params["import_time"]}

    ### create output directories
    start = time.perf_counter()
    dir_path= params["input_csvpath"].replace(params["input_csvpath"].split("/")[-1],"")
    struct_path= os.path.join(dir_path,"structures_files/")
    os.makedirs(struct_path,exist_ok=True)
    metrics_path= os.path.join(dir_path,"metrics/")
    os.makedirs(metrics_path,exist_ok=True)
    gnu_instructions_path=os.path.join(dir_path,"GNU_instructions/")
    os.makedirs(gnu_instructions_path,exist_ok=True)
    metrics_dict["dir_creation_time"]=time.perf_counter()-start

    ### load dataframe 
    start = time.perf_counter()
    df=pd.read_csv(params["input_csvpath"])
    metrics_dict["csv_loading_time"]=time.perf_counter()-start

    ###  load 3D template structure from pdb and sanitize it 
    start = time.perf_counter()
    start_mol=Chem.MolFromSmiles(params["input_smiles"])
    template_3D= Chem.MolFromPDBFile(params["input_pdbpath"])
    template_3D=AllChem.AssignBondOrdersFromTemplate(start_mol,template_3D)
    Chem.RemoveHs(template_3D)
    Chem.Kekulize(template_3D)
    metrics_dict["3D_template_time"]=time.perf_counter()-start

    ###  make pharmacophore from 3D template 
    start = time.perf_counter()
    feature_factory = AllChem.BuildFeatureFactory(str(Path(RDConfig.RDDataDir) / "BaseFeatures.fdef"))
    features = feature_factory.GetFeaturesForMol(template_3D)
    pcophore= Pharmacophore.Pharmacophore(features)
    metrics_dict["pharmacophore_gen_time"]=time.perf_counter()-start

    ### 3D struct generation and pharmacophore filtering
    start = time.perf_counter()
    path_list= [os.path.join(struct_path,f"mol_{index}.mol") for index in df["index"]]
    pair_list= zip(df["unique_smiles"],path_list)
    with mp.Pool(mp.cpu_count()) as pool:
        pool.map(partial(structure_3D_pipeline,
        template_3D=template_3D,
        pcophore= pcophore),
        pair_list)
        pool.close()
        pool.join()
    metrics_dict["structures_gen_time"]=time.perf_counter()-start

    ### convert mol file to pdbqt format 
    start = time.perf_counter()
    path_list=[path for path in path_list if os.path.exists(path)]
    for path in path_list:
        molfile2pdbqtfile(path,params["pH"])
    metrics_dict["mol2pdbqt_time"]=time.perf_counter()-start

    ### generate GNU instruction 
    filename= params['input_csvpath'].split("/")[-1].replace("subset","task").replace(".csv",".txt")
    instructions_file_path=os.path.join(gnu_instructions_path,filename)
    with open(instructions_file_path,"w") as instructions_file:
        for path in path_list:
            instructions_file.write(f"python scripts/docking_pipeline.py --ligand_path {path.replace('.mol','.pdbqt')}\n")
        instructions_file.close()
    
    ### generate slurm file if it does not exist
    slurm_file_path=os.path.join(script_dir,"../slurm_docking.sh")
    os.makedirs("logs/parallel_jobs/", exist_ok=True)
    if not os.path.exists(slurm_file_path):

        ### compute time allocation + add 20% safety margin and convert to minutes
        mem2allocate_per_task= params['cores']*256 
        duration= (len(df)/params['cores'])*params['dt_mol'] + params['dt_lib']
        duration += 0.5*duration
        duration= int(np.ceil(duration/60))

        with open(slurm_file_path,"w") as slurm_file:
            slurm_file.write("#!/bin/sh\n")
            slurm_file.write("#SBATCH --account=def-jeromew\n")
            slurm_file.write(f"#SBATCH --time=00:{duration:02d}:00\n")
            slurm_file.write("#SBATCH --job-name=GNU_docking\n")
            slurm_file.write("#SBATCH --output=/home/retienne/projects/def-jeromew/retienne/blocks_SELFIES/logs/GNU_docking_%A.out\n")
            slurm_file.write("#SBATCH --error=/home/retienne/projects/def-jeromew/retienne/blocks_SELFIES/logs/GNU_docking_%A.err\n")
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
            slurm_instruction="parallel --jobs $SLURM_CPUS_PER_TASK --joblog logs/parallel_jobs/tasks_${SLURM_ARRAY_TASK_ID}.log < " + f"{gnu_instructions_path}" + "task_${SLURM_ARRAY_TASK_ID}.txt\n"
            slurm_file.write(slurm_instruction)
            slurm_file.close()
    else:
        pass

if __name__ == '__main__':
    params=gen_3D_struct_parser()
    params["start_time"]= start
    params["import_time"]=duration
    main(params)




