import time
start = time.perf_counter()
import pandas as pd
import json
import os
import sys
from rdkit import Chem
from rdkit.Chem import RDConfig
sys.path.append(os.path.join(RDConfig.RDContribDir, 'SA_Score'))
import sascorer
from rdkit.Chem import QED,AllChem
import multiprocessing as mp
from functools import partial
duration=time.perf_counter()-start
print(f'It took {duration:.2f}s to load libraries')

### set working directory
script_dir = os.path.dirname(os.path.realpath(__file__))
if __name__ == '__main__':
    sys.path.append(os.path.join(script_dir, '..'))
from parsers import sanitize_smiles_parser


def compute_neutralisation_reaction(patterns_pair_list):
    return [(Chem.MolFromSmarts(x), Chem.MolFromSmiles(y, False)) for x, y in patterns_pair_list]


def NeutraliseCharges(smiles, reactions):

    ### instantiate mol object and initiate variable replaced
    mol = Chem.MolFromSmiles(smiles)
    replaced = False

    ### check if anything is to be neutralized 
    for reactant, product in reactions:
        while mol.HasSubstructMatch(reactant):
            replaced = True
            rms = AllChem.ReplaceSubstructs(mol, reactant, product)
            mol = rms[0]
    
    ### return new smiles if something changed or else return input smiles
    if replaced:
        return Chem.MolToSmiles(mol, True)
    else:
        return smiles


def smiles_quality_control(smiles,max_ringsize=6,SA_threshold=4.5,QED_threshold=0.5):
    
    ### compte 2D strucutre
    mol=Chem.MolFromSmiles(smiles,)

    ### compute metrics 
    sa_score = sascorer.calculateScore(mol)
    qed_score=QED.qed(mol)
    ssr=Chem.GetSymmSSSR(mol)

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
    if sa_score <= SA_threshold and qed_score >= QED_threshold  and mol_max_ringsize <= max_ringsize:
        return [Chem.MolToSmiles(mol,isomericSmiles=False,kekuleSmiles=True),number_rings,mol_max_ringsize,sa_score,qed_score,True]
    else:
        return [Chem.MolToSmiles(mol,isomericSmiles=False,kekuleSmiles=True),number_rings,mol_max_ringsize,sa_score,qed_score,False]
          
def sanitize_smiles(smiles,params,reactions):

    ### remove the charge from SMILES
    neutral_smiles = NeutraliseCharges(smiles,reactions)

    ### assess smiles quality 
    return smiles_quality_control(neutral_smiles,max_ringsize=params['max_ringsize'],SA_threshold=params['SA_threshold'],QED_threshold=params['QED_threshold'])  

def main(params):

    ### initiate metrics dict
    metrics_dict={}

    ###loading smiles from csv  
    start = time.perf_counter()
    dataset = pd.read_csv(params['input_libpath'])
    metrics_dict["data_loading_time"]=time.perf_counter()-start

    ### compute reactions for Neutralisation
    start = time.perf_counter()
    with open(params['input_pattern_path'],'r') as json_file:
        chemical_patterns_pair_list= json.load(json_file)
        json_file.close()
    reactions= compute_neutralisation_reaction(chemical_patterns_pair_list)
    metrics_dict["reaction_init_time"]=time.perf_counter()-start
    
    ### clean the smiles and convert them to selfies 
    start = time.perf_counter()
    with mp.Pool(mp.cpu_count()) as pool:
        results = pool.map(partial(sanitize_smiles,params=params,reactions=reactions), dataset['smiles'])
        pool.close()
        pool.join()
    metrics_dict["quality_control_time"]=time.perf_counter()-start

    ### aggregate results as a Dataframe 
    start = time.perf_counter()
    df= pd.DataFrame(results,columns=["sanitized_smiles","number_rings","max_ringsize","SA_score","QED_score","SMILES_status"])
    df= pd.concat([dataset,df], axis=1)
    original_len= len(df)
    metrics_dict["results_aggregation_time"]=time.perf_counter()-start
    
    ### write results as csv files
    start = time.perf_counter()
    filename= params['input_libpath'].split("/")[-1]
    output_path= os.path.join(params["output_dir"],filename.replace(".csv","_metrics.csv"))
    df.to_csv(output_path,index=False)
    df= df[df["SMILES_status"]==True]
    metrics_dict["perct_good_mols"]= len(df) /original_len *100
    output_path= os.path.join(params["output_dir"],filename.replace(".csv","_sanitized.csv"))
    df.to_csv(output_path,columns=["sanitized_smiles"],index=False)
    metrics_dict["write_results_time"]=time.perf_counter()-start
    
    ### write metrics dictionnary as csv file 
    output_path=os.path.join(params["output_dir"],filename.replace(".csv","_experience_metrics.csv"))
    metrics_dict["total_time"]= time.perf_counter() - params["start_time"]
    df= pd.DataFrame(metrics_dict,index=[0])
    df.to_csv(output_path,index=False)

if __name__ == '__main__':
    params=sanitize_smiles_parser()
    params["libraries_import_time"]= duration
    params["start_time"]=start
    main(params)