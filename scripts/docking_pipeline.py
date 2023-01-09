import time 
start = time.perf_counter()
import os 
import pandas as pd
import sys
import numpy as np
import subprocess
import shutil

### change the working directory 
script_dir = os.path.dirname(os.path.realpath(__file__))
if __name__ == '__main__':
    sys.path.append(os.path.join(script_dir, '..'))
from parsers import docking_parser
duration=time.perf_counter()-start

def main(params):

    ### initiate metrics dictionnary 
    metrics_dict={
        'filename': params["ligand_path"].split("/")[-1],
        'import_time': params["duration"],
        'docking_duration': 0.,
        'mean_score' : 0.,
        'docking_status': False
    }
    ### create sub_directory 
    dir_path= params["ligand_path"].replace(".pdbqt","/")
    os.makedirs(dir_path, exist_ok=True)

    ### docking prep
    new_ligand_path= os.path.join(dir_path,metrics_dict["filename"])
    shutil.move(params["ligand_path"],new_ligand_path)
    output_path= new_ligand_path.replace(".pdbqt","_out.pdbqt")
    smina_logpath= new_ligand_path.replace(".pdbqt","_log.txt")
    
    ### docking (make this part as function, when implementing RNA docking)
    try:
        start= time.perf_counter()
        ### prepare and execute subprocess command 
        cmd = f'{params["smina_path"]} --receptor {params["receptor_path"]} --ligand {new_ligand_path}' \
            f' --out {output_path} --autobox_ligand {params["ref_ligand_path"]}'\
                f' --exhaustiveness {params["exhaustiveness"]} --scoring {params["scoring_function"]} '\
                f'--log {smina_logpath} --cpu {params["nb_core"]}'#--quiet
        subprocess.run(cmd.split() ,timeout=params["docking_walltime"]) # ,stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL
        metrics_dict['docking_duration']= time.perf_counter()- start
        
        ### retrieve scores
        with open(output_path, 'r') as f:
            lines = f.readlines()
            slines = [l for l in lines if l.startswith('REMARK minimizedAffinity')]
            values = [l.split() for l in slines]
            ### In each split string, item with index 2 should be the kcal/mol energy.
            score_list = [float(v[2]) for v in values]
            metrics_dict["mean_score"]= np.mean(score_list)
            metrics_dict["docking_status"]= True
    except:
        pass

    ### write metrics as csv file
    output_path=output_path.replace("_out.pdbqt","_metrics.csv")
    metrics_dict["total_time"]= time.perf_counter()-params["start_time"]
    df= pd.DataFrame(metrics_dict,index=[0])
    df.to_csv(output_path,index=False)

if __name__ == '__main__':
    params= docking_parser()
    params["duration"]=duration
    params["start_time"]= start
    main(params)

