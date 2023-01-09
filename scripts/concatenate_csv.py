import os
import sys
import glob
import pandas as pd 
import multiprocessing as mp
script_dir = os.path.dirname(os.path.realpath(__file__))
if __name__ == '__main__':
    sys.path.append(os.path.join(script_dir, '..'))

my_path="results/mol_mutations/celecoxib/celecoxib_3D_structure/metrics"
files = glob.glob(my_path + '/*.csv', recursive=True)


def load_df_from_csv(path):
    return pd.read_csv(path)

with mp.Pool(mp.cpu_count()) as pool:
    df_list=pool.map(load_df_from_csv,files)
    pool.close()
    pool.join()

df= pd.concat(df_list,ignore_index= True)

df.to_csv("results/mol_mutations/celecoxib/gen_3D_struct_metrics.csv",index= False)