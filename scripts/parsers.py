import argparse

def parser2dict(parser=argparse.ArgumentParser()):
    args, _ = parser.parse_known_args()
    return vars(args)



def sanitize_smiles_parser():

    ### instantiate  argument parser object 
    parser = argparse.ArgumentParser()

    ### files & dir
    parser.add_argument('--input_libpath', help="the input file should be a csv" ,default="data/whole_MOSES_smiles.csv", type=str)
    parser.add_argument('--input_pattern_path', help="chemical pattern for neutralisation, should be a list in a json",
    default="data/chemical_patterns.json", type=str)
    parser.add_argument('--output_dir',help="the output dir for metrics and results " ,default='data/', type=str)
    
    ### SMILES filtering 
    parser.add_argument("--max_ringsize",help="maximum size of a single ring", type=int,default=6)
    parser.add_argument("--SA_threshold",help="upper bound threshold for Synthetic Accesibility", type=float,default=4.5)
    parser.add_argument("--QED_threshold",help="lower bound threshold for Quantitative Estimate of Drug-Likeness", type=float,default=0.6)
    
    ### convert arguments to dict 
    return  parser2dict(parser=parser)

def prep_blocks_generation_parser():
    parser = argparse.ArgumentParser()   

    ### files & dir  
    parser.add_argument('--input_path', help="the input file should be a csv" ,default="data/whole_MOSES_smiles_sanitized.csv", type=str)
    parser.add_argument('--output_dirpath',help="output dir will be created if it doesn't exist" ,default='data/SMILES_data/subsets_dir/', type=str)
    
    ### multiprocess 
    parser.add_argument('--processes', help="number of processes for multiprocessing", type=int,default=20)

    ### Slurm parameter for job array
    parser.add_argument('--nb_random_SMILES', help="number of random SMILES to try for block generation", type=int,default=2000)
    parser.add_argument('--nb_max_tokens', help="number of maximum tokens per blocks", type=int,default=25)
    parser.add_argument('--cores', help="number of cores used per job in job array", type=int,default=4)
    parser.add_argument('--array_task', help="number of task in job array", type=int,default=4000)
    parser.add_argument('--dt_smiles', help="estimated duration (secondes) to process 1 SMILES", type=float,default=5.)
    parser.add_argument('--dt_lib', help="estimated duration (secondes) to load all libraries", type=int,default=24)


    ### convert arguments to dict 
    args, _ = parser.parse_known_args()
    params=vars(args)
    return params 

def concatenate2lib_parser():
    parser = argparse.ArgumentParser()

    ### files & dir 
    parser.add_argument('--input_dirpath', help="the input path for dir containing the generated subsets folder" 
    ,default="data/SMILES_data/subsets_dir/", type=str)
    parser.add_argument('--output_dirpath', help="the output path for dir containing the metrics" ,default="data/SMILES_data/", type=str)
    parser.add_argument('--output_libpath', help="the output path for the generated library" ,default="data/block_libraires.csv", type=str)
 
    ### multiprocess 
    parser.add_argument('--processes', help="number of processes for multiprocessing", type=int,default=20)

    ### job parameter 
    parser.add_argument('--array_task', help="number of task in job array", type=int,default=4000)

    ### convert arguments to dict 
    args, _ = parser.parse_known_args()
    params=vars(args)
    return params 

def mutate_mol_parser():
    parser = argparse.ArgumentParser()

    ### files & dir 
    parser.add_argument('--input_libpath', help="the input path for the block library" ,default="data/block_libraires.csv", type=str)
    parser.add_argument('--input_blokcspath', help="the input path for the json file of the starting blocks",default="data/celecoxib_starting_blocks.json", type=str)
    parser.add_argument('--output_dirpath', help="the output path for dir containing the metrics and mutated SMILES" 
    ,default="results/mol_mutations/celecoxib/", type=str)

    ### Mutation 
    parser.add_argument('--name', help="experiment name" ,default="celecoxib", type=str) 
    parser.add_argument('--nb_random_SMILES', help="number of random SMILES to try for block generation", type=int,default=10000)
    parser.add_argument('--nb_mutations',help="number of mutation to be tested on input molecule", default=1000000, type=int)

    ### truncation 
    parser.add_argument("--max_trunc_rate",help="maximum truncation rate allowed", type=float,default=0.)

    ### molecular filtering  
    parser.add_argument("--max_ringsize",help="maximum size of a single ring", type=int,default=6)
    parser.add_argument("--max_nb_chiral",help="maximum number of chiral center allowed", type=int,default=3)
    parser.add_argument("--SA_threshold",help="upper bound threshold for Synthetic Accesibility, upper bound is included", type=float,default=4.5)
    parser.add_argument("--QED_threshold",
    help="lower bound threshold for Quantitative Estimate of Drug-Likeness,lower bound is included", type=float,default=0.6)
    parser.add_argument("--tanimoto_lower_bound",help="lower bound for similarity interval, lower bound is included", type=float,default=0.4)
    parser.add_argument("--tanimoto_upper_bound",help="upper bound for similarity interval, upper bound is included", type=float,default=0.6)
    
    ### structure generation prep
    parser.add_argument('--cores', help="number of cores used per job in job array", type=int,default=4)
    parser.add_argument('--array_task', help="number of task in job array", type=int,default=2000)
    parser.add_argument('--dt_smiles', help="estimated duration (secondes) to process 1 SMILES", type=float,default=5.)
    parser.add_argument('--dt_lib', help="estimated duration (secondes) to load all libraries", type=int,default=24)
    
    ### convert arguments to dict 
    return  parser2dict(parser=parser)

def gen_3D_struct_parser():
    parser = argparse.ArgumentParser()

    ### files & dir 
    parser.add_argument('--input_csvpath', help="the input path for the csv cointaining the SMILES" ,
    default="results/mol_mutations/celecoxib/celecoxib_3D_structure/subset_0.csv", type=str)
    parser.add_argument('--input_smiles', help="the SMILES of the template",
    default="CC1=CC=C(C=C1)C2=CC(=NN2C3=CC=C(C=C3)S(=O)(=O)N)C(F)(F)F", type=str)
    parser.add_argument('--input_pdbpath', help="the input path for the pdb file of  3D template",
    default="data/celecoxib_COX2.pdb", type=str)
    parser.add_argument('--output_dirpath', help="the output path for dir containing the metrics and mutated SMILES" 
    ,default="results/mol_mutations/celecoxib/", type=str)

    ### structure generation 
    parser.add_argument('--pH', help="pH for protonation state estimation", type=float,default=7.)

    ### docking prep
    parser.add_argument('--cores', help="number of cores used per job in job array", type=int,default=4)
    parser.add_argument('--array_task', help="number of task in job array", type=int,default=2000)
    parser.add_argument('--dt_mol', help="estimated duration (secondes) to generate structure of 1 molecules", type=float,default=500)
    parser.add_argument('--dt_lib', help="estimated duration (secondes) to load all libraries", type=int,default=24)

    ### convert arguments to dict 
    return  parser2dict(parser=parser)

def docking_parser():
    parser = argparse.ArgumentParser()

    ### files & dir 
    parser.add_argument('--ligand_path', help="path to the ligand to be docked, should be a .pdbqt" ,
    default="results/mol_mutations/celecoxib/celecoxib_3D_structure/structures_files/mol_0.pdbqt", type=str)
    parser.add_argument('--smina_path', help="path to the docking software smina.static" ,
    default="../smina.static", type=str)
    parser.add_argument('--receptor_path', help="path to the receptor used for docking, should be a .pdbqt" ,
    default="data/COX2_chain_A.pdbqt", type=str)
    parser.add_argument('--ref_ligand_path', help="path to the ligand to build the box on, a .pdbqt is expected but mol2 should work too",
    default="data/celecoxib_chain_A.pdbqt", type=str)

    ### SMINA parameters 
    parser.add_argument('--exhaustiveness', help="how exhaustive the search should be for a ligand, smina default is 8", type=float,default=16)
    parser.add_argument('--scoring_function', help="scoring functions to be used for docking", type=str,default="vinardo")
    parser.add_argument('--nb_core', help="number of core that should be attributed for the docking, if using GNU, should always be 1",
    type=int,default=1)
    parser.add_argument('--docking_walltime', help="maximum duration allowed for docking subprocess, in secondes",
    type=int,default=600)
    

    ### convert arguments to dict 
    return  parser2dict(parser=parser)