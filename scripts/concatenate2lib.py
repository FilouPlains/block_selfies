import time
global_start = time.perf_counter()
from parsers import concatenate2lib_parser
from functools import partial
from multiprocessing import Pool
import shutil
import sys
import os
import pandas as pd
duration = time.perf_counter()-global_start

script_dir = os.path.dirname(os.path.realpath(__file__))
if __name__ == '__main__':
    sys.path.append(os.path.join(script_dir, '..'))



def mkdir(path):
    os.makedirs(path, exist_ok=True)

def file_exist(path):
    return os.path.isfile(path)

def move_file(input_path, output_dir):
    filename = input_path.split("/")[-1]
    new_location = os.path.join(output_dir, filename)
    shutil.move(input_path, new_location)

def assess_launch_status(path, output_dir):
    if  file_exist(path):
        move_file(path, output_dir)
        return False
    else:
        return True

def assess_completion_status(dir_path, output_dir):

    ### instantiate variable
    core_name = dir_path.split("/")[-1]
    blocks_path = os.path.join(dir_path, f"blocks_{core_name}.csv")
    metrics_path = os.path.join(dir_path, f"metrics_{core_name}.csv")
    timelog_path = os.path.join(dir_path, f"timelog_{core_name}.csv")

    ### check if task finished normally
    task_status = file_exist(blocks_path) * \
        file_exist(metrics_path)*file_exist(timelog_path)

    ### move the subset if task not completed
    if not task_status:
        subset_path = os.path.join(dir_path, f"{core_name}.csv")
        move_file(subset_path, output_dir=output_dir)
        return False
    else:
        return True

def make_csv_paths(dir_path,):
    PREFIX_LIST = ["blocks", "metrics", "timelog"]
    return [f"{dir_path}/{prefix}_{dir_path.split('/')[-1]}.csv" for prefix in PREFIX_LIST]

def make_csv_list(dirpath_list):
    return [make_csv_paths(dir_path) for dir_path in dirpath_list]

def csv2pandas(path):
    return pd.read_csv(path)

def read_results(triplet_path):
    return [ csv2pandas(path) for path in triplet_path]

def concat_df_list(df_list, output_path=None):
    df = pd.concat(df_list, ignore_index=True)
    if output_path:
        df.to_csv(output_path,index= False)
    else:
        return df

def main(params):

    ### create time dict
    times_dict = {
        'libraries_loading': params["import_time"]
    }

    ### make  output dir and fail dir if they don't already exist
    failed_dir = os.path.join(params["output_dirpath"], "failed_tasks")
    mkdir(failed_dir)

    ### check if some tasks did not launched  and move them to the failed_dir
    start = time.perf_counter()
    input_csv_list = [os.path.join(params['input_dirpath'], f"subset_{i}.csv") for i in range(params['array_task'])]
    # print(input_csv_list[:10])
    with Pool(processes=params['processes']) as pool:
        launch_status_list = pool.map(partial(assess_launch_status, output_dir=failed_dir), input_csv_list)
        pool.close()
        pool.join()
    times_dict['unlaunched_check'] = time.perf_counter()-start
    
    ### check if  some tasks did not finished and move them to the failed_dir
    start = time.perf_counter()
    launched_task_dirs_list = [
    f"{path[:-4]}" for i, path in enumerate(input_csv_list) if launch_status_list[i]]
    with Pool(processes=params['processes']) as pool:
        completion_status_list = pool.map(partial(assess_completion_status, output_dir=failed_dir), launched_task_dirs_list)
        pool.close()
        pool.join()
    times_dict['task_status_time'] = time.perf_counter()-start

    ### read results
    start = time.perf_counter()
    successful_task_dirs_list = [dir_path for i, dir_path in enumerate(launched_task_dirs_list) if completion_status_list[i]]
    csv_triplet_list = make_csv_list(successful_task_dirs_list)
    with Pool(processes=params['processes']) as pool:
        results_df_list = pool.map(read_results,csv_triplet_list)
        pool.close()
        pool.join()
    times_dict['read_results_time'] = time.perf_counter()-start

    ###  results processing
    start = time.perf_counter() 
    blocks_df_list,metrics_df_list,timelogs_df_list= zip(*results_df_list)
    blocks_df=concat_df_list(blocks_df_list)
    blocks_df= pd.DataFrame(list(set(blocks_df["blocks"])),columns=["blocks"]) 
    mkdir(params['output_libpath'].replace(params['output_libpath'].split("/")[-1],""))
    blocks_df.to_csv(params['output_libpath'],index=False)
    output_path= os.path.join(params['output_dirpath'], "metrics.csv")
    concat_df_list(metrics_df_list,output_path=output_path)
    output_path= os.path.join(params['output_dirpath'], "experiments_timelog.csv")
    concat_df_list(timelogs_df_list,output_path=output_path)
    times_dict['read_results_time'] = time.perf_counter()-start

    ### write final timelog 
    times_dict['total_time']= time.perf_counter() - params["global_start"]
    timelog_df= pd.DataFrame(times_dict,index=[0])
    output_path= os.path.join(params['output_dirpath'], "timelog.csv")
    timelog_df.to_csv(output_path,index=False)

if __name__ == '__main__':
    params = concatenate2lib_parser()
    params["import_time"] = duration
    params["global_start"]= global_start
    main(params)
