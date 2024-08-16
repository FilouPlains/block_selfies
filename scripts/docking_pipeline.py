""" This scripts handle the docking experiment and data collection of results
"""
import time
import os
import subprocess
import argparse
import csv
import statistics


def main(params: dict) -> None:
    """ This function handle the docking pipeline
    Args:
        params (dict): this dictionnary contains all the argument for the function
        for details please refers to the argparser documentation with --help flag 
        '
    """

    # declare local variables
    first_time = time.perf_counter()
    metrics_dict = {
        'filename': params["ligand_path"].split("/")[-1],
        'docking_duration': 0.,
        'data_handling_duration': 0.,
        'average_score': None,
        'std_score': None,
    }
    scores_list = []

    # create sub_directory
    dir_path = params["ligand_path"].replace(".pdbqt", "/")
    os.makedirs(dir_path, exist_ok=True)

    # docking prep
    output_path = os.path.join(
        dir_path, metrics_dict["filename"].replace(".pdbqt", "_out.pdbqt"))
    smina_logpath = os.path.join(
        dir_path, metrics_dict["filename"].replace(".pdbqt", "_log.txt"))

    # docking (make this part as function, when implementing RNA docking)

    start = time.perf_counter()
    # prepare and execute subprocess command
    cmd = f'{params["smina_path"]} --receptor {params["receptor_path"]} '
    cmd += f'--ligand { params["ligand_path"]}'
    cmd += f' --out {output_path} --autobox_ligand {params["ref_ligand_path"]}'
    cmd += f' --exhaustiveness {params["exhaustiveness"]} '
    cmd += f'--scoring {params["scoring_function"]} '
    cmd += f'--log {smina_logpath} --cpu {params["nb_core"]}'  # --quiet
    if params["score_only"]:
        cmd+=" --score_only"
    subprocess.run(cmd.split(),
                   timeout=params["docking_walltime"],
                   check=True)
    metrics_dict['docking_duration'] = time.perf_counter() - start

    # retrieve scores
    start = time.perf_counter()
    with open(output_path, 'r', encoding="utf-8") as file:
        for line in file:
            if line.startswith('REMARK minimizedAffinity'):
                scores_list.append(float(line.split()[2]))

    # compute statistics on score
    if not params["score_only"]:
        metrics_dict["average_score"] = statistics.fmean(scores_list)
        metrics_dict["std_score"] = statistics.stdev(
            scores_list)  # N-1 degree of freedom

    # write results as csv file
    output_path = output_path.replace("_out.pdbqt", "_results.csv")
    with open(output_path, "w", encoding="utf-8") as file:
        csv_writer = csv.writer(file, delimiter=',',
                                quotechar='"', quoting=csv.QUOTE_MINIMAL)

        # write column name
        csv_writer.writerow(['filename',
                             'pose_number',
                             'score_function',
                             'score'])
        # write score to file
        for i, score in enumerate(scores_list):
            csv_writer.writerow([metrics_dict["filename"],
                                 i,
                                 params["scoring_function"],
                                 score])
    metrics_dict["data_handling_duration"] = time.perf_counter() - start

    # write metrics
    output_path = output_path.replace("_results.csv", "_metrics.csv")
    metrics_dict["total_time"] = time.perf_counter()-first_time
    with open(output_path, "w", encoding="utf-8") as file:
        csv_writer = csv.DictWriter(file, fieldnames=metrics_dict.keys())
        csv_writer.writeheader()
        csv_writer.writerow(metrics_dict)


if __name__ == '__main__':

    # create argument parser
    parser = argparse.ArgumentParser()

    # files & dir
    parser.add_argument('--ligand_path',
                        help="path to the ligand to be docked, should be a .pdbqt",
                        default="data/celecoxib_chain_A.pdbqt",
                        type=str)
    parser.add_argument('--receptor_path', help="path to the receptor used for docking",
                        default="data/COX2_chain_A.pdbqt", type=str)
    parser.add_argument('--ref_ligand_path', help="path to the ligand to build the box on," +
                        " a .pdbqt is expected but mol2 should work too",
                        default="data/celecoxib_chain_A.pdbqt", type=str)

    # SMINA parameters
    parser.add_argument('--smina_path', help="path to the docking software smina.static",
                        default="./smina.static", type=str)
    parser.add_argument('--exhaustiveness',
                        help=" exhausstiveness of the conformational space search " +
                        " for ligand, smina default is 8",
                        type=int,
                        default=8)
    parser.add_argument('--scoring_function', help="scoring functions to be used for docking",
                        type=str,
                        default="vinardo")
    parser.add_argument('--score_only', help="use score only mode",
                        default=False,
                        action="store_true")
    parser.add_argument('--nb_core', help="number of core used for the docking,"
                        + " if using GNU, should always be 1",
                        type=int,
                        default=1)
    parser.add_argument('--docking_walltime',
                        help="maximum duration allowed for docking subprocess, in secondes",
                        type=int, default=600)

    # parse arguments
    args, _ = parser.parse_known_args()
    args_dict = vars(args)

    # execute docking
    main(args_dict)
