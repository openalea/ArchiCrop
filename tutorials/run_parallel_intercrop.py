from __future__ import annotations

import sys
import time as t
from math import ceil
from multiprocessing import Pool
from pathlib import Path

sys.path.append('../data')
from archi_dict import archi_sorghum, archi_maize as archi_1, archi_2

from openalea.archicrop.growth import demand_dist_bis
from openalea.archicrop.simulation import define_archicrop_parameters_IC, run_archicrop_and_light_parallel
from openalea.archicrop.stics_io import read_csv_file_intercrop

if __name__ == '__main__':

    path = Path("D:/ArchiCrop_for_STICS-IC_light/0-data/workspace_v11_gen/")

    plant_1 = "sorgum"
    plant_2 = "maize"

    # Define the inputs for the simulation
    tec_files_1=list(path.glob("sorghum_*_tec.xml"))[:1]
    plant_file_1='D:/ArchiCrop_for_STICS-IC_light/0-data/workspace_v11_gen/plant/sorgho_imp_M_v10_plt.xml'
    tec_files_2=list(path.glob("maize_*_tec.xml"))[:1]
    plant_file_2='D:/ArchiCrop_for_STICS-IC_light/0-data/workspace_v11_gen/plant/corn_LI_step2_MANT_plt.xml'

    file_csv = "D:/ArchiCrop_for_STICS-IC_light/2-outputs/simulations_stics_intercrops.csv"
    weather_file = '../data/ntarla_corr.2018'
    location = {  
    'longitude': 3.87,
    'latitude': 12.58,
    'altitude': 800,
    'timezone': 'Europe/Paris'}

    d_outputs = read_csv_file_intercrop(file_csv)

    print("Nb CPU : ")
    n_cpu = int(input())
    id_sim = list(range(n_cpu))

    pot_factor_lai = 3
    pot_factor_height = 5
    distribution_function = demand_dist_bis
    light_inter = False
    direct = True

    param_sets = {}

    for i, t1, t2 in enumerate(zip(tec_files_1, tec_files_2)):

        usm = f"usm_{i+1}"

        for algo in ["Beer", "2.5D"]:

            param_sets_1 = define_archicrop_parameters_IC(archi_params = archi_1, 
                                                    tec_file = t1, 
                                                    plant_file = plant_file_1, 
                                                    d_outputs=d_outputs, usm=usm, algo=algo, plant=plant_1,
                                                    pot_factor_lai = pot_factor_lai,
                                                    pot_factor_height = pot_factor_height)
            
            param_sets_2 = define_archicrop_parameters_IC(archi_params = archi_2, 
                                                    tec_file = t2, 
                                                    plant_file = plant_file_1, 
                                                    d_outputs=d_outputs, usm=usm, algo=algo, plant=plant_2,
                                                    pot_factor_lai = pot_factor_lai,
                                                    pot_factor_height = pot_factor_height)
            
            param_sets[usm][algo][plant_1] = param_sets_1
            param_sets[usm][algo][plant_2] = param_sets_2


    # Config spatiale
    # density_1 = 5
    # density_2 = 5
    # inter_row_1 = 0.70
    # inter_row_2 = 0.70
    # width = 8
    # length = 2
    # nb_rows_1 = 1
    # nb_rows_2 = 1        

    keys = list(param_sets.keys())
    chunk_size = ceil(len(keys)*2 / n_cpu)
    params_sets_split = {
        i: {k: param_sets[k] for k in keys[i*chunk_size:(i+1)*chunk_size]}
        for i in id_sim
    }
    print(f"Run {len(param_sets)} simulations on {n_cpu} CPU")

    with Pool(n_cpu) as p:
        start_time = t.time()
        p.starmap_async(run_archicrop_and_light_parallel, 
                        [(id, params_sets_split[id], tec_file, plant_file, stics_output_file, weather_file, location, distribution_function, light_inter, direct) 
                        for id in id_sim]).get()
        end_time = t.time()
        elapsed_time = (end_time - start_time)/3600
        print(f"Simulation time: {elapsed_time:.2f} hours for {len(param_sets)*2} simulations on {n_cpu} CPU")


