from __future__ import annotations

import sys
import time as t
from math import ceil
from multiprocessing import Pool

sys.path.append('../data')
from archi_dict import archi_sorghum_angles as archi

from openalea.archicrop.growth import dev_dist, demand_dist, demand_dist_bis
from openalea.archicrop.simulation import define_archicrop_parameters, run_archicrop_and_light_parallel

if __name__ == '__main__':

    # Define the inputs for the simulation
    tec_file='../data/02NT18SorgV2D1_tec.xml'
    plant_file='../data/sorgho_imp_M_v10_plt.xml'
    stics_output_file='../data/mod_s02NT18SorgV2D1.sti'
    stics_output_file_no_stress='../data/mod_s02NT18SorgV2D1_no_stress.sti'
    weather_file = '../data/ntarla_corr.2018'
    location = {  
    'longitude': 3.87,
    'latitude': 12.58,
    'altitude': 800,
    'timezone': 'Europe/Paris'}

    print("Nb CPU : ")
    n_cpu = int(input())
    id_sim = list(range(n_cpu))
    pot_factor_lai = 1.4
    pot_factor_height = 3
    distribution_function = demand_dist_bis
    end = -1
    light_inter = True
    direct = False

    param_sets = define_archicrop_parameters(archi_params = archi, 
                                            tec_file = tec_file, 
                                            plant_file = plant_file, 
                                            pot_dynamics_file = stics_output_file_no_stress,
                                            pot_factor_lai = pot_factor_lai,
                                            pot_factor_height = pot_factor_height,
                                            end=end)

    keys = list(param_sets.keys())
    chunk_size = ceil(len(keys) / n_cpu)
    params_sets_split = {
        i: {k: param_sets[k] for k in keys[i*chunk_size:(i+1)*chunk_size]}
        for i in id_sim
    }
    print(f"Run {len(param_sets)} simulations on {n_cpu} CPU")

    with Pool(n_cpu) as p:
        start_time = t.time()
        p.starmap_async(run_archicrop_and_light_parallel, 
                        [(id, params_sets_split[id], tec_file, plant_file, stics_output_file, weather_file, location, distribution_function, end, light_inter, direct) 
                        for id in id_sim]).get()
        end_time = t.time()
        elapsed_time = (end_time - start_time)/3600
        print(f"Simulation time: {elapsed_time:.2f} hours for {len(param_sets)} simulations on {n_cpu} CPU")


