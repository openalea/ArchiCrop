from __future__ import annotations

import sys
import time as t
from math import ceil
from multiprocessing import Pool

sys.path.append('../data')
from archi_dict import archi_sorghum_angles as archi

from openalea.archicrop.simulation import define_archicrop_parameters, run_archicrop_and_light_parallel

if __name__ == '__main__':

    # Define the inputs for the simulation
    tec_file_xml='../data/02NT18SorgV2D1_tec.xml'
    plt_file_xml='../data/sorgho_imp_M_v10_plt.xml'
    stics_output_file='../data/mod_s02NT18SorgV2D1.sti'
    weather_file = '../data/ntarla_corr.2018'
    location = {  
    'longitude': 3.87,
    'latitude': 12.58,
    'altitude': 800,
    'timezone': 'Europe/Paris'}

    n_cpu = 3
    id_sim = list(range(n_cpu))
    light_inter = True

    daily_dynamics, param_sets, density, inter_row = define_archicrop_parameters(archi_params = archi, 
                                                                                tec_file = tec_file_xml, 
                                                                                plant_file = plt_file_xml, 
                                                                                dynamics_file = stics_output_file,
                                                                                pot_factor = 1.5)

    keys = list(param_sets.keys())
    chunk_size = ceil(len(keys) / n_cpu)

    params_sets_split = {
        i: {k: param_sets[k] for k in keys[i*chunk_size:(i+1)*chunk_size]}
        for i in id_sim
    }

    with Pool(n_cpu) as p:
        start_time = t.time()
        p.starmap_async(run_archicrop_and_light_parallel, 
                        [(id, params_sets_split[id], daily_dynamics, density, weather_file, location, inter_row, light_inter) 
                        for id in id_sim]).get()
        end_time = t.time()
        elapsed_time = (end_time - start_time)/3600
        print(f"Simulation time: {elapsed_time:.2f} hours for {len(param_sets)} simulations on {n_cpu} CPU")


