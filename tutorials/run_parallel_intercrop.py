from __future__ import annotations

import math
import sys
import time as t
from math import ceil
from multiprocessing import Pool
from pathlib import Path

sys.path.append('../data')
from archi_dict import archi_sorghum, archi_maize as archi_1, archi_2

from openalea.archicrop.growth import demand_dist_bis
from openalea.archicrop.simulation import define_archicrop_parameters_IC, run_archicrop_and_light_parallel_IC
from openalea.archicrop.stics_io import read_csv_file_IC, read_doe_intercrop

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

    d_outputs = read_csv_file_IC(file_csv)

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

            param_sets_1, density_1 = define_archicrop_parameters_IC(archi_params = archi_1, 
                                                    tec_file = t1, 
                                                    plant_file = plant_file_1, 
                                                    d_outputs=d_outputs, usm=usm, algo=algo, plant=plant_1,
                                                    pot_factor_lai = pot_factor_lai,
                                                    pot_factor_height = pot_factor_height)
            
            param_sets_2, density_2 = define_archicrop_parameters_IC(archi_params = archi_2, 
                                                    tec_file = t2, 
                                                    plant_file = plant_file_1, 
                                                    d_outputs=d_outputs, usm=usm, algo=algo, plant=plant_2,
                                                    pot_factor_lai = pot_factor_lai,
                                                    pot_factor_height = pot_factor_height)
            
            param_sets[usm][algo][plant_1] = (param_sets_1, density_1)
            param_sets[usm][algo][plant_2] = (param_sets_2, density_2)
  

    row_orientation_values = {
    "N-S" : 0,
    "E-W" : math.pi / 2
    }

    interrow_distance_per_species = {
    "sorghum" : {
        "high" : 0.8,
        "middle" : 0.4,
        "low" : 0.2
    },
    "maize" : {
        "high" : 0.8,
        "middle" : 0.4,
        "low" : 0.2
    }
    }

    n_rows_per_species = {
    "sorghum" : {
        "one" : 1, # For non-strip, this is just one row
        "high" : 6,
        "middle" : 4,
        "low" : 2
    },
    "maize" : {
        "one" : 1,
        "high" : 6,
        "middle" : 4,
        "low" : 2
    }
    }

    intrarow_distance_per_species = {
    "sorghum" : {
        "high" : 0.8,
        "middle" : 0.4, # ~6 plants per m2 with 0.4m interrow distance, gives 0.41m intrarow distance
        "low" : 0.2
    },
    "maize" : {
        "high" : 0.8,
        "middle" : 0.4,
        "low" : 0.2
    }
    }

    doe_file = "D:/ArchiCrop_for_STICS-IC_light/2-outputs/doe.csv"
    doe = read_doe_intercrop(doe_file)

    doe_adapt = {}

    for usm,spat_conf in doe.items():
        spat_conf["row_orientation"] = row_orientation_values[spat_conf["row_orientation"]]
        spat_conf["interrow_distance_principal"] = interrow_distance_per_species[spat_conf["species_principal"]][spat_conf["interrow_distance_principal"]]
        spat_conf["interrow_distance_secondary"] = interrow_distance_per_species[spat_conf["species_secondary"]][spat_conf["interrow_distance_secondary"]]
        spat_conf["n_rows_principal"] = n_rows_per_species[spat_conf["species_principal"]][spat_conf["n_rows_principal"]]
        spat_conf["n_rows_secondary"] = n_rows_per_species[spat_conf["species_secondary"]][spat_conf["n_rows_secondary"]]
        spat_conf["intrarow_distance"] = intrarow_distance_per_species[spat_conf["species_principal"]][spat_conf["intrarow_distance"]]

        for algo in param_sets[usm]:

            doe_adapt[usm][algo] = {
                "density_1" : param_sets[usm][algo][plant_1][1],
                "density_2" : param_sets[usm][algo][plant_2][1],
                "inter_row_1" : spat_conf["interrow_distance_principal"],
                "inter_row_2" : spat_conf["interrow_distance_secondary"],
                "width" : 2 * spat_conf["interrow_distance_principal"] if spat_conf["row_orientation"] == "intercrop mixed" else (spat_conf["n_rows_principal"]+1) * spat_conf["interrow_distance_principal"] + (spat_conf["n_rows_secondary"]-1) * spat_conf["interrow_distance_secondary"],
                "length" : 2 * spat_conf["intrarow_distance"] if spat_conf["row_orientation"] == "intercrop mixed" else spat_conf["intrarow_distance"],
                "nb_rows_1" : spat_conf["n_rows_principal"],
                "nb_rows_2" : spat_conf["n_rows_secondary"]
            }


    # First trial on usms 1 to 6  

    # keys = list(param_sets.keys())
    # chunk_size = ceil(len(keys)*2 / n_cpu)
    # params_sets_split = {
    #     i: {k: param_sets[k] for k in keys[i*chunk_size:(i+1)*chunk_size]}
    #     for i in id_sim
    # }
    # print(f"Run {len(param_sets)} simulations on {n_cpu} CPU")

    with Pool(n_cpu) as p:
        start_time = t.time()
        p.starmap_async(run_archicrop_and_light_parallel_IC, 
                        [(usm, param_sets[usm], d_outputs[usm], doe_adapt[usm], weather_file, location, distribution_function, light_inter, direct) 
                        for usm in id_sim]).get()
        end_time = t.time()
        elapsed_time = (end_time - start_time)/3600
        print(f"Simulation time: {elapsed_time:.2f} hours for {len(param_sets)*2} simulations on {n_cpu} CPU")


