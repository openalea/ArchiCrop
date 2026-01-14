from __future__ import annotations

import sys
import time as t

sys.path.append('../data')
from archi_dict import archi_sorghum_angles as archi

from openalea.archicrop.simulation import run_simulations, write_netcdf

seed = 18

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

start_time = t.time()

# Run the simulation
daily_dynamics, params_sets, pot_la, pot_h, realized_la, realized_h, nrj_per_plant, mtgs, density = run_simulations(
    archi_params=archi, 
    tec_file=tec_file_xml, 
    plant_file=plt_file_xml, 
    dynamics_file=stics_output_file, 
    weather_file=weather_file,
    location=location,
    n_samples=3,
    pot_factor=1.5,
    latin_hypercube=False,
    light_inter=True,
    direct=False,
    seed=seed)

end_time = t.time()

elapsed_time = (end_time - start_time)/3600
print(f"Simulation time: {elapsed_time:.2f} hours for {len(realized_la)} simulations")  # noqa: T201

write_netcdf("results_light_inter", daily_dynamics, params_sets, pot_la, pot_h, realized_la, realized_h, nrj_per_plant, density, seed)

print(f"Simulations saved")  # noqa: T201
