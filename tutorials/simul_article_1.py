import sys
import time as t
from datetime import date

sys.path.append("../data") 
from archi_dict import archi_sorghum_angles as archi

from openalea.archicrop.simulation import (
    run_simulations,
    write_netcdf,
)

tec_file_xml='../data/02NT18SorgV2D1_tec.xml'
plt_file_xml='../data/sorgho_imp_M_v10_plt.xml'
stics_output_file='../data/mod_s02NT18SorgV2D1.sti'
weather_file = '../data/ntarla_corr.2018'
location = {  
'longitude': 3.87,
'latitude': 12.58,
'altitude': 800,
'timezone': 'Europe/Paris'}

seed = 18

daily_dynamics, params_sets, pot_la, pot_h, realized_la, realized_h, nrj_per_plant, mtgs, density = run_simulations(
    archi_params=archi, 
    tec_file=tec_file_xml, 
    plant_file=plt_file_xml, 
    dynamics_file=stics_output_file, 
    weather_file=weather_file,
    location=location,
    n_samples=4,
    pot_factor=1.4,
    latin_hypercube=False,
    light_inter=True,
    direct=False,
    seed=seed)

write_netcdf("results_light_inter", daily_dynamics, params_sets, pot_la, pot_h, realized_la, realized_h, nrj_per_plant, density, seed)