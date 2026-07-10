from __future__ import annotations

from random import *  # noqa: F403
import numpy as np
import matplotlib.pyplot as plt
from openalea.caribu.data_samples import data_path

from openalea.archicrop.archicrop import ArchiCrop
from openalea.archicrop.display import build_scene
from openalea.archicrop.light_it import compute_light_inter
from openalea.archicrop.simulation import run_simulations
from openalea.archicrop.stand import agronomic_plot
from openalea.archicrop.stics_io import get_stics_data
from openalea.plantgl.all import Color3, Material

seed(18)  # noqa: F405

# Retrieve STICS management and senescence parameters
tec_file='../data/02NT18SorgV2D1_tec.xml'
plant_file='../data/plant/sorgho_imp_M_v10_plt.xml'
stics_output_file='../data/mod_s02NT18SorgV2D1.sti'
stics_output_file_no_stress='../data/mod_s02NT18SorgV2D1_no_stress.sti'
weather_file = '../data/ntarla_corr.2018'
location = {  
'longitude': 3.87,
'latitude': 12.58,
'altitude': 800,
'timezone': 'Europe/Paris'}

density, stics_output_data, lifespan, lifespan_early, inter_row = get_stics_data(
    file_tec_xml=tec_file,  # Path to the STICS management XML file
    file_plt_xml=plant_file,  # Path to the STICS plant XML file
    stics_output_file=stics_output_file  # Path to the STICS output file
)

# Retrieve STICS growth and senescence dynamics
time = [value["Thermal time"] for value in stics_output_data.values() if value is not None]
LA_stics = [value["Plant leaf area"] for value in stics_output_data.values() if value is not None]
sen_LA_stics = [value["Plant senescent leaf area"] for value in stics_output_data.values() if value is not None]
height_stics = [value["Plant height"] for value in stics_output_data.values() if value is not None]
par_stics = [value["Absorbed PAR"] for value in stics_output_data.values() if value is not None]


# Set ArchiCrop parameters
archi = {
        'height': 2*max(height_stics),
        'leaf_area': 1.5*max(LA_stics),
        'leaf_duration': 1.6,
        'nb_phy': 16,
        'nb_short_phy': 4,
        'wl': 0.12,
        'diam_base': 4.0,
        'diam_top': 1.5,
        'insertion_angle': 60,
        'scurv': 0.7,
        'curvature': 90,
        'klig': 0.6,
        'swmax': 0.55,
        'f1': 0.64,
        'f2': 0.92,
        'stem_q': 1,
        'rmax': [0.5,0.95],
        'skew': [0.0005,1],
        'phyllotactic_angle': 180,
        'phyllotactic_deviation': 20,
        'phyllochron': [30,80],
        'nb_tillers': 0,
        'tiller_delay': 1,
        'tiller_angle': 20,
        'reduction_factor': 1,
}


# Run the simulation
daily_dynamics, params_sets, pot_la, pot_h, realized_la, realized_h, nrj_per_plant, mtgs, density, inter_row = run_simulations(
    archi_params=archi, 
    tec_file=tec_file, 
    plant_file=plant_file, 
    dynamics_file=stics_output_file, 
    pot_dynamics_file=stics_output_file_no_stress,
    weather_file=weather_file,
    location=location,
    n_samples=1,
    pot_factor_lai=1.5,
    pot_factor_height=2,
    latin_hypercube=False,
    light_inter=True,
    direct=False,
    zenith=False,
    seed=seed)

# # Instanciate ArchiCrop object
# plant = ArchiCrop(daily_dynamics=stics_output_data, **archi)
# # Generate a potential plant
# plant.generate_potential_plant()
# # Simulate growth and senescence of this plant according to the STICS dynamics
# growing_plant = plant.grow_plant()

# # Sky
# zenith = str(data_path('zenith.light'))

# nice_green=Color3((50,100,0))

# nplants, positions, domain, domain_area, unit = agronomic_plot(length=1, width=1, density=density, inter_row=inter_row, noise=0.1)

# inter_plant = 1/density/inter_row
# x_pattern = inter_plant/2
# y_pattern = inter_row/2
# pattern = (-x_pattern, -y_pattern, x_pattern, y_pattern)
# position = (0,0,0)

# scenes = {}
# for k,v in mtgs[0].items():
#     scene, nump = build_scene(v, position, leaf_material=Material(nice_green), stem_material=Material(nice_green), senescence=True)
#     scenes[k] = scene

# par_caribu = []
# for scene in scenes.values():
#     par_caribu.append(compute_light_inter(scene, zenith, pattern)[0])

nrj_crop = {id:[np.nansum([leaf for leaf in t])*density for t in nrj] for id,nrj in nrj_per_plant.items()}

# plt.plot(time, nrj_crop, alpha=0.5, linestyle='--')  # Plot each curve (optional for visualization)
# plt.plot(time, par_stics, color="black", label="STICS")

# # Labels and legend
# plt.xlabel("Thermal time")
# plt.ylabel("% of absorbed PAR")
# plt.title("Absorbed PAR: 3D canopy vs. STICS")
# plt.legend()
# plt.show()

# %gui qt
# %run test_use_archicrop.py