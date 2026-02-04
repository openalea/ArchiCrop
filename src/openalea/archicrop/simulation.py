from __future__ import annotations

import itertools
import math
import os
from datetime import date

import numpy as np
import pandas as pd
import xarray as xr
from scipy.stats import qmc

from .archicrop import ArchiCrop
from .growth import dev_dist, demand_dist, demand_dist_bis
from .light_it import light_interception
from .stics_io import get_stics_data
from .viability import compute_viable_params



def constrain_archi_params(archi_params: dict, daily_dynamics: dict, 
                           lifespan: float = None, lifespan_early: float = None, 
                           pot_factor_lai: float = 1.5, pot_factor_height: float = 3) -> dict:
    """Complete ArchiCrop parameters with values from constraint."""

    if daily_dynamics is not None:
        leaf_area_plant = [value["Plant leaf area"] for value in daily_dynamics.values()]
        height_canopy = [value["Plant height"] for value in daily_dynamics.values()]
        archi_params["height"] = pot_factor_height*max(height_canopy) 
        archi_params["leaf_area"] = pot_factor_lai*max(leaf_area_plant) 

    if lifespan is not None:
        archi_params["leaf_lifespan"] = lifespan
        archi_params["leaf_lifespan_early"] = lifespan_early

    return archi_params



def param_sampling(archi_params: dict, n_samples: int, seed: int = 42, latin_hypercube: bool = False) -> dict:
    """Generate samples for ArchiCrop parameters."""

    fixed_params = {}
    sampled_params = []

    l_bounds = []
    u_bounds = []

    values_lists = []

    for key, value in archi_params.items():
        if isinstance(value, (int, float)):  # Fixed parameter
            fixed_params[key] = value
        # Parameter distribution in Latin Hypercube
        elif isinstance(value, list) and key not in ["rmax", "skew", "phyllochron"]:  # Range to sample

            if latin_hypercube:
                l_bounds.append(min(value))
                u_bounds.append(max(value))
            else:
                # Discretize each parameter
                # value_list = np.linspace(min(value), max(value), n_samples)
                # values_lists.append(value_list)
                values_lists.append(value)

            sampled_params.append(key)

        elif key in ["rmax", "skew", "phyllochron"]:
            fixed_params[key] = value


    if latin_hypercube:
        # Create a Latin Hypercube sampler
        sampler = qmc.LatinHypercube(d=len(l_bounds), seed=seed)  # d = number of parameters
        
        # Generate normalized samples (0 to 1)
        lhs_samples = sampler.random(n=n_samples)
        
        # Scale samples to parameter bounds
        samples = qmc.scale(lhs_samples, l_bounds=l_bounds, u_bounds=u_bounds)
    else:

        samples = list(itertools.product(*values_lists))

    # Create parameter sets
    param_sets = {}
    for i,sample in enumerate(samples):
        sampled_dict = {
            key: int(value) if key in {"nb_phy", "nb_tillers", "nb_short_phy"} else round(value, 4)
            for key, value in zip(sampled_params, sample)
        }
        param_sets[i] = {**fixed_params, **sampled_dict}
    
    return param_sets


def cumsum_organs(g, thermal_time: list, starts_ends: list, property: str) -> list:
    '''Compute the sum of all organs dimension, i.e. the plant-scale dimension at each thermal time step'''
    cumsum = []
    for t in thermal_time:
        sum_organs = 0
        for (vid, s,e) in starts_ends:
            if s <= t < e:
                sum_organs += (t-s)/(e-s) * g.properties()[property][vid]
            elif t >= e:
                sum_organs += g.properties()[property][vid]
        cumsum.append(sum_organs)
    return cumsum



def compute_pot_growth(param_sets: dict, daily_dynamics: dict):
    '''
    Compute theretical growth of ArchiCrop plants.
    Parameters:
        param_sets (dict): List of parameter sets to filter.
        daily_dynamics (dict): Dictionary containing daily dynamics data.
    Returns:
        pot_la (list): List of potential leaf area for each parameter set.
        pot_h (list): List of potential height for each parameter set.
    '''
    thermal_time = [value["Thermal time"] for value in daily_dynamics.values()]

    pot_la = {}
    pot_h = {}
    mtgs = {}
    
    for id, params in param_sets.items():
        plant = ArchiCrop(daily_dynamics=daily_dynamics, **params)
        plant.generate_potential_plant()
        g = plant.g
        mtgs[id] = g

        stem_starts_ends = []
        leaf_starts_ends = []
        for vid in g.properties()["start_tt"]:
            start =  g.node(vid).start_tt
            end = g.node(vid).end_tt
            if g.node(vid).label.startswith("Stem"):
                stem_starts_ends.append((vid, start, end))
            elif g.node(vid).label.startswith("Leaf"):
                leaf_starts_ends.append((vid, start, end))
        # Order by start time
        stem_starts_ends.sort(key=lambda x: x[1])  
        leaf_starts_ends.sort(key=lambda x: x[1])

        pot_la[id] = cumsum_organs(g, thermal_time, leaf_starts_ends, 'leaf_area')
        pot_h[id] = cumsum_organs(g, thermal_time, stem_starts_ends, 'mature_length')

    return pot_la, pot_h, mtgs
    

def simulate_plant_growth(param_sets: dict, daily_dynamics: dict, distribution_function):
    """ Simulate ArchiCrop plant growth constrained by daily dynamics.
    Parameters:
        param_sets (dict): List of parameter sets.
        daily_dynamics (dict): Dictionary containing daily dynamics data.
    Returns:
        realized_la (dict): List of realized leaf area for each parameter set.
        realized_h (dict): List of realized height for each parameter set.
        mtgs (dict): List of MTG objects for each parameter set.
    """
    
    realized_la = {}
    realized_h = {}
    mtgs = {}

    for id, params in param_sets.items():
        plant = ArchiCrop(daily_dynamics=daily_dynamics, **params)
        plant.generate_potential_plant()
        growing_plant = plant.grow_plant(distribution_function=distribution_function)
        growing_plant_mtg = list(growing_plant.values())
        mtgs[id] = growing_plant_mtg

        realized_la[id] = [sum(la - sen 
                            for la, sen in zip(gp.properties()["visible_leaf_area"].values(), gp.properties()["senescent_area"].values())) 
                            for gp in growing_plant_mtg]
        realized_h[id] = [sum([vl 
                                for vid, vl in gp.properties()["visible_length"].items() 
                                if gp.node(vid).label.startswith("Stem")]) 
                                for gp in growing_plant_mtg]
        
    return realized_la, realized_h, mtgs



def define_archicrop_parameters(archi_params: dict, 
             tec_file: str, plant_file: str, pot_dynamics_file: str,
             n_samples: int = 10, seed: int = 42, latin_hypercube: bool = False,
             pot_factor_lai: float = 1.5, pot_factor_height: float = 3, end: int = -1):

    # Retrieve STICS management and senescence parameters
    _, pot_daily_dynamics, lifespan, lifespan_early, _ = get_stics_data(
        file_tec_xml=tec_file,  # Path to the STICS management XML file
        file_plt_xml=plant_file,  # Path to the STICS plant XML file
        stics_output_file=pot_dynamics_file,  # Path to the STICS output file
        end=end
    )

    # Complete ArchiCrop parameters with values from constraint
    archi_params = constrain_archi_params(archi_params=archi_params, daily_dynamics=pot_daily_dynamics, 
                                          lifespan=lifespan, lifespan_early=lifespan_early, 
                                          pot_factor_lai=pot_factor_lai, pot_factor_height=pot_factor_height)
    
    # Generate samples for ArchiCrop parameters 
    param_sets = param_sampling(archi_params=archi_params, n_samples=n_samples, seed=seed, latin_hypercube=latin_hypercube)

    # Compute viable parameters wrt the dynamic constraint of LAI and height
    param_sets = compute_viable_params(
        params_sets=param_sets, 
        daily_dynamics=pot_daily_dynamics, 
        pot_factor_lai=pot_factor_lai, pot_factor_height=pot_factor_height)

    return param_sets


def run_archicrop_and_light(param_sets: dict, 
             tec_file: str, plant_file: str, dynamics_file: str,
             weather_file: str, location: dict,
             distribution_function,
             end: int = -1,
             light_inter: bool = True, zenith: bool = False, direct : bool = False, save_scenes: bool = False):

    # Retrieve STICS management and senescence parameters
    density, daily_dynamics, _, _, inter_row = get_stics_data(
        file_tec_xml=tec_file,  # Path to the STICS management XML file
        file_plt_xml=plant_file,  # Path to the STICS plant XML file
        stics_output_file=dynamics_file,  # Path to the STICS output file
        end=end
    )

    # Simulate plant growth with fitting parameters
    pot_la, pot_h, _ = compute_pot_growth(param_sets=param_sets, daily_dynamics=daily_dynamics)
    realized_la, realized_h, mtgs = simulate_plant_growth(param_sets=param_sets, daily_dynamics=daily_dynamics, distribution_function=distribution_function)

    # Compute light interception on 3D scenes through time
    if light_inter:
        nrj_per_plant = light_interception(
            weather_file=weather_file, 
            daily_dynamics=daily_dynamics, 
            density=density, 
            location=location, 
            mtgs=mtgs, 
            zenith=zenith, 
            direct=direct,
            save_scenes=save_scenes, 
            inter_row=inter_row
        )
    else:
        nrj_per_plant = {k : [None] * len(realized_la[k]) for k in realized_la}

    return daily_dynamics, density, inter_row, pot_la, pot_h, realized_la, realized_h, nrj_per_plant, mtgs

def run_archicrop_and_light_parallel(id_sim, param_sets: dict, tec_file: str, plant_file: str, dynamics_file: str,
             weather_file: str, location: dict, distribution_function, end: int = -1, light_inter: bool = True, direct: bool = False):
    
    daily_dynamics, density, inter_row, pot_la, pot_h, realized_la, realized_h, nrj_per_plant, _ = run_archicrop_and_light(param_sets=param_sets, 
                                                                                       tec_file=tec_file, plant_file=plant_file, 
                                                                                       dynamics_file=dynamics_file, weather_file=weather_file, 
                                                                                       location=location, 
                                                                                       distribution_function=distribution_function,
                                                                                       end=end, light_inter=light_inter, direct=direct)
    
    # first_key = list(param_sets.keys())[0]
    # filename = f"results_light_inter_{param_sets[first_key]['nb_phy']}"
    filename = "results_light_inter"
    write_netcdf(filename, daily_dynamics, param_sets, 
                 pot_la, pot_h, realized_la, realized_h, nrj_per_plant, 
                 density, inter_row, id_sim, 
                 dir = "../simulations_ArchiCrop/")
    

def run_simulations(archi_params: dict, 
             tec_file: str, plant_file: str, pot_dynamics_file: str, dynamics_file: str, weather_file: str, location: dict,
             distribution_function,
             n_samples: int = 10, seed: int = 42, latin_hypercube: bool = False, 
             pot_factor_lai: float = 1.5, pot_factor_height: float = 3,
             light_inter: bool = True, zenith: bool = False, direct : bool = False, save_scenes: bool = False):

    param_sets = define_archicrop_parameters(archi_params=archi_params, tec_file=tec_file, plant_file=plant_file, 
                                            pot_dynamics_file=pot_dynamics_file, 
                                            n_samples=n_samples, seed=seed, latin_hypercube=latin_hypercube, 
                                            pot_factor_lai=pot_factor_lai, pot_factor_height=pot_factor_height)

    
    # print(len(param_sets), "simulations")

    daily_dynamics, density, inter_row, pot_la, pot_h, realized_la, realized_h, nrj_per_plant, mtgs = run_archicrop_and_light(param_sets=param_sets, 
                                                                                       tec_file=tec_file, plant_file=plant_file, 
                                                                                       dynamics_file=dynamics_file, weather_file=weather_file, 
                                                                                       location=location, 
                                                                                       distribution_function=distribution_function,
                                                                                       end=-1, light_inter=light_inter, zenith=zenith, 
                                                                                       direct=direct)

    return daily_dynamics, param_sets, pot_la, pot_h, realized_la, realized_h, nrj_per_plant, mtgs, density, inter_row



def compute_extinction_coef(nrj_crop: dict, par_incident: list, leaf_area_plant:list, density: float, conv_unit: float = 100):
    '''Compute daily extinction coefficient of a crop.'''

    k_dict = {}
    for id,nrj_time_series in nrj_crop.items():
        k_per_sim = []
        for i,(nrj,par) in enumerate(zip(nrj_time_series, par_incident)):
            fapar = nrj/par
            if fapar <= 1:
                lai = leaf_area_plant[i]*density/(conv_unit**2)
                k = -1/lai * math.log(1-fapar)
            else:
                k = 1
            k_per_sim.append(k)
        k_dict[id] = k_per_sim
    return k_dict



def write_netcdf(filename: str, daily_dynamics: dict, params_sets: dict, 
                 pot_la: dict, pot_h: dict, realized_la: dict, realized_h: dict, nrj_per_leaf: dict, 
                 density: float, inter_row: float, seed: int, 
                 dir: str = "../simulations_ArchiCrop/"):
    
    # Prepare the data for xarray Dataset
    daily_dyn = {k: [v[k] for v in daily_dynamics.values()]
                for k in daily_dynamics[1]}
    dates = pd.to_datetime(daily_dyn["Date"])

    # ds_archi = (
    #     pd.DataFrame.from_dict(params_sets, orient="index")
    #     .to_xarray()
    #     .rename({"index": "id"})
    # )

    # ds_archi = (
    #     pd.DataFrame.from_dict(params_sets, orient="index")
    #     .rename_axis("id")      # rename pandas index
    #     .to_xarray()
    # )

    df_archi = pd.DataFrame.from_dict(params_sets, orient='index')
    # print(df_archi)
    ds_archi = df_archi.to_xarray().rename({'index':'id'})

    # mtgs_string = {k: [write_mtg(g) if g is not None else None for g in mtg] for k, mtg in mtgs.items()}

    # Per-id time series
    df_realized_la = pd.DataFrame.from_dict(realized_la, orient="index", columns=dates)
    df_realized_h = pd.DataFrame.from_dict(realized_h, orient="index", columns=dates)
    df_pot_la = pd.DataFrame.from_dict(pot_la, orient="index", columns=dates)
    df_pot_h = pd.DataFrame.from_dict(pot_h, orient="index", columns=dates)
    # df_nrj_per_plant = pd.DataFrame.from_dict(nrj_per_plant, orient="index", columns=dates)

    # print(df_realized_la)
    # print(df_pot_la)
    # print(df_nrj_per_plant)

    # # 3D array for nrj per leaf 
    ids = list(nrj_per_leaf.keys())
    n_id = len(ids)
    n_time = len(dates)
    n_leaf = max([len(nrj_per_leaf[id][-1]) for id in ids])  # length of inner list for a plant with all leaves


    nrj_leaf = np.zeros((n_id, n_time, n_leaf), dtype=float)

    for i, plant_id in enumerate(ids):
        for j, t in enumerate(nrj_per_leaf[plant_id]):
            nrj_leaf[i, j, : len(t)] = t

    ds = xr.Dataset(
        data_vars=dict(
            thermal_time=("time", daily_dyn["Thermal time"]),
            lai_stics=("time", daily_dyn["Plant leaf area"]),
            sen_lai_stics=("time", daily_dyn["Plant senescent leaf area"]),
            height_stics=("time", daily_dyn["Plant height"]),
            inc_par=("time", daily_dyn["Incident PAR"]),
            abs_par_stics=("time", daily_dyn["Absorbed PAR"]),
            density=density,
            inter_row=inter_row,
            realized_la=(("id", "time"), df_realized_la),
            realized_h=(("id", "time"), df_realized_h),
            pot_la=(("id", "time"), df_pot_la),
            pot_h=(("id", "time"), df_pot_h),
            nrj_leaf=(("id", "time", "leaf"), nrj_leaf)
        ),
        coords=dict(
            id=df_realized_la.index,
            time=dates,
            leaf=np.arange(n_leaf) + 1
        )
    )

    ds = xr.merge([ds, ds_archi]) 

    # Save the dataset to a NetCDF file
    today_str = date.today().strftime("%Y-%m-%d")
    os.makedirs(dir+f"{today_str}", exist_ok=True)  # noqa: PTH103
    ds.to_netcdf(dir+f"{today_str}/{filename}_{seed}.nc")


def concat_ds(files: list):
    '''Concatenate a list of datasets, ensuring the ids stay unique'''
    for i,file in enumerate(files):
        ds = xr.open_dataset(file)
        if i > 0:
            ds_new = ds.assign_coords(id=ds.id + offset)
            ds_tot = xr.concat([ds_tot, ds_new], dim="id")
            offset = int(ds_new.id.max()) + 2
        else:
            ds_tot = ds
            offset = int(ds.id.max()) + 2
    return ds_tot

