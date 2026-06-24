"""Some macro / utilities to run light simpulation on pgl/lpy virtual scene"""

from __future__ import annotations

import numpy as np
import pandas as pd
from itertools import chain
from openalea.caribu.CaribuScene import CaribuScene

from openalea.archicrop.display import build_scene
from openalea.archicrop.stand import compute_domain, config_crop_intercrop
from openalea.archicrop.stics_io import stics_weather_3d, stics_weather_3d_bis
from openalea.astk.sky_irradiance import sky_irradiance
from openalea.astk.sky_sources import caribu_light_sources, sky_sources

# Einc is the total incident energy received on the domain
# Eabs (float): the surfacic density of energy absorbed (MJ m⁻² s⁻¹)
# raint (float): from STICS, Photosynthetic Active Radiation intercepted by the canopy 	(MJ m⁻²) --> faPAR !!!

def illuminate(scene, light=None, domain=None, scene_unit='cm', labels=None, direct=True): # north=0):
    """Illuminate scene

    Args:
        scene: the scene (plantgl)
        light: lights. If None a vertical light is used
        pattern: the toric canopy pattern. If None, no pattern is used
        scene_unit: string indicating length unit in the scene (eg 'cm')
        north: the angle (deg, positive clockwise) from X+ to
         North (default: 0)

    Returns:

    """
    infinite = False
    if domain is not None:
        infinite = True
    # if light is not None:
    #     light = light_sources(*light, orientation=north)
    cs = CaribuScene(scene, light=light,scene_unit=scene_unit, pattern=domain, soil_mesh=1)
    raw, agg = cs.run(simplify=True, infinite=infinite, direct=direct) # direct=True
    df = pd.DataFrame(agg)  # noqa: PD901
    if labels is not None:
        labs = pd.DataFrame(labels)
        df = labs.merge(df.reset_index(), left_index=True, right_index=True)  # noqa: PD901
    return cs, raw['Eabs'], df


def mean_leaf_irradiance(df):
    df['Energy'] = df['Eabs'] * df['area']
    agg = df.loc[(df.is_green) & (df.label=='Leaf'),('plant','Energy','area')].groupby('plant').agg('sum')
    agg['Irradiance'] = agg['Energy'] / agg['area']
    return agg


def toric_canopy_pattern(dx=80, dy=5, density=None):
    if density is not None:
        if dx is not None:
            dy = 1. / density / (dx / 100.) * 100
        elif dy is not None:
            dx = 1. / density / (dy / 100.) * 100
        else:
            raise ValueError('At least one grid dimension (dx, dy) should be specified')  # noqa: EM101
    return (-0.5 * dx, -0.5 * dy,
             0.5 * dx, 0.5 * dy)


def compute_light_inter(scene, sky, pattern):
    
    # sky = str(data_path('Turtle16soc.light'))
    # zenith = str(data_path('zenith.light'))
    # opts = str(data_path('par_sorghum.opt'))
    # opts = list(map(str, [data_path('par_sorghum.opt'), data_path('nir_sorghum.opt')]))
    materials = {'par': (0.15, 0.06, 0.16, 0.03)}
    
    # complete set of files
    cs = CaribuScene(scene=scene, # give mtg rather than scene
                     light=sky, 
                     pattern=pattern, # for toric scene
                     scene_unit='cm',
                     soil_mesh=1, # 1
                     opt=materials) 
    
    raw,agg=cs.run(simplify=True, infinite=True)

    Qi, Qem, Einc = cs.getIncidentEnergy()
    # Qi is the total horizontal irradiance emitted by sources (m-2)
    # Qem is the sum of the normal-to-the-sources irradiance emitted by sources (m-2)
    # Einc is the total incident energy received on the domain

    # Qi_soil, Einc_soil = cs.getSoilEnergy()
    
    scene_eabs,values_eabs = cs.plot(raw['Eabs'],display=False) 

    v99_eabs = np.percentile(values_eabs, 99)
    nvalues_eabs = np.array(values_eabs)
    mask_eabs = nvalues_eabs > v99_eabs  
    nvalues_eabs[mask_eabs] = v99_eabs 
    sum_values_eabs = np.sum(nvalues_eabs)

    # Division of sums
    result = sum_values_eabs / Einc if Einc != 0 else 0 
    # result = (sum_values_ei - Einc_soil) / Einc if Einc != 0 else 0 

    return result  # noqa: RET504


# Viewer.display(scene_eabs, group_by_color=False, property=values_eabs)
# PlantGL(scene, group_by_color=False, property=values)


# location = {  
#     'longitude': 3.87,
#     'latitude': 45,
#     'altitude': 56,
#     'timezone': 'Europe/Paris'}

def light_interception(weather_file, daily_dynamics, density, location, mtgs, zenith=False, direct=False, save_scenes=False, inter_row=0.7):
    '''Compute light interception on plants 
    Args:
        weather_file: path to the weather file
        daily_dynamics: daily dynamics from STICS
        density: sowing density
        location: dictionary with location parameters (longitude, latitude, altitude, timezone)
        mtgs: list of MTGs for each plant
        zenith: if True, use zenith light sources
        save_scenes: if True, save the scenes as images
    Returns:
        nrj_per_plant: list of energy per plant
    '''
    # Read weather data
    df_weather = stics_weather_3d(filename=weather_file, daily_dynamics=daily_dynamics)

    # Sowing pattern
    domain = compute_domain(density = density, inter_row = inter_row) # cm

    # Define incident PAR
    incident_rads = list(df_weather.itertuples())

    # Compute light interception for each plant at each time step
    nrj_per_plant = {}
    # For each plant
    for k, mtgs_plant in mtgs.items():
        if mtgs_plant[0] is None:
            nrj_per_plant[k] = [None] * len(incident_rads)
        else:
            nrj_per_leaf = []
            # For each time step
            for i,(mtg, incident_rad) in enumerate(zip(mtgs_plant, incident_rads)):
                par = incident_rad.rad * 0.48 * 0.95
                # Compute light sources
                if zenith:
                    lights = [(par,(0,0,-1))]
                else:
                    irr = sky_irradiance(daydate=incident_rad.daydate, day_ghi=par, **location)
                    sun, sky = sky_sources(sky_type='clear_sky', sky_irradiance=irr, scale='global') #, source_irradiance='horizontal')
                    lights = caribu_light_sources(sun, sky)
                # Build and illuminate scene
                scene, labels = build_scene(mtg, senescence=True)
                cs, raw, agg = illuminate(scene=scene, light=lights, labels=labels, domain=domain, direct=direct) # --> cf PARaggregators in caribu scene node
                # Compute energy per leaf
                df_mod = mean_leaf_irradiance(agg) 
                nrj_per_leaf.append(agg.loc[(agg.is_green) & (agg['label'] == 'Leaf')]['Energy'].values)
                # Save scene if required
                if save_scenes:
                    scene_tmp = cs.plot(raw, display=False)[0]
                    scene_tmp.save(f'scene_{i}.png') # not as images !!!

            # nrj_per_plant[k] = [sum(nrj) for nrj in nrj_per_leaf]
            nrj_per_plant[k] = nrj_per_leaf

    '''
    # Calculate energy per leaf and irradiance per plant
    nrj_per_leaf = []
    # irr_per_plant = []

    for illuminated_plant in par_caribu:
        nrj_tmp = []
        # irr_tmp = []
        for df_scene in illuminated_plant:
            df_mod = mean_leaf_irradiance(df_scene)  # noqa: F841
            nrj_tmp.append(df_scene.loc[df_scene['label'] == 'Leaf']['Energy'].values)
            # irr_tmp.append(df_mod['Irradiance'].values[0])  
        nrj_per_leaf.append(nrj_tmp)
        # irr_per_plant.append(irr_tmp)

    nrj_per_plant = [[sum(growing_plant) for growing_plant in plant] for plant in nrj_per_leaf] # to dict !!!!
    '''

    return nrj_per_plant



def light_interception_IC(weather_file, daily_dynamics, stand, location, mtgs, zenith=False, direct=False, save_scenes=False, conv_coef=100):
    '''Compute light interception on plants 
    Args:
        weather_file: path to the weather file
        daily_dynamics: daily dynamics from STICS
        density: sowing density
        location: dictionary with location parameters (longitude, latitude, altitude, timezone)
        mtgs: list of MTGs for each plant
        zenith: if True, use zenith light sources
        save_scenes: if True, save the scenes as images
    Returns:
        nrj_per_plant: list of energy per plant
    '''

    # Compute light interception for each plant at each time step
    nrj_per_plant = {}
    list_of_graphs = {}
    positions_crop = {}
    domain_areas = {}
    # For each plant
    for a, usm in mtgs.items():
        nrj_per_plant[a] = {}
        list_of_graphs[a] = {}
        positions_crop[a] = {}
        domain_areas[a] = {}
        for b, algo in usm.items():
            plants = list(mtgs[a][b].keys())
            ids = list(mtgs[a][b][plants[0]].keys())
            nrj_per_plant[a][b] = {plants[0]:{},
                                   plants[1]:{}}
            nrj_per_leaf = {plants[0]:{},
                            plants[1]:{}}

            dates1 = list(mtgs[a][b][plants[0]][ids[0]].keys())
            dates2 = list(mtgs[a][b][plants[1]][ids[0]].keys())
            dates = sorted(list(dict.fromkeys(dates1 + dates2)))
            
            list_of_graphs[a][b], positions_crop[a][b] = config_crop_intercrop(dates, growing_plant_1=mtgs[a][b][plants[0]][ids[0]], growing_plant_2=mtgs[a][b][plants[1]][ids[0]], **stand[a][b])
            # print(a, positions_crop[a][b])

            domain = ((0, 0), (stand[a][b]["width"] * conv_coef, stand[a][b]["length"] * conv_coef))
            domain_area = (abs(domain[1][0] - domain[0][0])
                            / conv_coef
                            * abs(domain[1][1] - domain[0][1])
                            / conv_coef
                            )
            domain_areas[a][b] = domain_area

            # Read weather data
            df_weather = stics_weather_3d_bis(filename=weather_file, dates=dates)
            # Define incident PAR
            incident_rads = list(df_weather.itertuples())

            # For each time step
            for i,incident_rad in zip(dates,incident_rads):
                
                # Compute light sources
                par = incident_rad.rad * 0.48 * 0.95
                if zenith:
                    lights = [(par,(0,0,-1))]
                else:
                    irr = sky_irradiance(daydate=incident_rad.daydate, day_ghi=par, **location)
                    sun, sky = sky_sources(sky_type='clear_sky', sky_irradiance=irr, scale='global') #, source_irradiance='horizontal')
                    lights = caribu_light_sources(sun, sky)

                # print("lights :", lights)
                # print("mtg :", list_of_graphs[a][b][i])
                # print("positions :", positions_crop[a][b])
                # print("domain :", domain)

                # Build and illuminate scene
                scene, labels = build_scene(mtg=list_of_graphs[a][b][i], position=positions_crop[a][b], senescence=True)
                cs, raw, agg = illuminate(scene=scene, light=lights, labels=labels, domain=domain, direct=direct) # --> cf PARaggregators in caribu scene node
                
                # Compute energy per leaf
                df_mod = mean_leaf_irradiance(agg) 

                if stand[a][b]["nb_rows_1"] == 0 or stand[a][b]["nb_rows_2"] == 0:
                    for plant in [0,1]:
                        nrj_per_leaf[plants[plant]] = {i : (sum(list(agg.loc[(agg.is_green) & (agg['label'] == 'Leaf') & (agg['plant'] == plant)]['Energy'].values)))}
                else:
                    nrj_per_leaf[plants[0]] = {i : (sum(list(agg.loc[(agg.is_green) & (agg['label'] == 'Leaf') & (agg['plant'].isin(range(0,stand[a][b]["nb_rows_1"])))]['Energy'].values)) / stand[a][b]["nb_rows_1"])}
                    nrj_per_leaf[plants[1]] = {i : (sum(list(agg.loc[(agg.is_green) & (agg['label'] == 'Leaf') & (agg['plant'].isin(range(stand[a][b]["nb_rows_1"],stand[a][b]["nb_rows_2"]+stand[a][b]["nb_rows_1"])))]['Energy'].values)) / stand[a][b]["nb_rows_2"])}
                    #  nrj_per_leaf[plants[0]].append(agg.loc[(agg.is_green) & (agg['label'] == 'Leaf') & (agg['plant'].isin(range(0,stand[a][b]["nb_rows_1"])))]['Energy'].values)
                    #  nrj_per_leaf[plants[1]].append(agg.loc[(agg.is_green) & (agg['label'] == 'Leaf') & (agg['plant'].isin(range(stand[a][b]["nb_rows_1"],stand[a][b]["nb_rows_2"])))]['Energy'].values)
                
                # Save scene if required
                if save_scenes:
                    values = list(chain.from_iterable(raw.values()))
                    v99 = np.percentile(values, 99)
                    nvalues=np.array(values)
                    nvalues[nvalues>v99]=v99
                    raw_new = nvalues.tolist()
                    scene_tmp = cs.plot(raw, maxval=max(raw_new), display=False)[0]
                    scene_tmp.save(f'scene_{i}.png') # not as images !!!

            # nrj_per_plant[k] = [sum(nrj) for nrj in nrj_per_leaf]
            nrj_per_plant[a][b][plants[0]][ids[0]] = nrj_per_leaf[plants[0]]
            nrj_per_plant[a][b][plants[1]][ids[0]] = nrj_per_leaf[plants[1]]

    return nrj_per_plant, domain_areas
