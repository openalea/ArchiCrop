
import numpy as np

def leaf_mass(SLA, leaf_area):
    '''Compute leaf mass from its area and SLA'''
    return leaf_area / SLA

def stem_mass(vol_mass, length, diameter):
    '''Compute stem mass from its volumetric mass and volume, computed with length and diameter'''
    volume = length * (diameter/2)**2 * np.pi
    return volume * vol_mass

def leaf_n_demand(SLA, N_conc_leaf, d_leaf_area):
    '''Compute N demand of a leaf from its area, SLA and N concentration'''
    return d_leaf_area * N_conc_leaf / SLA

def stem_n_demand(vol_mass, N_conc_stem, length, d_length, diameter, d_diameter):
    '''Compute N demand of a stem element from its volumetric mass, N concentration 
    and volume, computed with length and diameter'''
    d_volume_height = d_length * (diameter/2)**2 * np.pi
    d_volume_width = (length - d_length) * np.pi * ((diameter/2)**2 - ((diameter - d_diameter)/2)**2)
    return (d_volume_height + d_volume_width) * vol_mass * N_conc_stem

def leaf_n_content(SLA, N_conc_leaf, leaf_area):
    '''Compute N content of a leaf from its area, SLA and N concentration'''
    return leaf_area * N_conc_leaf / SLA

def stem_n_content(vol_mass, N_conc_stem, length, diameter):
    '''Compute N content of a stem element from its volumetric mass, N concentration 
    and volume, computed with length and diameter'''
    volume = length * (diameter/2)**2 * np.pi
    return volume * vol_mass * N_conc_stem

def compute_n_demand(g, dates):
    '''Compute nitrogen demand of all organs in an MTG'''

    organs = g.vertices(scale=g.max_scale())

    # Store plant-scale daily N demand
    plant_n_demand = {date: 0.0 for date in dates}
    plant_n_content = {date: 0.0 for date in dates}
    plant_leaf_mass = {date: 0.0 for date in dates}
    plant_stem_mass = {date: 0.0 for date in dates}

    # Loop over all organs in the MTG
    for organ in organs:
        # Get organ type
        n = g.node(organ)
        organ_type = n.label

        if organ_type == "Leaf":
            n.n_demands = [0.0]+[leaf_n_demand(n.SLA, n.N_conc_leaf, dla) 
                           for dla in np.diff(np.array(n.leaf_areas))]
        
        elif organ_type.startswith("Stem"):
            n.n_demands = [0.0]+[stem_n_demand(n.vol_mass, n.N_conc_stem, l, dl, d, dd)
                           for l, dl, d, dd in zip(np.array(n.stem_lengths),
                                                   np.diff(np.array(n.stem_lengths)), 
                                                   np.array(n.stem_diameters),
                                                   np.diff(np.array(n.stem_diameters)))]

        for i, nd in enumerate(n.n_demands):
            plant_n_demand[dates[i]] += nd
            if organ_type == "Leaf":
                plant_n_content[dates[i]] += leaf_n_content(n.SLA, n.N_conc_leaf, n.leaf_areas[i])
                plant_leaf_mass[dates[i]] += leaf_mass(n.SLA, n.leaf_areas[i])
            elif organ_type.startswith("Stem"):
                plant_stem_mass[dates[i]] += stem_mass(n.vol_mass, n.stem_lengths[i], n.stem_diameters[i])
                plant_n_content[dates[i]] += stem_n_content(n.vol_mass, n.N_conc_stem, n.stem_lengths[i], n.stem_diameters[i])

    return plant_n_demand, plant_n_content, plant_leaf_mass, plant_stem_mass


