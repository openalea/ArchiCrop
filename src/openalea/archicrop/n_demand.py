
import numpy as np

def leaf_n_demand(SLA, N_conc_leaf, d_leaf_area):
    '''Compute N demand of a leaf from its area, SLA and N concentration'''
    return d_leaf_area * N_conc_leaf / SLA

def stem_n_demand(vol_mass, N_conc_stem, length, d_length, diameter, d_diameter):
    '''Compute N demand of a stem element from its volumetric mass, N concentration 
    and volume, computed with length and diameter'''
    d_volume_height = d_length * (diameter/2)**2 * np.pi
    d_volume_width = (length - d_length) * np.pi * ((diameter/2)**2 - ((diameter - d_diameter)/2)**2)
    return (d_volume_height + d_volume_width) * vol_mass * N_conc_stem

def compute_organ_n_demand(g, dates):
    '''Compute nitrogen demand of all organs in an MTG'''

    organs = g.vertices(scale=g.max_scale())

    # Store plant-scale daily N demand 
    plant_n_demand = {date: 0.0 for date in dates}


    # Loop over all organs in the MTG
    for organ in organs:
        # Get organ type
        n = g.node(organ)
        organ_type = n.label

        if organ_type == "Leaf":
            n.n_demands = [leaf_n_demand(n.SLA, n.N_conc_leaf, dla) 
                           for dla in np.diff(np.array(n.leaf_areas))]

        elif organ_type.startswith("Stem"):
            n.n_demands = [stem_n_demand(n.vol_mass, n.N_conc_stem, l, dl, d, dd)
                           for l, dl, d, dd in zip(np.array(n.stem_lengths),
                                                   np.diff(np.array(n.stem_lengths)), 
                                                   np.array(n.stem_diameters),
                                                   np.diff(np.array(n.stem_diameters)))]

        for i, nd in enumerate(n.n_demands):
            plant_n_demand[dates[i]] += nd

    return plant_n_demand


