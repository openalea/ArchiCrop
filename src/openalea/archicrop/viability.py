import math

import numpy as np
from scipy.interpolate import splev, splrep

from .cereal_axis import bell_shaped_dist, geometric_dist
from .cereal_leaf import growing_leaf_area
from .stics_io import get_pheno


def leaf_area_plant(g):
    '''Compute leaf area of a plant from its MTG, at any growing stage'''
    S = 0
    for k,leaf in g.properties()["shape"].items():
        S += growing_leaf_area(leaf, g.properties()["visible_length"][k], g.properties()["mature_length"][k], g.properties()["shape_max_width"][k])
    return S


def resolve_organ_growth(N, ligul_factor, la_ends):
    '''Matrix equation AS=B to resolve leaf area S assuming linear leaf growth 
    and following constrained plant growth B, 
    for a single stem cereal'''

    # Initialize matrices
    A = np.zeros((N, N))
    B = np.zeros(N)

    # Build matrices A and B
    for x in range(0, N-1):
        A[x, x] = 1
        A[x, x+1] = (ligul_factor - 1) / ligul_factor # upper diagonal = advancement of leaf x+1 at end of growth of leaf x
        for i in range(x):
            A[x, i] = 1 # lower triangle = fully grown leaves 0 to x
        B[x] = la_ends[x] 

    # Last line of A and B
    A[N-1, N-1] = 1
    for i in range(N-1):
        A[N-1, i] = 1
    B[N-1] = la_ends[N-1] 

    # Resolution
    S = np.linalg.solve(A, B)

    return S


def compute_skew(leaf_areas, rank, nb_phy, rmax):
    '''Comute skew parameter of bell shape, so that the curve passes through leaf area of rank rank.'''
    return math.exp(math.log(leaf_areas[rank-1]/max(leaf_areas))/(2*(rank/nb_phy - rmax)**2 + (rank/nb_phy - rmax)**3))


def compute_viable_params(params_sets: dict, daily_dynamics: dict, pot_factor_lai: float = 1.5, pot_factor_height: float = 3) -> dict:
    """Compute viable parameters wrt the dynamic constraint of LAI and height."""

    # Dynamics of vegetative phase
    index_end_veg, end_veg, _ = get_pheno(daily_dynamics)
    thermal_time = [value["Thermal time"] for value in daily_dynamics.values()][:index_end_veg+2]
    leaf_area_plant = [value["Plant leaf area"] for value in daily_dynamics.values()][:index_end_veg+2]
    
    # Viable values for parameters
    new_params_sets = {}
    for id,params in params_sets.items(): 
        # Define leaf elongation duration and nb of phytomers
        leaf_duration = params["leaf_duration"] 
        nb_phy = params["nb_phy"]
        # Computes phyllochron
        # phyllochron = (end_veg-thermal_time[0])/(nb_phy + leaf_duration)
        phyllochron = (end_veg)/(nb_phy - 1 + leaf_duration)

        if min(params["phyllochron"]) <= phyllochron <= max(params["phyllochron"]): # Check that the value of skew is within an acceptable range for the species
            params_sets[id]["phyllochron"] = phyllochron

            # Compute organ development
            # starts = [i * phyllochron + thermal_time[0] for i in range(nb_phy)]
            starts = [i * phyllochron for i in range(nb_phy)]
            ends = [start + phyllochron * leaf_duration for start in starts]

            # Interpolate growth dynamics and compute it at given points (instead of daily)
            spl_la = splrep(thermal_time, leaf_area_plant)
            la_ends = splev(ends, spl_la)

            # Numeric solution for minimal constrained growth (H: linear organ growth)
            min_leaf_areas = resolve_organ_growth(nb_phy, leaf_duration, la_ends)

            # Find viable rmax interval
            for i, la in enumerate(min_leaf_areas):
                if la == max(min_leaf_areas):
                    id_max = i+1
                    break
            rmax_int = np.linspace(max(0,max((id_max-1)/nb_phy, params['rmax'][0])), min(1,min((id_max+1)/nb_phy, params['rmax'][1])), 10)
        
            # Find viable (rmax, skew) pairs
            leaf_areas_norm = [la/max(min_leaf_areas) for la in min_leaf_areas]
            skews_rmax_ok = []
            for rmax in rmax_int: 
                for rank in range(1, nb_phy+1): 
                    if rmax != rank/nb_phy and min_leaf_areas[rank-1] > 0:
                        skew = compute_skew(rank=rank, nb_phy=nb_phy, rmax=rmax, leaf_areas=leaf_areas_norm) # Compute skew parameter, so that the bell shape curve passes through leaf area of rank rank
                        if params['skew'][0] < skew < params['skew'][1]: # Check that the value of skew is within an acceptable range
                            ok = True
                            bell_shaped_leaf_areas = bell_shaped_dist(leaf_area_plant[-1]*pot_factor_lai, nb_phy, rmax, skew) # Compute the bell shape with given (rmax, skew) and constrained plant-scale leaf area
                            # Verify that the bell shape is greater or equal to the minimal solution 
                            for i,bs in enumerate(bell_shaped_leaf_areas):
                                if bs <= min_leaf_areas[i]: 
                                    ok = False
                                    break
                            if ok:
                                skews_rmax_ok.append((skew, rmax))
                    
            # Build new parameter set with updated 'skew' and 'rmax'
            for (s,r) in skews_rmax_ok[:1]: #################### All realized leaf area distributions are the same no matter the (rmax,skew) given
                new_param = {}
                for key, value in params.items():
                    if key == "skew":
                        new_param[key] = s
                    elif key == "rmax":
                        new_param[key] = r
                    else:
                        new_param[key] = value
                new_params_sets[len(new_params_sets)+1] = new_param

    return new_params_sets
