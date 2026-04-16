import math

import numpy as np
from scipy.interpolate import splev, splrep
from scipy.optimize import nnls
from openalea.mtg import MTG

from .cereal_axis import bell_shaped_dist, geometric_dist
from .cereal_leaf import growing_leaf_area
from .stics_io import get_pheno


def leaf_area_plant(g):
    '''Compute leaf area of a plant from its MTG, at any growing stage'''
    S = 0
    for k,leaf in g.properties()["shape"].items():
        S += growing_leaf_area(leaf, g.properties()["visible_length"][k], g.properties()["mature_length"][k], g.properties()["shape_max_width"][k])
    return S

def generate_tillers(nb_phy, phyllochron, duration, nb_tillers, tiller_delay=1):
    '''
    Generates a MTG of a plant with main stem and tillers, considering the delay of appearance of tillers and the timing of appearance and elongation of the phytomers.
    Return a dict with the order and the delay of each tiller.
    '''
    ranks = range(1, nb_phy + 1)
    ntop = max(ranks) - np.array(ranks) + 1
    
    tt = 0
    dtt = phyllochron * duration
    
    tiller_points = []
    rank_tillers = []
    
    g = MTG()
    # Add a root vertex for the plant
    vid_plant = g.add_component(g.root, label="Plant", edge_type="/") 
    # Add a plant vertex for the main axis
    vid_axis = g.add_component(vid_plant, label="MainAxis", edge_type="/")
 
    # iterate over the number of phytomers of the main stem
    first = True
    for rank in ranks:
        if first:
            vid = g.add_component(vid_axis)
            first = False
        else:
            vid = g.add_child(vid, edge_type="<")
        g.node(vid).start_tt = tt 
        g.node(vid).end_tt = tt + dtt
        tt += phyllochron
        tiller_points.append((vid, tt + tiller_delay*phyllochron))
    
    for i in range(nb_tillers):
        # add a tiller
        vid, time = tiller_points.pop(0)
    
        tillers = []
        axis_id = g.complex(vid)
        r = g.Rank(vid) + 1  # Number of edges from the root of the axis
        rank_tillers.append(r)
        n = len(g.Axis(vid))
        len_tiller = n - r - 1  # we remove the parent that does not belong to the tiller 
        nb_phy = len_tiller
    
        ranks = range(1, nb_phy + 1)
        ntop = max(ranks) - np.array(ranks) + 1
    
        tt = time
    
        tid = g.add_child(parent=axis_id, edge_type='+', label='Axis')
    
        first = True
        for rank in ranks:
            if first:
                vid, tid2 = g.add_child_and_complex(parent=vid, complex=tid, edge_type='+')
                first = False
            else:
                vid = g.add_child(parent=vid, edge_type='<')
            g.node(vid).start_tt = tt
            g.node(vid).end_tt = tt + dtt
            tt += phyllochron
            tillers.append((vid, tt + tiller_delay*phyllochron))
    
        tiller_points.extend(tillers)
        tiller_points.sort(key=lambda x: x[1]) # sorted by time
        tiller_points = tiller_points[:nb_tillers-i]

    axes = g.vertices(scale=2)  # scale=2: axes
    T = {i+1: [g.order(axis), g.order(axis)+r] for i,(axis,r) in enumerate(zip(axes[1:],rank_tillers))} # {id: (order, delay)}

    return g, T 

'''
def resolve_organ_growth(N, ligul_factor, la_ends):
    """Matrix equation AS=B to resolve leaf area S assuming linear leaf growth 
    and following constrained plant growth B, 
    for a single stem cereal"""

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
'''


def resolve_organ_growth(N, duration, la_ends, nb_tillers=0, tillers=None, reduction_factor=1):
    '''
    A: matrix of the advancement of elongation of phytomers i at the end of elongation of phytomers j on the main stem
    T: dict of matrices of the advancement of elongation of phytomers i at the end of elongation of phytomers j on tillers
    B: vector of constrained plant-scale dynamics evaluated at end of elongation of phytomers of the main stem
    Resolve (A+∑T)S=B
    '''
    
    A = np.zeros((N, N)) 
    B = np.zeros(N)
    if nb_tillers > 0:
        T = {i: np.zeros((N, N)) for i in range(1,nb_tillers+1)}
    else:
        T = {}

    for x in range(0, N):
        A[x, min(x+1,N-1)] = (duration - 1) / duration
        for i in range(x+1):
            A[x, i] = 1
        B[x] = la_ends[x] 
        
        if nb_tillers > 0:
            for id,t in T.items():
                delay = tillers[id][1]
                if x >= delay:
                    t[x, min(x+1-delay,N-1)] = (duration - 1) / duration
                    for j in range(x+1): # range(delay+1, x+1) for other hypothesis
                        t[x, max(0,j-delay)] = 1

    # Sum matrices of main stem and tillers 
    P = A
    if nb_tillers > 0:
        for id,t in T.items():
            t *= reduction_factor**tillers[id][0]
            P += t

    # S = np.linalg.solve(P, B)
    S = nnls(P, B)[0]

    return S, P, A, T



def compute_skew_old(leaf_areas, rank, nb_phy, rmax):
    '''Comute skew parameter of bell shape, so that the curve passes through leaf area of rank rank.'''
    return math.exp(math.log(leaf_areas[rank-1]/max(leaf_areas))/(2*(rank/nb_phy - rmax)**2 + (rank/nb_phy - rmax)**3))

def compute_skew(leaf_areas, rank, nb_phy, rmax):
    '''Comute skew parameter of bell shape, so that the curve passes through leaf area of rank rank.'''
    return math.log(leaf_areas[rank-1]/max(leaf_areas))/((rank/nb_phy - rmax)**2)


def compute_viable_params(params_sets: dict, daily_dynamics: dict, pot_factor_lai: float = 1.5, pot_factor_height: float = 3, end_veg_indefinite=None) -> dict:
    """Compute viable parameters wrt the dynamic constraint of LAI and height."""

    # Dynamics of vegetative phase
    thermal_time = [value["Thermal time"] for value in daily_dynamics.values()]
    if end_veg_indefinite:
        index_end_veg = end_veg_indefinite
        end_veg = thermal_time[index_end_veg]
    else:
        index_end_veg, end_veg, _ = get_pheno(daily_dynamics)
    thermal_time = [value["Thermal time"] for value in daily_dynamics.values()][:index_end_veg+2]
    leaf_area_plant = [value["Plant leaf area"] for value in daily_dynamics.values()][:index_end_veg+2]
    
    # Viable values for parameters
    new_params_sets = {}
    for id,params in params_sets.items(): 
        # Define leaf elongation duration and nb of phytomers
        leaf_duration = params["leaf_duration"] 
        nb_phy = params["nb_phy"]
        nb_tillers = params_sets[id]["nb_tillers"]
        reduction_factor = params_sets[id]["reduction_factor"]
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

            g, tillers = generate_tillers(nb_phy, phyllochron, leaf_duration, nb_tillers) # Dict for tillers: {id: [order, delay]}

            # Numeric solution for minimal constrained growth (H: linear organ growth)
            min_leaf_areas,_,_,_ = resolve_organ_growth(nb_phy, leaf_duration, la_ends, nb_tillers, tillers, reduction_factor)

            # Find viable rmax interval
            for i, la in enumerate(min_leaf_areas):
                if la == max(min_leaf_areas):
                    id_max = i+1
                    break
            rmax_int = np.linspace(max(0,min(max((id_max-1)/nb_phy, params['rmax'][0]), params['rmax'][1])), min(1,max(min((id_max+1)/nb_phy, params['rmax'][1]), params['rmax'][0])), 10)
            print(rmax_int)

            # Find viable (rmax, skew) pairs
            leaf_areas_norm = [la/max(min_leaf_areas) for la in min_leaf_areas]
            skews_rmax_ok = []
            for rmax in rmax_int: 
                for rank in range(1, nb_phy+1): 
                    if rmax != rank/nb_phy and min_leaf_areas[rank-1] > 0:
                        skew = compute_skew(rank=rank, nb_phy=nb_phy, rmax=rmax, leaf_areas=leaf_areas_norm) # Compute skew parameter, so that the bell shape curve passes through leaf area of rank rank
                        if min(params['skew']) <= skew <= max(params['skew']): # Check that the value of skew is within an acceptable range
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
