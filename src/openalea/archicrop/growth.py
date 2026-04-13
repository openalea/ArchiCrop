from __future__ import annotations

from scipy.interpolate import splev

from openalea.mtg.traversal import pre_order2_with_filter

from .cereal_plant import CerealsVisitor
from .geometry import compute_growing_organ_geometry
from .turtle import CerealsTurtle


def init_visible_variables(g, daily_dynamics):
    """Initialize leaf and stem visible lengths"""

    for k in g.properties()["visible_length"]:
        g.properties()["visible_length"][k] = 0.0
    for k in g.properties()["leaf_lengths"]:
        g.properties()["leaf_lengths"][k] = [0.0]*len(daily_dynamics)
        g.properties()["leaf_areas"][k] = [0.0]*len(daily_dynamics)
        g.properties()["senescent_lengths"][k] = [0.0]*len(daily_dynamics)
    # for k in g.properties()["visible_leaf_area"]:
    #     g.properties()["visible_leaf_area"][k] = 0.0
    for k in g.properties()["leaf_lengths"]:
        g.properties()["stem_lengths"][k] = [0.0]*len(daily_dynamics)
    for k in g.properties()["stem_diameter"]:
        g.properties()["stem_diameter"][k] = 0.0

    return g


def equal_dist(increment, growing_organs):
    '''return initial distribution per organ (H: same amount for all organs)'''
    return {vid: increment/len(growing_organs) for vid in growing_organs}


def demand_dist(increment, growing_organs):
    '''return distribution of increment per organ proportionnal to potential dimension, i.e. length or area'''
    sum_growing_organs = sum([value["potential"] for value in growing_organs.values() if value["axis_order"]==0])
    return {vid: increment*values["potential"]/sum_growing_organs 
            for vid, values in growing_organs.items()}


def demand_dist_bis(increment, growing_organs):
    '''return distribution of increment per organ proportionnal to potential dimension, i.e. length or area'''
    sum_growing_organs = sum([value["potential"] - value["visible"] for value in growing_organs.values() if value["axis_order"]==0])
    return {vid: increment*(values["potential"] - values["visible"])/sum_growing_organs 
            for vid, values in growing_organs.items()}


def dev_dist(increment, growing_organs):
    '''return distribution of increment per organ proportionnal to development advancement'''
    sum_growing_organs = sum([1/value["age"] if value["age"] > 0. else 1 for value in growing_organs.values() if value["axis_order"]==0])
    return {vid: increment*(1/values["age"]/sum_growing_organs) if len(growing_organs) > 1 else increment
            if values["age"] > 0. else increment
            for vid, values in growing_organs.items()}


def age_dist(increment, senescing_organs):
    '''return distribution of increment per organ proportionnal to age'''
    sum_senescing_organs = sum([value["age"] for value in senescing_organs.values()])
    return {vid: increment*values["age"]/sum_senescing_organs
            for vid, values in senescing_organs.items()}


def get_growing_and_senescing_organs(g, time, prev_time):
    """Identify growing organs and their potential dimension at a given time"""

    growing_internodes = {}
    growing_leaves = {}
    senescing_leaves = {}

    # For each organ in the plant
    for vid, ml in g.property("mature_length").items(): 
        n = g.node(vid)
        # Update organ age
        n.age = max(0.0, time - n.start_tt)
        # Update stem diameter 
        n_stem = n.parent() if n.label.startswith("Leaf") else n
        n.stem_diameter = min(n.mature_stem_diameter/2 * (0.5+0.5*(time - n_stem.start_tt) / (n_stem.end_tt - n_stem.start_tt)), n.mature_stem_diameter) 
        # If it is a stem element / internode
        if n.label.startswith("Stem"): 
            # If the internode is growing, or has finished growing in the last time step, according to development
            if (n.start_tt < time <= n.end_tt or prev_time < n.end_tt <= time) and n.visible_length < ml: 
                growing_internodes[vid] = {"potential": ml, "visible": n.visible_length, "age": n.age, "axis_order": g.order(g.parent(vid))}
        # If it is a leaf
        elif n.label.startswith("Leaf"):
            # Update leaf inclination 
            n.inclination = min(1.5*(time - n.start_tt) / (n.end_tt - n.start_tt), 1)
            # If the leaf is growing, or has finished growing in the last time step, according to development
            if (n.start_tt < time <= n.end_tt or prev_time < n.end_tt <= time) and n.visible_leaf_area < n.leaf_area:
                growing_leaves[vid] = {"potential": n.leaf_area, "visible": n.visible_leaf_area, "age": n.age, "axis_order": 0}
            # If the leaf is senescing, according to development
            if n.senescence <= time and n.visible_leaf_area > n.senescent_area: # not n.dead: # and n.srt > 0: 
                senescing_leaves[vid] = {"potential": n.visible_leaf_area, "visible": n.senescent_area, "age": n.age, "axis_order": 0}
        
    return growing_internodes, growing_leaves, senescing_leaves


def distribute_among_organs(g, growing_organs, day, time, time_increment, increment_to_distribute, distribution_function, grow_function):
    """Distribute increment among growing organs up to potential of each organ"""

    # Create dict {id : 0.0 for id in growing organs}
    increment_for_each_organ = dict.fromkeys(growing_organs, 0.0)

    # While there are growing organs (i.e. organs that have not reached their theoretical dimension),
    # and while there is still growth increment to distribute among growing organs
    while len(growing_organs) > 0 and increment_to_distribute > 1e-5: 
        # Compute dict {id : growth increment for id in growing organs},
        # growth increment for each organ being computed according to distribution function from plant-scale growth increment
        incr_temp = distribution_function(increment_to_distribute, growing_organs)
        # Re-initialize increment to distribute at each iteration
        increment_to_distribute = 0.0
        # Iterate over growing organs, in theory from lowest to highest ids
        for vid, val in incr_temp.items():
            # For each growing organ
            if vid in growing_organs:
                n = g.node(vid)
                if n.label == 'Leaf' or n.label == 'Stem':
                    # Set potential increment for organ growth to the minimal value between 
                    # the increment computed for this organ and the remaining bit to growth to reach potential
                    potential_increment = min(val, growing_organs[vid]["potential"] - growing_organs[vid]["visible"])
                    # Store the growth increment of the organ in a dict
                    increment_for_each_organ[vid] += potential_increment
                    # Update the dict of growing organs, adding the increment to the visible (i.e. current) dimension
                    growing_organs[vid]["visible"] += potential_increment
                    # Increment to distribute during the next iteration among the organs that have not reached their potential yet
                    # Not equal to zero if the organ has reached its potential
                    increment_to_distribute += val - potential_increment
                # elif n.label == 'Stem' and vid in g.components(main_axis):

        # Filter growing_organs to keep only the ones that have not reached their potential dimension
        growing_organs = {vid: {"potential": values["potential"], "visible": values["visible"], "age": values["age"]} 
                          for vid, values in growing_organs.items() 
                          if values["visible"] < values["potential"]}

    # Update organ geometry in MTG  
    for vid, incr in increment_for_each_organ.items():
        n = g.node(vid)
        grow_function(n, incr, day, time_increment)


def grow_organs(g, day, daily_dynamics, rate=False, distribution_function=demand_dist):
    """Grow organs according to plant-scale growth increment on a given day"""

    time = daily_dynamics["Thermal time"]
    thermal_time_incr = daily_dynamics["Thermal time increment"]
    prev_time = max(0.0, time - thermal_time_incr)

    # Get growing and senescent organs in a dict of dict : {id : dict(potential, visible, age)}
    growing_internodes, growing_leaves, senescing_leaves = get_growing_and_senescing_organs(g, time, prev_time)

    # Distribute daily plant-scale growth increment among organs
    if daily_dynamics["Leaf area increment"] > 0.0 and len(growing_leaves) > 0:
        update_cereal_leaf_growth = update_cereal_leaf_growth_rate if rate else update_cereal_leaf_growth_area
        distribute_among_organs(g=g, 
                                growing_organs=growing_leaves, 
                                day=day,
                                time=time,
                                time_increment=thermal_time_incr,
                                increment_to_distribute=daily_dynamics["Leaf area increment"], 
                                distribution_function=distribution_function, 
                                grow_function=update_cereal_leaf_growth)

    if daily_dynamics["Height increment"] > 0.0 and len(growing_internodes) > 0:
        update_stem_growth = update_stem_growth_rate if rate else update_stem_growth_height
        distribute_among_organs(g=g, 
                                growing_organs=growing_internodes, 
                                day=day,
                                time=time,
                                time_increment=thermal_time_incr,
                                increment_to_distribute=daily_dynamics["Height increment"], 
                                distribution_function=distribution_function,
                                grow_function=update_stem_growth)

    if daily_dynamics["Senescent leaf area increment"] > 0.0 and len(senescing_leaves) > 0:
        distribute_among_organs(g=g, 
                                growing_organs=senescing_leaves, 
                                day=day,
                                time=time,
                                time_increment=thermal_time_incr,
                                increment_to_distribute=daily_dynamics["Senescent leaf area increment"], 
                                distribution_function=age_dist,
                                grow_function=update_leaf_senescence_area)

'''
Detailed process of leaf elongation:
-	At the beginning of the simulation: 
        - From the normalized shape of the blade along the midrib scaled for each leaf according to the potential area, 
          compute points on that curve, to have a correspondence between S and l.
	    - Integrate l=f(S) with python module scipy.integrate.cumulated_simpson.
	    - Fit a spline of degree k=3 using python module scipy.interpolate.splrep.
-	When growing or senescing: 
        - Evaluate normalized length l for a normalized leaf area S+dS (i.e. percent of potential leaf area reached).
        - Scale l according to potential total length L.
'''

def update_cereal_leaf_growth_area(n, LA_for_this_leaf, day, time_increment=None):
    """Update the visible leaf area and length of a cereal leaf, 
    under a potential leaf area."""
    if n.visible_leaf_area + LA_for_this_leaf < n.leaf_area:
        n.visible_leaf_area += LA_for_this_leaf
        relative_visible_area = n.visible_leaf_area / n.leaf_area
        n.visible_length = float(splev(x=relative_visible_area, tck=n.tck)) * n.mature_length
    else: # warning
        n.visible_leaf_area = n.leaf_area
        n.visible_length = n.mature_length

    for d in range(day-1, len(n.leaf_areas)):
        n.leaf_areas[d] = n.visible_leaf_area
        n.leaf_lengths[d] = n.visible_length
    

def update_cereal_leaf_growth_rate(n, LA_for_this_leaf, day, time_increment):
    """Update the visible leaf area and length of a cereal leaf, 
    under a potential growth rate."""
    rate = LA_for_this_leaf / time_increment
    if rate <= n.potential_growth_rate:
        n.visible_leaf_area += LA_for_this_leaf
    else:
        new_LA_for_this_leaf = n.potential_growth_rate * time_increment
        n.visible_leaf_area += new_LA_for_this_leaf
    relative_visible_area = n.visible_leaf_area / n.leaf_area
    n.visible_length = float(splev(x=relative_visible_area, tck=n.tck)) * n.mature_length

    for d in range(day-1, len(n.leaf_areas)):
        n.leaf_areas[d] = n.visible_leaf_area
        n.leaf_lengths[d] = n.visible_length


def update_stem_growth_height(n, height_for_this_internode, day, time_increment=None):
    """Update the visible length of a stem internode, 
    under a potential height."""
    if n.visible_length + height_for_this_internode < n.mature_length:
        n.visible_length += height_for_this_internode
    else:
        n.visible_length = n.mature_length

    for d in range(day-1, len(n.stem_lengths)):
        n.stem_lengths[d] = n.visible_length


def update_stem_growth_rate(n, height_for_this_internode, day, time_increment):
    """Update the visible length of a stem internode, 
    under a potential growth rate."""
    rate = height_for_this_internode / time_increment
    if rate <= n.potential_growth_rate:
        n.visible_length += height_for_this_internode
    else:
        new_height_for_this_internode = n.potential_growth_rate * time_increment
        n.visible_length += new_height_for_this_internode

    for d in range(day-1, len(n.stem_lengths)):
        n.stem_lengths[d] = n.visible_length


def update_leaf_senescence_area(n, sen_LA_for_this_leaf, day, time_increment=None):
    """Update the senescent area and length of a cereal leaf, 
    that cannot be greater than the visible leaf area."""
    n.grow = False
    if n.senescent_area + sen_LA_for_this_leaf < n.visible_leaf_area:
        n.senescent_area += sen_LA_for_this_leaf
        relative_senescent_area = n.senescent_area / n.visible_leaf_area
        n.srt = 1 - float(splev(x=relative_senescent_area, tck=n.tck))
        n.senescent_length = (1 - n.srt) * n.visible_length
    else:
        n.senescent_area = n.visible_leaf_area
        n.senescent_length = n.visible_length
        n.dead = True

    for d in range(day, len(n.senescent_lengths)):
        n.senescent_lengths[d] = n.senescent_length


def mtg_turtle_time_with_constraint(g, day, daily_dynamics, rate=True, distribution_function=demand_dist, update_visitor=None):  # noqa: ARG001
    """Compute the geometry on each node of the MTG using Turtle geometry.

    Update_visitor is a function called on each node in a pre order (parent before children).
    This function allow to update the parameters and state variables of the vertices.

    """

    g.properties()["geometry"] = {}
    g.properties()["_plant_translation"] = {}

    max_scale = g.max_scale()

    time = daily_dynamics["Thermal time"]
    grow_organs(g, day, daily_dynamics, rate=rate, distribution_function=distribution_function)
    
    cereal_visitor = CerealsVisitorGrowth(False)
    
    def traverse_with_turtle_time(g, vid, time, visitor=cereal_visitor):
        turtle = CerealsTurtle()

        def push_turtle(v):
            n = g.node(v)
            # if 'Leaf' in n.label:
            #    return False
            try:
                start_tt = n.start_tt
                if start_tt > time:
                    return False
            except:  # noqa: E722
                pass
            if g.edge_type(v) == "+":
                turtle.push()
            return True

        def pop_turtle(v):  # noqa: RET503
            n = g.node(v)
            try:
                start_tt = n.start_tt
                if start_tt > time:
                    return False
            except:  # noqa: E722
                pass
            if g.edge_type(v) == "+":
                turtle.pop()

        if g.node(vid).label.startswith("Leaf") or g.node(vid).label.startswith("Stem"):  # noqa: SIM102
            if g.node(vid).start_tt <= time:
                visitor(g, vid, turtle) 
                # turtle.push()

        for v in pre_order2_with_filter(g, vid, None, push_turtle, pop_turtle):
            if v == vid:
                continue
            # Done for the leaves
            if g.node(v).start_tt > time:
                continue
            visitor(g, v, turtle) 

        # scene = turtle.getScene()
        return g


    for plant_id in g.component_roots_at_scale_iter(g.root, scale=max_scale):
        g = traverse_with_turtle_time(g, plant_id, time)

    return g



class CerealsVisitorGrowth(CerealsVisitor):
    def __init__(self, classic):
        super().__init__(classic)

    def __call__(self, g, v, turtle): 
        
        # 1. retrieve the node
        geoms_senesc = g.property("geometry_senescent")
        n = g.node(v)

        # Go to plant position if first plant element
        if n.parent() is None:
            turtle.move(0, 0, 0)
            # initial position to be compatible with canMTG positioning
            turtle.setHead(0, 0, 1, -1, 0, 0)

        # Manage inclination of tiller
        if g.edge_type(v) == "+" and not n.label.startswith("Leaf"):
            # axis_id = g.complex(vid); g.property('insertion_angle')
            # TODO : vary as function of age and species(e.g. rice)
            angle = 2*n.tiller_angle if g.order(v) == 1 else n.tiller_angle 
            turtle.down(angle)
            turtle.elasticity = n.gravitropism_coefficient 
            turtle.tropism = (0, 0, 1)


        # incline turtle at the base of stems,
        if n.label.startswith("Stem"):
            azim = float(n.azimuth) if n.azimuth else 0.0
            if azim:
                turtle.rollL(azim)

        # Try to build the mesh from GC rather than Cylinder
        if n.label.startswith("Leaf") or n.label.startswith("Stem"):
            # update geometry of elements
            mesh = compute_growing_organ_geometry(n, self.classic)
            if mesh:  # To DO : reset to None if calculated so ?
                n.geometry = turtle.transform(mesh)
                n.anchor_point = turtle.getPosition()
                n.heading = turtle.getHeading()
            if v in geoms_senesc and geoms_senesc[v] is not None:
                geoms_senesc[v] = turtle.transform(geoms_senesc[v])

        # 3. Update the turtle and context
        turtle.setId(v)
        if n.label.startswith("Stem"):
            if n.visible_length > 0:
                turtle.f(n.visible_length)
            turtle.context.update({"top": turtle.getFrame()})
        if n.label.startswith("Leaf"):  # noqa: SIM102
            if n.lrolled > 0:
                turtle.f(n.lrolled)
                turtle.context.update({"top": turtle.getFrame()})


