from __future__ import annotations

import math

import numpy as np


def geometric_dist(height, nb_phy, q, u0):
    """
    Calculates the heights of leaves'ligules along an axis, according to a geometric series,
    starting from an offset height (u0) = pseudostem height.

    Parameters:
    - height (float): Total height of the plant.
    - nb_phy (int): Number of phytomers (leaves).
    - q (float): Geometric progression factor (controls spacing between leaves).
    - u0 (float): Offset height to start the geometric distribution.

    Returns:
    - List[float]: Normalized distances of leaves along the plant's height.
    """
    if nb_phy <= 0:
        msg = "Number of phytomers (nb_phy) must be greater than zero."
        raise ValueError(msg)
    if q <= 0:
        msg = "Geometric progression factor (q) must be positive."
        raise ValueError(msg)
    if u0 >= height:
        msg = "Offset height (u0) must be less than the total height."
        raise ValueError(msg)

    # Calculate the height available for geometric distribution
    # remaining_height = height - u0

    # Generate table of heights for geometric distribution
    heights = np.array([i*float(height) / nb_phy if q == 1 else height * (1 - q**i) / (1 - q**nb_phy) for i in range(1, nb_phy + 1)])
    
    # Add the offset height (u0) to each leaf's position
    # heights = u0 + heights_without_u0
    
    return heights.tolist()


# def bell_shaped_dist(max_leaf_length, nb_phy, rmax, skew):
#     """returns relative leaf area of individual leaves along bell shaped model"""
#     k = -np.log(skew) * rmax
#     r = np.linspace(1.0 / nb_phy, 1, nb_phy)
#     relative_area = np.exp(-k / rmax * (2 * (r - rmax) ** 2 + (r - rmax) ** 3))
#     return relative_area * max_leaf_length

def bell_shaped_dist(leaf_area, nb_phy, rmax, skew):
    """returns relative leaf area of individual leaves along bell shaped model, so that the sum equals leaf_area"""
    k = -np.log(skew) * rmax
    r = np.linspace(1.0 / nb_phy, 1, nb_phy)
    relative_area = np.exp(-k / rmax * (2 * (r - rmax) ** 2 + (r - rmax) ** 3))
    total_area = sum(relative_area)
    normalized_leaf_areas = [area / total_area for area in relative_area]
    return [leaf_area * la for la in normalized_leaf_areas]



