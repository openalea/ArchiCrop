from __future__ import annotations

import os
from datetime import date

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.lines import Line2D


def plot_constained_vs_pot(dates, pot_la, pot_h, leaf_area_plant, height_canopy, density, stics_color="orange", archicrop_color="green"):
    
    # conversion factor
    cf_cm = 100

    fig, axes = plt.subplots(2, 1, figsize=(12, 6), sharex=True)  # 1 row, 2 columns
    for la in pot_la.values():
        if la[0] is not None:
            axes[0].plot(dates, [a*density/cf_cm**2 for a in la]) # , color=archicrop_color, alpha=0.6)
    axes[0].plot(dates, [a*density/cf_cm**2 for a in leaf_area_plant], color=stics_color)
    axes[0].set_xticks(np.arange(0, len(dates)+1, (len(dates)+1)/8))
    axes[0].set_ylabel("LAI (m²/m²)", fontsize=16, fontname="Times New Roman")

    legend_elements_lai = [
        Line2D([0], [0], color=stics_color, alpha=0.9, lw=2, label='LAI STICS'),
        Line2D([0], [0], color=archicrop_color, alpha=0.6, lw=2, label='LAI potential morphotypes')
    ]
    axes[0].legend(handles=legend_elements_lai, loc=2, prop={'family': 'Times New Roman', 'size': 12})

    for height in pot_h.values():
        if height[0] is not None:
            axes[1].plot(dates, [h/cf_cm for h in height]) #, color=archicrop_color, alpha=0.6)
    axes[1].plot(dates, [h/cf_cm for h in height_canopy], color=stics_color)
    axes[1].set_xticks(np.arange(0, len(dates)+1, (len(dates)+1)/8))
    axes[1].set_xlabel("Date", fontsize=16, fontname="Times New Roman")
    axes[1].set_ylabel("Crop height (m)", fontsize=16, fontname="Times New Roman")

    legend_elements_height = [
        Line2D([0], [0], color=stics_color, alpha=0.9, lw=2, label='Height STICS'),
        Line2D([0], [0], color=archicrop_color, alpha=0.6, lw=2, label='Height potential morphotypes')
    ]
    axes[1].legend(handles=legend_elements_height, loc=2, prop={'family': 'Times New Roman', 'size': 12})

    # Adjust layout
    plt.tight_layout()

    # Save figure
    today_str = date.today().strftime("%Y-%m-%d")
    os.makedirs(f"D:/PhD_Oriane/simulations_ArchiCrop/{today_str}", exist_ok=True)  # noqa: PTH103
    plt.savefig(f"D:/PhD_Oriane/simulations_ArchiCrop/{today_str}/plot_constrained_vs_pot.png")

    # Show the plot
    plt.show()


def plot_constrainted_vs_realized(dates, LA_archicrop, height_archicrop, leaf_area_plant, sen_leaf_area_plant, height_canopy, density, stics_color="orange", archicrop_color="green"):

    # conversion factor
    cf_cm = 100

    fig, axes = plt.subplots(2, 1, figsize=(12, 6), sharex=True)  # 1 row, 2 columns

    axes[0].plot(dates, [(la-sen)*density/cf_cm**2 for la, sen in zip(leaf_area_plant, sen_leaf_area_plant)], color=stics_color, linewidth=6)
    for result in LA_archicrop.values():
        if result[0] is not None:
            axes[0].plot(dates, [r*density/cf_cm**2 for r in result], color=archicrop_color) #, alpha=0.6)
    axes[0].set_ylabel("LAI (m²/m²)", fontsize=16, fontname="Times New Roman")
    axes[0].set_xticks(np.arange(0, len(dates)+1, (len(dates)+1)/9))
    # axes[0].set_title("Leaf Area: 3D canopy vs. STICS")
    # axes[0].legend(loc=2)

    legend_elements_lai = [
        Line2D([0], [0], color=stics_color, alpha=0.9, lw=6, label='STICS'),
        Line2D([0], [0], color=archicrop_color, alpha=0.6, lw=2, label='ArchiCrop morphotypes')
    ]
    axes[0].legend(handles=legend_elements_lai, loc=2, prop={'family': 'Times New Roman', 'size': 12})


    axes[1].plot(dates, [h/cf_cm for h in height_canopy], color=stics_color, linewidth=6)
    for result in height_archicrop.values():
        if result[0] is not None:
            axes[1].plot(dates, [r/cf_cm for r in result], color=archicrop_color) #, alpha=0.6)
    axes[1].set_xlabel("Date", fontsize=16, fontname="Times New Roman")
    axes[1].set_ylabel("Crop height (m)", fontsize=16, fontname="Times New Roman")
    axes[0].set_xticks(np.arange(0, len(dates)+1, (len(dates)+1)/9))
    # axes[1].set_title("Plant height: 3D canopy vs. STICS")

    legend_elements_height = [
        Line2D([0], [0], color=stics_color, alpha=0.9, lw=6, label='STICS'),
        Line2D([0], [0], color=archicrop_color, alpha=0.6, lw=2, label='ArchiCrop morphotypes')
    ]
    axes[1].legend(handles=legend_elements_height, loc=2, prop={'family': 'Times New Roman', 'size': 12})

    # Adjust layout
    plt.tight_layout()

    # Save figure
    today_str = date.today().strftime("%Y-%m-%d")
    os.makedirs(f"D:/PhD_Oriane/simulations_ArchiCrop/{today_str}", exist_ok=True)  # noqa: PTH103
    plt.savefig(f"D:/PhD_Oriane/simulations_ArchiCrop/{today_str}/plot_constrainted_vs_realized.png")

    # Show the plot
    plt.show()


def plot_faPAR(dates, nrj_per_plant, par_incident, par_stics, density, stics_color="orange", archicrop_color="green"):
    # curves_array = np.array(nrj_per_plant)

    # # Calculate the envelope: min and max values for each time point
    # min_values = curves_array.min(axis=0)
    # max_values = curves_array.max(axis=0)

    # Plotting the envelope along with individual curves for context
    fig, ax = plt.subplots(figsize=(12, 6))
    for curve in nrj_per_plant.values():
        # ax.plot(dates, [nrj*density/par for nrj,par in zip(curve, par_incident)]) #, color=archicrop_color, alpha=0.4, label="ArchiCrop x Caribu")
        ax.plot(dates, [nrj/par for nrj,par in zip(curve, par_incident)]) #, color=archicrop_color, alpha=0.4, label="ArchiCrop x Caribu")
        # ????????????????

    # ax.fill_between(time_points, min_values, max_values, color="skyblue", alpha=0.4)
    # ax.plot(time_points, min_values, color="blue", linestyle="--", label="Min 3D")
    # ax.plot(time_points, max_values, color="red", linestyle="--", label="Max 3D")
    ax.plot(dates, par_stics, color=stics_color, label="STICS")

    # Labels and legend
    ax.set_xticks(np.arange(0, len(dates)+1, (len(dates)+1)/8))
    ax.set_xlabel("Dates", fontsize=16, fontname="Times New Roman") 
    ax.set_ylabel("Fraction of absorbed PAR", fontsize=16, fontname="Times New Roman")
    ax.set_title("Fraction of absorbed PAR: 3D canopy vs. STICS", fontsize=16, fontname="Times New Roman")
    ax.legend()

    # Save figure
    today_str = date.today().strftime("%Y-%m-%d")
    os.makedirs(f"D:/PhD_Oriane/simulations_ArchiCrop/{today_str}", exist_ok=True)  # noqa: PTH103
    plt.savefig(f"D:/PhD_Oriane/simulations_ArchiCrop/{today_str}/plot_faPAR.png")

    plt.show()


def plot_PAR(dates, nrj_per_plant, par_incident, par_stics, density, stics_color="orange", archicrop_color="green"):
    # curves_array = np.array(nrj_per_plant)

    # # Calculate the envelope: min and max values for each time point
    # min_values = curves_array.min(axis=0)
    # max_values = curves_array.max(axis=0)

    # Plotting the envelope along with individual curves for context
    fig, ax = plt.subplots(figsize=(12, 6))
    for curve in nrj_per_plant.values():
        # ax.plot(dates, [nrj*density for nrj in curve]) #, color=archicrop_color, alpha=0.4, label="ArchiCrop x Caribu")
        ax.plot(dates, [nrj for nrj in curve]) #, color=archicrop_color, alpha=0.4, label="ArchiCrop x Caribu")
        # ????????????????

    # ax.fill_between(time_points, min_values, max_values, color="skyblue", alpha=0.4)
    # ax.plot(time_points, min_values, color="blue", linestyle="--", label="Min 3D")
    # ax.plot(time_points, max_values, color="red", linestyle="--", label="Max 3D")
    ax.plot(dates, [abs*inc for abs,inc in zip(par_stics, par_incident)], color=stics_color, label="STICS")

    # Labels and legend
    ax.set_xticks(np.arange(0, len(dates)+1, (len(dates)+1)/8))
    ax.set_xlabel("Dates", fontsize=16, fontname="Times New Roman") 
    ax.set_ylabel("Absorbed PAR", fontsize=16, fontname="Times New Roman")
    ax.set_title("Absorbed PAR: 3D canopy vs. STICS", fontsize=16, fontname="Times New Roman")
    ax.legend()

    # Save figure
    today_str = date.today().strftime("%Y-%m-%d")
    os.makedirs(f"D:/PhD_Oriane/simulations_ArchiCrop/{today_str}", exist_ok=True)  # noqa: PTH103
    plt.savefig(f"D:/PhD_Oriane/simulations_ArchiCrop/{today_str}/plot_PAR.png")

    plt.show()


def plot_extinction_coef(extinP_stics, extin_coefs, dates):
    fig, ax = plt.subplots(figsize=(12, 6))
    for curve in extin_coefs.values():
        ax.plot(dates, curve)
    ax.plot(dates, [extinP_stics]*len(dates), color="black", label="STICS")
    ax.set_xticks(np.arange(0, len(dates)+1, (len(dates)+1)/8))
    ax.set_xlabel("Dates", fontsize=16, fontname="Times New Roman") 
    ax.set_ylabel("Extinction coefficient", fontsize=16, fontname="Times New Roman")
    ax.legend()
    # Save figure
    today_str = date.today().strftime("%Y-%m-%d")
    os.makedirs(f"D:/PhD_Oriane/simulations_ArchiCrop/{today_str}", exist_ok=True)  # noqa: PTH103
    plt.savefig(f"D:/PhD_Oriane/simulations_ArchiCrop/{today_str}/plot_extin_coef.png")
    plt.show()