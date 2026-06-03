from __future__ import annotations

import os
from datetime import date

import math
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.lines import Line2D
from sklearn.metrics import r2_score
from sklearn.linear_model import LinearRegression
from scipy.interpolate import splev, splrep

from openalea.archicrop.archicrop import ArchiCrop
from openalea.archicrop.cereal_axis import bell_shaped_dist
from openalea.archicrop.viability import resolve_organ_growth, compute_skew, generate_tillers, remove_reps
from openalea.archicrop.stics_io import stics_output



def plot_parameter_space_N_phi(plant_name, d, N, phi):

    fig = plt.figure(figsize=(8,6)) 
    ax = fig.add_subplot()

    ax.plot(N, list(phi.values()), linewidth=10, color='gray', alpha=0.3)
    ax.plot([n for n in d if d[n] != []], [phi[n] for n in d if d[n] != []], linewidth=10, color='green', alpha=0.3)

    ax.set_xlabel(r"$N$", fontsize=30)
    ax.set_ylabel(r"$\phi$", fontsize=30)
    ax.tick_params(labelsize=30)

    plt.savefig(f"../figures/{plant_name}_morphospace_N_phi.png", dpi=1000)
    plt.show()


def plot_parameter_space_N_rmax(plant_name, d, N, rmax, rmax_bounds):

    fig = plt.figure(figsize=(8,6)) 
    ax = fig.add_subplot()

    ax.fill_between(N, 
                    [min(rmax_bounds)]*len(N), 
                    [max(rmax_bounds)]*len(N), 
                    color='gray', alpha=0.3)

    ax.fill_between([n for n in d if d[n] != []], 
                    [min(rmax[n]) for n in d if d[n] != []], 
                    [min(max(rmax_bounds),max(rmax[n])) for n in d if d[n] != []], 
                    color='green', alpha=0.3)

    ax.set_xlabel(r"$N$", fontsize=30)
    ax.set_ylabel(r"$r_{max}$", fontsize=30)
    ax.tick_params(labelsize=30)

    plt.savefig(f"../figures/{plant_name}_morphospace_N_rmax.png", dpi=1000)
    plt.show()


def leaf_profiles(plant_name, N, params, crop_file_stress, crop_file_no_stress, tec_file, plant_file):

    fig, ax = plt.subplots(1, len(N), figsize=(6*len(N), 6))

    first = True

    for k, nb_phy in enumerate(N):
        mins = np.zeros(nb_phy)
        maxs = np.zeros(nb_phy)
        for pot_factor in [1.0,1.25,1.5,1.75,2]:
            for pot in [True, False]:
                if pot:
                    stics_output_file=crop_file_no_stress
                    color_stics = "cornflowerblue"
                    label_stics = "STICS unstressed"
                    color_archicrop = "green"
                    label_archicrop = "ArchiCrop unstressed"
                    label_viable = "Viable unstressed"
                else:
                    stics_output_file=crop_file_stress
                    color_stics = "orange"
                    label_stics = "STICS stressed"
                    color_archicrop = "coral"
                    label_archicrop = "ArchiCrop stressed"
                    label_viable = "Viable stressed"
                density, daily_dynamics, lifespan, lifespan_early, thermal_time, leaf_area_plant, sen_leaf_area_plant, height_canopy, end_juv, end_veg, index_end_veg = stics_output(tec_file, plant_file, stics_output_file)

                if pot:
                    ligul_factor = params['leaf_duration']
                    nb_tillers = params['nb_tillers']
                    reduction_factor = params['reduction_factor']

                    pot_leaf_area = pot_factor * max(leaf_area_plant)
                    phyllochron = (end_veg)/(nb_phy - 1 + ligul_factor)
                    ligulochron = phyllochron * ligul_factor
                
                    leaf_area_plant = leaf_area_plant[:index_end_veg+2]
                    thermal_time = thermal_time[:index_end_veg+2]
                    # end_veg = thermal_time[index_end_veg+1]
                    
                    time_unique, leaf_area_unique = remove_reps(thermal_time, leaf_area_plant)

                    # Interpolate growth dynamics
                    spl_la = splrep(time_unique, leaf_area_unique)
                 
                    starts = [i * phyllochron for i in range(nb_phy)]
                    ends = [start + phyllochron * ligul_factor for start in starts]

                    # Compute growth dynamics at given points (instead of daily)
                    la_ends = splev(ends, spl_la)

                    g, tillers = generate_tillers(nb_phy, phyllochron, ligul_factor, nb_tillers) # Dict for tillers: {id: [order, delay]}

                    # Numeric solution for minimal constrained growth (H: linear organ growth)
                    min_leaf_areas,_,_,_ = resolve_organ_growth(nb_phy, ligul_factor, la_ends, nb_tillers, tillers, reduction_factor)
                    # ax[k].plot(min_leaf_areas, range(1, nb_phy + 1), linewidth=2, color='black')
        
                    S = []
                    tt = np.arange(max(ends))
                    for t in tt:
                        sum_temp = 0
                        for i,(s,e) in enumerate(zip(starts,ends)):
                            if s <= t < e:
                                sum_temp += (t-s)/(e-s) * min_leaf_areas[i]
                            elif t >= e:
                                sum_temp += min_leaf_areas[i]
                        S.append(sum_temp)

                    
                    # Find viable rmax interval
                    for i, la in enumerate(min_leaf_areas):
                        if la == max(min_leaf_areas):
                            id_max = i+1
                            break
                    # rmax_int = np.linspace(max(0,min(max((id_max-1)/nb_phy, min(params['rmax'])), max(params['rmax']))), min(1,max(min((id_max+1)/nb_phy, max(params['rmax'])), min(params['rmax']))), 10)
                    rmax_int = np.linspace(max(min(params['rmax']),(id_max-1)/nb_phy), min(max(params['rmax']),(id_max+1)/nb_phy), 10)

                    # Find viable (rmax, skew) pairs
                    leaf_areas_norm = [la/max(min_leaf_areas) for la in min_leaf_areas]
                    
                    skews_rmax = {}
                    id = 0
                    for rank in range(1, nb_phy+1):
                        for rmax in rmax_int:
                            if rmax != rank/nb_phy and min_leaf_areas[rank-1] > 0:
                                skew = compute_skew(rank=rank, nb_phy=nb_phy, rmax=rmax, leaf_areas=leaf_areas_norm)
                                if -10 < skew < 0:
                                    skews_rmax[id] = (skew, rmax)
                                    id += 1
        
                    shapes = []
                    for id, (skew,rmax) in list(skews_rmax.items()):
                        bell_shaped_leaf_areas = bell_shaped_dist(pot_leaf_area, nb_phy, rmax, skew)
                        viable = True
                        for i,bs in enumerate(bell_shaped_leaf_areas):
                            if bs < min_leaf_areas[i]: 
                                viable = False
                                del skews_rmax[id]
                                break
                        if viable:
                            shapes.append(bell_shaped_leaf_areas)
    
                    if len(shapes) > 0:
                        shapes = np.array(shapes)
                        mins = np.array([new if new < old or old == 0 else old for old, new in zip(mins, shapes.min(axis=0))])
                        maxs = np.array([new if new > old or old == 0 else old for old, new in zip(maxs, shapes.max(axis=0))])
            
                mtgs = {}
                for id,(skew,rmax) in skews_rmax.items():
                    plant = ArchiCrop(daily_dynamics=daily_dynamics, leaf_area=pot_leaf_area, height=max(height_canopy)*pot_factor, 
                                        nb_phy=nb_phy, phyllochron=phyllochron, leaf_duration=ligul_factor, rmax=rmax, skew=skew)
                    plant.generate_potential_plant()
                    growing_plant = plant.grow_plant(rate=False)
                    growing_plant_mtg = list(growing_plant.values())
                    mtgs[id] = growing_plant_mtg
    
                for id in mtgs:
                    visible_leaf_areas = []
                    for leaf in mtgs[id][79].properties()["visible_leaf_area"]:
                        visible_leaf_areas.append(mtgs[id][79].properties()["visible_leaf_area"][leaf])
                    if first and len(visible_leaf_areas) > 0:
                        ax[k].plot(visible_leaf_areas, range(1, nb_phy + 1), linewidth=2, color=color_archicrop, label=label_archicrop)
                        first = False
                    else:
                        ax[k].plot(visible_leaf_areas, range(1, nb_phy + 1), linewidth=2, color=color_archicrop)

        if first:
            ax[k].fill_betweenx(range(1,len(mins)+1), mins, maxs, color="green", alpha=0.3, label=label_viable)
        else:
            ax[k].fill_betweenx(range(1,len(mins)+1), mins, maxs, color="green", alpha=0.3)
        
        ax[k].set_yticks(np.arange(0,nb_phy+1,4))
        ax[k].tick_params(labelsize=30)
            
        
    # ax[-1].set_xlabel("Leaf surface (cm²)", fontsize=24)
    # ax[2].set_ylabel("Leaf rank", fontsize=24)
    # ax[len(N)-1].legend(fontsize=14, loc="lower right")

    fig.tight_layout()  # otherwise the right y-label is slightly clipped
    
      
    plt.savefig(f"../figures/{plant_name}_leaf_profiles.png", dpi=1000)
    plt.show()


def plot_LAI_STICS_and_3D(plant_name, growing_plant, leaf_area_plant, sen_leaf_area_plant, density):

    realized_la = np.array([sum(la - sen 
                    for la, sen in zip(gp.properties()["visible_leaf_area"].values(), gp.properties()["senescent_area"].values())) 
                    for gp in growing_plant.values()])

    dae = range(len(realized_la))

    fig = plt.figure(figsize=(14, 7))
    plt.scatter(dae, [(la-sen)*density/10000 for la,sen in zip(leaf_area_plant, sen_leaf_area_plant)], color="orange", linewidth=6, label='STICS stressed')
    plt.plot(dae, realized_la*density/10000, color="coral", linewidth=4, label='ArchiCrop stressed')
    plt.xlabel("Days after emergence", fontsize=30)
    plt.ylabel("LAI", fontsize=30)
    plt.tick_params(labelsize=30)

    plt.savefig(f"../figures/{plant_name}_LAI_STICS_and_3D.png", dpi=1000)
    plt.show()


def plot_LAI_stressed_and_not(thermal_time, leaf_area_plant, sen_leaf_area_plant, leaf_area_plant_pot, sen_leaf_area_plant_pot, density):    
    dae = range(len(thermal_time))

    fig = plt.figure(figsize=(14, 7))
    plt.scatter(dae, [(la-sen)*density/10000 for la,sen in zip(leaf_area_plant_pot, sen_leaf_area_plant_pot)], color="cornflowerblue", linewidth=6, label='STICS unstressed')
    plt.plot(dae, [(la-sen)*density/10000 for la,sen in zip(leaf_area_plant_pot, sen_leaf_area_plant_pot)], color="green", linewidth=3, label='ArchiCrop unstressed')
    plt.scatter(dae, [(la-sen)*density/10000 for la,sen in zip(leaf_area_plant, sen_leaf_area_plant)], color="orange", linewidth=6, label='STICS stressed')
    plt.plot(dae, [(la-sen)*density/10000 for la,sen in zip(leaf_area_plant, sen_leaf_area_plant)], color="coral", linewidth=3, label='ArchiCrop stressed')
    plt.xlabel("Days after emergence", fontsize=24)
    plt.ylabel("LAI", fontsize=24)
    plt.yticks([0,1,2,3])
    plt.legend(fontsize=16)
    plt.tick_params(labelsize=15)

    plt.savefig(f"../figures/result_2_LAI_stressed_and_not.png", dpi=1000)
    plt.show()



def plot_LAI_and_pot(d, thermal_time, thermal_time_pot, leaf_area_plant, leaf_area_plant_pot, sen_leaf_area_plant, index_end_veg, index_end_veg_pot, end_veg, density, no_stress=True):

    time_unique, leaf_area_unique = remove_reps(thermal_time, leaf_area_plant)
    spl = splrep(time_unique, leaf_area_unique)

    fig, ax1 = plt.subplots(figsize=(12, 9))

    first = True
    mins = np.zeros(len(np.arange(0,thermal_time_pot[index_end_veg_pot+1])))
    maxs = np.zeros(len(np.arange(0,thermal_time_pot[index_end_veg_pot+1])))

    for n,l in d.items():
        for p,r in l:
            ligulochron = p * 1.6
            starts = []
            ends = []
            for i in range(n):
                start = i * p
                end = start + ligulochron
                starts.append(start)
                ends.append(end)
            
            la_starts = splev(starts, spl)
            la_ends = splev(ends, spl)
            
            leaf_areas,_,_,_ = resolve_organ_growth(n, 1.6, la_ends)
            leaf_areas_norm = [la/max(leaf_areas) for la in leaf_areas]
        
            skews = []
            for rank in range(1, n+1):
                if r != rank/n and leaf_areas[rank-1] > 0:
                    skew = compute_skew(rank=rank, nb_phy=n, rmax=r, leaf_areas=leaf_areas_norm)
                    if -10 < skew < 0:
                        skews.append(skew)
        
            for pot_factor in [1.0,1.25,1.5,1.75,2]:
                
                shapes = []
                for skew in skews:
                    bell_shaped_leaf_areas = bell_shaped_dist(leaf_area_plant_pot[-1] * pot_factor, n, r, skew)
                    viable = True
                    for i,bs in enumerate(bell_shaped_leaf_areas):
                        if bs < leaf_areas[i]: 
                            viable = False
                            break
                    if viable:
                        shapes.append(bell_shaped_leaf_areas)
                
                for j,shape in enumerate(shapes):
                    S = []
                    tt = np.arange(max(ends))
                    # print(len(tt))
                    for k,t in enumerate(tt):
                        sum_temp = 0
                        for i,(s,e) in enumerate(zip(starts,ends)):
                            if s <= t < e:
                                sum_temp += (t-s)/(e-s) * shape[i]
                            elif t >= e:
                                sum_temp += shape[i]
                        if sum_temp < mins[k] or mins[k] == 0:
                            mins[k] = sum_temp
                        if sum_temp > maxs[k]:
                            maxs[k] = sum_temp
                        S.append(sum_temp)
                
                    # if first:
                    #     ax1.plot(tt, S, color='cornflowerblue', linestyle='dotted', alpha=0.2, label="Viable")
                    #     first = False
                    # else:
                    #     ax1.plot(tt, S, color='cornflowerblue', linestyle='dotted', alpha=0.2)

    ax1.fill_between(np.concatenate((tt,np.array(thermal_time_pot[index_end_veg_pot+2:]))), 
                    [m*density/100**2 for m in mins]+[max(mins)*density/100**2]*len(thermal_time_pot[index_end_veg_pot+2:]), 
                    [m*density/100**2 for m in maxs]+[max(maxs)*density/100**2]*len(thermal_time_pot[index_end_veg_pot+2:]), 
                    alpha = 0.3, color='green')

    dae = range(len(thermal_time))

    if no_stress:
        ax1.scatter(dae, [la*density/100**2 for la in leaf_area_plant_pot], linewidth=8, alpha=0.8, color="cornflowerblue")
        ax1.plot(dae, [la*density/100**2 for la in leaf_area_plant_pot], linewidth=6, color="green")
    ax1.scatter(dae, [(la-sen)*density/100**2 for la,sen in zip(leaf_area_plant, sen_leaf_area_plant)], color="orange", linewidth=8, alpha=0.8)
    ax1.plot(dae, [(la-sen)*density/100**2 for la,sen in zip(leaf_area_plant, sen_leaf_area_plant)], color="coral", linewidth=6)

    ax1.set_xlabel("Days after emergence", fontsize=30)
    ax1.set_ylabel("LAI", fontsize=30)
    ax1.set_xticks([0,end_veg,max(thermal_time)], [0,index_end_veg,dae[-1]])
    # ax1.legend(fontsize=15)
    ax1.tick_params(labelsize=30)
    
    plt.savefig(f"../figures/result_2_viable_morphospace_dynamics.png", dpi=1000)
    plt.show()




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
    plt.savefig(f"../figures/plot_constrained_vs_pot.png", dpi=1000)

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
    plt.savefig(f"../figures/plot_constrainted_vs_realized.png", dpi=1000)

    # Show the plot
    plt.show()


def plot_faPAR(dates, nrj_per_plant, par_incident, par_stics, density, stics_color="orange", archicrop_color="green"):

    # Plotting the envelope along with individual curves for context
    fig, ax = plt.subplots(figsize=(12, 6))
    for curve in nrj_per_plant.values():
        ax.plot(dates, [nrj*density/par for nrj,par in zip(curve, par_incident)]) #, color=archicrop_color, alpha=0.4, label="ArchiCrop x Caribu")

    ax.plot(dates, par_stics, color=stics_color, label="STICS")

    # Labels and legend
    ax.set_xticks(np.arange(0, len(dates)+1, (len(dates)+1)/8))
    ax.set_xlabel("Dates", fontsize=16, fontname="Times New Roman") 
    ax.set_ylabel("Fraction of absorbed PAR", fontsize=16, fontname="Times New Roman")
    ax.set_title("Fraction of absorbed PAR: 3D canopy vs. STICS", fontsize=16, fontname="Times New Roman")
    ax.legend()

    # Save figure
    plt.savefig(f"../figures/plot_faPAR.png", dpi=1000)

    plt.show()


def plot_faPAR_variance(dates, nrj_per_plant, par_incident, par_stics, parameter, parameter_name, archicrop_color, stics_color):
    start = 20
    end = len(dates)- 1 # 35
    fig, ax = plt.subplots(figsize=(12, 6))

    amplitude = []
    variance = []
    diff = []
    diff_bis = []
    mean = []
    for i in range(end):
        faPAR_values = []
        d = []
        for k, curve in nrj_per_plant.items():
            faPAR = curve[i] / par_incident[i]
            faPAR_values.append(faPAR)
            d.append((par_stics[i] - faPAR) / par_stics[i])
        amp = (max(faPAR_values) - min(faPAR_values)) / par_stics[i]
        var = np.std(faPAR_values)
        m = np.mean(faPAR_values)
        d_bis = np.abs(m-par_stics[i])
        # print(f"Day {i}: std={var}, |mean-STICS|={d}, mean={m}")
        amplitude.append(amp)
        variance.append(var)
        diff.append(d)
        diff_bis.append(d_bis)
        mean.append(m)

    std_mean_diff = [d/v for d,v in zip(diff_bis,variance)]
    mean_diff = [np.mean(d) for d in diff]
    ax1 = ax.twinx()
    # ax2 = ax.twinx()
    # ax2.plot(range(end), diff, color="red", label="|mean 3D - STICS|", linestyle="dotted") 
    # ax1.plot(range(start,end), std_mean_diff[start:], color=(0.6,0.5,0), label="Standardized mean difference", linestyle="dotted") 
    ax1.plot(range(end), mean_diff, color=(0,0,0), label="Relative mean difference", linestyle="dotted") 
    # ax2.plot(range(end), variance, color="black", label="relative std 3D", linestyle="dotted")  
    # ax1.plot(range(end), [0.05 for i in range(end)], color="grey", linestyle="--", label="5% threshold")
    # ax1.plot(range(end), [d/v for v,d in zip(variance,diff)], color="grey")
    ax1.set_yticks([-1,0,1])

    # print(f'Std: [{min(rel_var)},{max(rel_var)}]')

    ax.plot(range(end), par_stics[:end], color=stics_color, linewidth=6, label="STICS")

    param_values = np.array([parameter[k] for k in nrj_per_plant])
    norm = plt.Normalize(param_values.min(), param_values.max())
    cmap = plt.get_cmap('summer')  # You can choose any matplotlib colormap

    first = True
    for k, curve in nrj_per_plant.items():
        color = cmap(norm(parameter[k]))
        label = f"{parameter_name}={parameter[k]:.2f}" if parameter[k] == param_values.min() or parameter[k] == param_values.max() else None
            
        if first:
            first = False
            ax.plot(range(end), [nrj/par for nrj, par in zip(curve[:end], par_incident[:end])], color=color, alpha=0.7)
        else:
            ax.plot(range(end), [nrj/par for nrj, par in zip(curve[:end], par_incident[:end])], color=color, alpha=0.7)


    # ax.fill_between(time_points, min_values, max_values, color="skyblue", alpha=0.4)
    # ax.plot(time_points, min_values, color="blue", linestyle="--", label="Min 3D")
    # ax.plot(time_points, max_values, color="red", linestyle="--", label="Max 3D")
    # ax.plot(range(end), par_stics[:end], color=stics_color, linewidth=3, label="STICS")

    # Add colorbar legend
    # sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    # sm.set_array([])
    # cbar = plt.colorbar(sm, ax=ax, orientation='vertical', pad=0.02)
    # cbar.set_label(r"$\theta_{leaf}$", fontsize=14)
    # cbar.set_ticks([param_values.min(), param_values.max()])
    # cbar.set_ticklabels([f"min: {param_values.min()}", f"max: {param_values.max()}"], fontsize=14)


    # Labels and legend
    # ax.set_xticks(np.arange(0, len(dates)+1, (len(dates)+1)/8))
    ax.set_yticks([0.0,0.2,0.4,0.6,0.8])
    ax.set_xlabel("Days after emergence", fontsize=24) 
    ax.set_ylabel("faPAR", fontsize=24)
    ax1.set_ylabel("Relative mean difference", fontsize=24)
    ax.tick_params(axis='both', labelsize=30)
    ax1.tick_params(axis='both', labelsize=30)


    # ax.legend(fontsize=14, loc='upper left')
    # ax1.legend(fontsize=20, loc='lower right')


    # Save figure
    plt.savefig(f"../figures/plot_faPAR_variance.png", dpi=1000)

    plt.show()


def plot_faPAR_parameter_only_archicrop(dates, nrj_per_plant, par_incident, par_stics, density, archi, parameter, parameter_name, archicrop_color, stics_color):
    end = len(dates)-1 #34
    fig, ax = plt.subplots(figsize=(12, 6))
    # for k,curve in nrj_per_plant.items():
    #     # ax.plot(dates, [nrj*density/par for nrj,par in zip(curve, par_incident)]) #, color=archicrop_color, alpha=0.4, label="ArchiCrop x Caribu")
    #     label = f"{parameter_name}={parameter[k]:.2f}" if parameter[k] == min(parameter) or parameter[k] == max(parameter) else None
    #     alpha = (parameter[k]-min(archi[parameter_name]))/(max(parameter)-min(archi[parameter_name]))
    #     ax.plot(dates, [nrj/par for nrj,par in zip(curve, par_incident)], color=archicrop_color, alpha=alpha, label=label) #, label=f"{parameter[k]}") #, color=archicrop_color, alpha=0.4, label="ArchiCrop x Caribu")

    # Normalize parameter values to [0, 1] for colormap
    param_values = np.array([parameter[k] for k in nrj_per_plant])
    norm = plt.Normalize(param_values.min(), param_values.max())
    cmap = plt.get_cmap('summer')  

    for k, curve in nrj_per_plant.items():
        color = cmap(norm(parameter[k]))
        label = f"{parameter_name}={parameter[k]:.2f}" if parameter[k] == param_values.min() or parameter[k] == param_values.max() else None
        ax.plot(range(end), [nrj/par for nrj, par in zip(curve[:end], par_incident[:end])], color=color, alpha=0.7) #, label=label)

    ax.plot(range(end), par_stics[:end], color="orange", linewidth=3, label="STICS")
    
    # Labels and legend
    # ax.set_xticks(np.arange(0, len(dates)+1, (len(dates)+1)/8))
    ax.set_xlabel("Days after emergence", fontsize=24, fontname="Times New Roman") 
    ax.set_ylabel("faPAR", fontsize=24, fontname="Times New Roman")
    # ax.set_title(f"faPAR as a function of {parameter_name}", fontsize=24, fontname="Times New Roman")

    # Add colorbar legend
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax, orientation='vertical', pad=0.02)
    # cbar.set_label(f"{parameter_name}")
    cbar.set_ticks([param_values.min(), param_values.max()])
    cbar.set_ticklabels([f"min: {param_values.min()}", f"max: {param_values.max()}"], fontsize=14)


    # ax.legend(fontsize=14)

    # Save figure
    plt.savefig(f"../figures/plot_faPAR_without_crop_{parameter_name}.png", dpi=1000)

    plt.show()



def plot_faPAR_parameter(dates, nrj_per_plant, par_incident, par_stics, density, archi, nb_phy, parameter, parameter_name, archicrop_color, stics_color):
    end = len(dates)- 1 # 35
    # fig, ax = plt.subplots(1,len(archi['nb_phy']),figsize=(5*len(archi['nb_phy']), 6))
    fig, ax = plt.subplots(1,2,figsize=(12, 6))


    # for k,curve in nrj_per_plant.items():
    #     # ax.plot(dates, [nrj*density/par for nrj,par in zip(curve, par_incident)]) #, color=archicrop_color, alpha=0.4, label="ArchiCrop x Caribu")
    #     label = f"{parameter_name}={parameter[k]:.2f}" if parameter[k] == min(parameter) or parameter[k] == max(parameter) else None
    #     alpha = (parameter[k]-min(archi[parameter_name]))/(max(parameter)-min(archi[parameter_name]))
    #     ax.plot(dates, [nrj/par for nrj,par in zip(curve, par_incident)], color=archicrop_color, alpha=alpha, label=label) #, label=f"{parameter[k]}") #, color=archicrop_color, alpha=0.4, label="ArchiCrop x Caribu")

    # Normalize parameter values to [0, 1] for colormap
    param_values = np.array([parameter[k] for k in nrj_per_plant])
    norm = plt.Normalize(param_values.min(), param_values.max())
    cmap = plt.get_cmap('summer')  # You can choose any matplotlib colormap

    # for i,n in enumerate(archi['nb_phy']):
    for i,n in enumerate([12,24]):
        ax[i].plot(range(end), par_stics[:end], color=stics_color, linewidth=6, label="STICS")
        ax[i].set_title(f"{n} leaves", fontsize=24)

        for k, curve in nrj_per_plant.items():
            if nb_phy[k] == n:
                color = cmap(norm(parameter[k]))
                label = f"{parameter_name}={parameter[k]:.2f}" if parameter[k] == param_values.min() or parameter[k] == param_values.max() else None
                ax[i].plot(range(end), [nrj/par for nrj, par in zip(curve[:end], par_incident[:end])], color=color) #, alpha=0.7) #, label=label)


    # Labels and legend
    # for i in range(1,len(archi['nb_phy'])):
    ax[1].set_yticks([])
    for i in range(2):
        ax[i].set_xlabel("Days after emergence", fontsize=24) 
        ax[i].tick_params(axis='both', labelsize=30)
    ax[0].set_ylabel("faPAR", fontsize=24)

    # # Add colorbar legend
    # sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    # sm.set_array([])
    # cbar = plt.colorbar(sm, ax=ax, orientation='vertical', pad=0.02)
    # cbar.set_label(r"$\theta_{leaf}$", fontsize=14)
    # cbar.set_ticks([param_values.min(), param_values.max()])
    # cbar.set_ticklabels([f"min: {param_values.min()}", f"max: {param_values.max()}"], fontsize=14)


    # ax[0].legend(fontsize=14)

    plt.tight_layout()

    # Save figure
    plt.savefig(f"../figures/plot_faPAR_{parameter_name}.png", dpi=1000)

    plt.show()


def plot_PAR(dates, nrj_per_plant, par_incident, par_stics, density, stics_color="orange", archicrop_color="green"):

    # Plotting the envelope along with individual curves for context
    fig, ax = plt.subplots(figsize=(12, 6))
    for curve in nrj_per_plant.values():
        ax.plot(dates, [nrj for nrj in curve]) #, color=archicrop_color, alpha=0.4, label="ArchiCrop x Caribu")

    ax.plot(dates, [abs*inc for abs,inc in zip(par_stics, par_incident)], color=stics_color, label="STICS")

    # Labels and legend
    ax.set_xticks(np.arange(0, len(dates)+1, (len(dates)+1)/8))
    ax.set_xlabel("Dates", fontsize=16, fontname="Times New Roman") 
    ax.set_ylabel("Absorbed PAR", fontsize=16, fontname="Times New Roman")
    ax.set_title("Absorbed PAR: 3D canopy vs. STICS", fontsize=16, fontname="Times New Roman")
    ax.legend()

    # Save figure
    plt.savefig(f"../figures/plot_PAR.png", dpi=1000)

    plt.show()


def cumsum_nrj_per_plant(dates, nrj_per_plant, par_incident, par_stics, parameter, parameter_name, archicrop_color, stics_color):
    end = len(dates)- 1 #34
    fig, ax = plt.subplots(figsize=(12, 6))

    cumsum_stics = np.cumsum([nrj*par for nrj, par in zip(par_stics, par_incident)])
    cumsum_nrj = {k:np.cumsum(v) for k,v in nrj_per_plant.items()}

    # ax.plot(range(len(cumsum_stics)), cumsum_stics, color=stics_color, linewidth=6, label="STICS")

    # Normalize parameter values to [0, 1] for colormap
    param_values = np.array([parameter[k] for k in nrj_per_plant])
    norm = plt.Normalize(param_values.min(), param_values.max())
    cmap = plt.get_cmap('summer')  


    for k, curve in cumsum_nrj.items():
        # if parameter[k] >= 0:
        color = cmap(norm(parameter[k]))
        label = f"{parameter_name}={parameter[k]:.2f}" if parameter[k] == param_values.min() or parameter[k] == param_values.max() else None
        ax.plot(range(len(curve)), curve, color=color, alpha=0.8) #, label=label)

    ax.plot(range(len(cumsum_stics)), cumsum_stics, color=stics_color, linewidth=6, label="STICS")

    # ax1 = ax.twinx()
    # for i in range(20,len(list(cumsum_nrj.values())[0])):
    #     par_values = [par[i] for k,par in cumsum_nrj.items()]
    #     amp = (max(par_values) - min(par_values)) / cumsum_stics[i]
    #     ax1.scatter(i, amp, s=10, color='black')
    # print(max(par_values) - min(par_values))
    # print(amp)


    # Labels and legend
    ax.set_xlabel("Days after emergence", fontsize=30) 
    ax.set_ylabel("Cum. aPAR (MJ/m²)", fontsize=30)
    ax.set_yticks([0,100,200,300,400,500])
    ax.tick_params(axis='both', labelsize=30)

    # Add colorbar legend
    # sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    # sm.set_array([])
    # cbar = plt.colorbar(sm, ax=ax, orientation='vertical', pad=0.02)
    # cbar.set_label(f"{parameter_name}", fontsize=14)
    # cbar.set_ticks([param_values.min(), param_values.max()])
    # cbar.set_ticklabels([f"min: {param_values.min()}", f"max: {param_values.max()}"], fontsize=14)


    # ax.legend(fontsize=14)

    # Save figure
    plt.savefig(f"../figures/plot_cumsum_PAR_N_{parameter_name}.png", dpi=1000)

    plt.show()


def sum_nrj_per_plant(dates, nrj_per_plant, par_incident, par_stics, parameter, parameter_name, archicrop_color, stics_color):
    end = len(dates)-1 #34
    fig, ax = plt.subplots(figsize=(6, 6))

    cumsum_stics = np.cumsum([nrj*par for nrj, par in zip(par_stics[:end], par_incident[:end])])
    cumsum_nrj = {k:np.cumsum(v) for k,v in nrj_per_plant.items()}


    list_sum = [curve[end-1] for curve in cumsum_nrj.values()]
    ax.boxplot(list_sum)
    ax.set_xticks(ticks = np.arange(1,3), labels = ['ArchiCrop x Caribu', 'STICS x Beer'])

    ax.scatter([2], [cumsum_stics[end-1]], color=stics_color, s=200)

    # Labels and legend
    ax.set_ylabel("Cumulated PAR", fontsize=24, fontname="Times New Roman")

    # Add colorbar legend
    # sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    # sm.set_array([])
    # cbar = plt.colorbar(sm, ax=ax, orientation='vertical', pad=0.02)
    # # cbar.set_label(f"{parameter_name}")
    # cbar.set_ticks([param_values.min(), param_values.max()])
    # cbar.set_ticklabels([f"min: {param_values.min()}", f"max: {param_values.max()}"], fontsize=14)


    # ax.legend(fontsize=14)

    # Save figure
    plt.savefig(f"../figures/plot_sum_PAR_{parameter_name}.png", dpi=1000)

    plt.show()


def plot_faPAR_per_leaf_parameter(lai, nrj_leaves, par_incident, par_stics, parameter, parameter_name, archicrop_color, stics_color):
    end = len(lai)-34
    fig, ax = plt.subplots(figsize=(12, 6))

    # Normalize parameter values to [0, 1] for colormap
    param_values = np.array([parameter[k] for k in nrj_leaves])
    norm = plt.Normalize(param_values.min(), param_values.max())
    cmap = plt.get_cmap('summer')  # You can choose any matplotlib colormap

    # for i,leaf in enumerate(nrj_leaves):
    #     for k, curve in leaf.items():
    #         color = cmap(norm(parameter[k]))
    #         label = f"{parameter_name}={parameter[k]:.2f}" if parameter[k] == param_values.min() or parameter[k] == param_values.max() else None
    #         ax[0].scatter([curve[40] / par_incident[40]], [i+1], color=color)
    #         ax[1].scatter([curve[end-1] / par_incident[end-1]], [i+1], color=color) #, label=label)

    cumsum_rank = {k: np.zeros(nb_phy[k]) for k in nrj_leaves}
    for k,plant in nrj_leaves.items():
        for i,leaf in enumerate(plant[end-1]):
            cumsum_rank[k][i:] += plant[end-1][i] if not math.isnan(plant[end-1][i]) else 0.0
    for k in nrj_leaves:
        color = cmap(norm(parameter[k]))
        label = f"{parameter_name}={parameter[k]:.2f}" if parameter[k] == param_values.min() or parameter[k] == param_values.max() else None
        ax.plot(cumsum_rank[k] / np.sum(cumsum_rank[k]), range(1,nb_phy[k]+1), color=color, alpha=0.3) 
        ax.scatter(cumsum_rank[k] / np.sum(cumsum_rank[k]), range(1,nb_phy[k]+1), color=color)

    # Labels and legend
    ax.set_ylabel("Leaf rank", fontsize=24, fontname="Times New Roman") 
    ax.set_xlabel("Cumulated aPAR", fontsize=24, fontname="Times New Roman")

    # Add colorbar legend
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax, orientation='vertical', pad=0.02)
    # cbar.set_label(f"{parameter_name}")
    cbar.set_ticks([param_values.min(), param_values.max()])
    cbar.set_ticklabels([f"min: {param_values.min()}", f"max: {param_values.max()}"], fontsize=14)


    # ax.legend(fontsize=14)

    # Save figure
    plt.savefig(f"../figures/plot_aPAR_per_leaf_{parameter_name}.png", dpi=1000)

    plt.show()


def plot_faPAR_lai_parameter(lai, nrj_per_plant, par_incident, par_stics, density, archi, parameter, parameter_name, archicrop_color, stics_color):
    end = len(lai)-1 #34
    fig, ax = plt.subplots(figsize=(12, 6))
    # for k,curve in nrj_per_plant.items():
    #     # ax.plot(dates, [nrj*density/par for nrj,par in zip(curve, par_incident)]) #, color=archicrop_color, alpha=0.4, label="ArchiCrop x Caribu")
    #     label = f"{parameter_name}={parameter[k]:.2f}" if parameter[k] == min(parameter) or parameter[k] == max(parameter) else None
    #     alpha = (parameter[k]-min(archi[parameter_name]))/(max(parameter)-min(archi[parameter_name]))
    #     ax.plot(dates, [nrj/par for nrj,par in zip(curve, par_incident)], color=archicrop_color, alpha=alpha, label=label) #, label=f"{parameter[k]}") #, color=archicrop_color, alpha=0.4, label="ArchiCrop x Caribu")

    # Normalize parameter values to [0, 1] for colormap
    param_values = np.array([parameter[k] for k in nrj_per_plant])
    norm = plt.Normalize(param_values.min(), param_values.max())
    cmap = plt.get_cmap('summer')  # You can choose any matplotlib colormap
    # cmap = (mpl.colors.ListedColormap(['green', 'orange', '#BF864D']))
    # bounds = [param_values.min(), param_values.min()+(param_values.max()-param_values.min())/3, param_values.min()+2*(param_values.max()-param_values.min())/3, param_values.max()]
    # norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

    par_stics = [1-math.exp(-0.7*la) for la in lai]
    ax.plot(lai[:end], par_stics[:end], color=stics_color, linewidth=6, label="STICS")

    for k, curve in nrj_per_plant.items():
        color = cmap(norm(parameter[k]))
        label = f"{parameter_name}={parameter[k]:.2f}" if parameter[k] == param_values.min() or parameter[k] == param_values.max() else None
        ax.plot(lai[:end], [nrj/par for nrj, par in zip(curve, par_incident)][:end], color=color, alpha=0.7) #, label=label)


    # ax.fill_between(time_points, min_values, max_values, color="skyblue", alpha=0.4)
    # ax.plot(time_points, min_values, color="blue", linestyle="--", label="Min 3D")
    # ax.plot(time_points, max_values, color="red", linestyle="--", label="Max 3D")

    # Labels and legend
    # ax.set_xticks(np.arange(0, len(dates)+1, (len(dates)+1)/8))
    ax.set_xlabel("LAI", fontsize=24) 
    ax.set_ylabel("faPAR", fontsize=24)
    ax.tick_params(axis='both', labelsize=30)
    # ax.set_title(f"faPAR as a function of {parameter_name}", fontsize=16, fontname="Times New Roman")

    # Add colorbar legend
    # sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    # sm.set_array([])
    # cbar = plt.colorbar(sm, ax=ax, orientation='vertical', pad=0.02)
    # # cbar.set_label(f"{parameter_name}")
    # cbar.set_ticks([param_values.min(), param_values.max()])
    # cbar.set_ticklabels([f"min: {param_values.min()}", f"max: {param_values.max()}"], fontsize=14)


    # ax.legend(fontsize=14)

    # Save figure
    plt.savefig(f"../figures/plot_faPAR_LAI_{parameter_name}.png", dpi=1000)

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
    plt.savefig(f"../figures/plot_extin_coef.png", dpi=1000)
    plt.show()


def plot_extinction_coef_dyn_parameter(extinP_stics, extin_coefs, dates, nrj_crop, archi, parameter, parameter_name, archicrop_color, stics_color):
    start = 20
    end = len(dates)-1 #34

    fig, ax = plt.subplots(figsize=(12, 6))
    # lai = [la*density/10000 for la in leaf_area_plant]
    # Normalize parameter values to [0, 1] for colormap
    param_values = np.array([parameter[k] for k in nrj_crop])
    norm = plt.Normalize(param_values.min(), param_values.max())
    cmap = plt.get_cmap('summer')

    ax.plot(range(end), [extinP_stics]*len(dates[:end]), color=stics_color, linewidth=6, label="STICS")

    # Plot each curve with color from colormap
    for k,curve in extin_coefs.items():
        color = cmap(norm(parameter[k]))
        ax.plot(range(start,end), curve[start:end], color=color, alpha=0.7)
        # label = f"{parameter_name}={parameter[k]:.2f}" if parameter[k] == min(parameter) or parameter[k] == max(parameter) else None
        # alpha = (parameter[k]-min(archi[parameter_name]))/(max(parameter)-min(archi[parameter_name]))
        # ax.plot(dates, curve, color=archicrop_color, alpha=alpha, label=label)
    
    # ax.set_xticks(np.arange(0, len(dates[:end])+1, (len(dates[:end])+1)/8))
    ax.set_xlabel("Days after emergence", fontsize=24, fontname="Times New Roman") 
    ax.set_ylabel("Extinction coefficient", fontsize=24, fontname="Times New Roman")
    # ax.set_title(f"Extinction coefficient as a function of {parameter_name}", fontsize=16, fontname="Times New Roman")

    # Add colorbar legend
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax, orientation='vertical', pad=0.02)
    # cbar.set_label(f"{parameter_name}")
    cbar.set_ticks([param_values.min(), param_values.max()])
    cbar.set_ticklabels([f"min: {param_values.min():.2f}", f"max: {param_values.max():.2f}"])


    ax.legend()

    # Save figure
    plt.savefig(f"../figures/plot_extin_coef_day_{parameter_name}.png", dpi=1000)

    plt.show()


def plot_extinction_coef_parameter(extinP_stics, extin_coefs, dates, nrj_crop, archi, parameter, parameter_name, insertion_angle, archicrop_color, stics_color):
    start = 20
    end = len(dates)-1 #34

    fig, ax = plt.subplots(figsize=(7, 6))
    # lai = [la*density/10000 for la in leaf_area_plant]
    # Normalize parameter values to [0, 1] for colormap
    param_values = np.array([insertion_angle[k] for k in nrj_crop])
    norm = plt.Normalize(param_values.min(), param_values.max())
    cmap = plt.get_cmap('summer')

    # Plot each curve with color from colormap
    if parameter_name == 'insertion_angle':
        ax.plot([param_values.min()-1,param_values.max()+1],[0.7,0.7], color='orange', linewidth=6, label='STICS')

        x = list(parameter.values())
        y = [np.mean(curve[start:end]) for curve in extin_coefs.values()]
        poly_values = np.polyfit(x, y, 2)
        print(f"k = {poly_values[0]:.6f} theta^2 + {poly_values[1]:.4f} theta + {poly_values[2]:.2f}")
        p = np.poly1d(poly_values)
        # ax.plot(np.arange(param_values.min(),param_values.max()+1), p(range(param_values.min(),param_values.max()+1)), color="black", linestyle='dashed')

        r2 = r2_score(y, p(x))
        print("R2 = ", r2)

        for k,curve in extin_coefs.items():
            color = cmap(norm(parameter[k]))
            ax.scatter(parameter[k], np.mean(curve[start:end]), color=color, s=200)
        ax.set_yticks([0.6,0.7,0.8,0.9,1.0])
        ax.set_xticks([10,30,50,70,90])
        ax.set_ylabel(r"\bar{k}", fontsize=24)
        ax.set_xlabel(r"$\theta_{leaf}$", fontsize=24)
    
    else:
        x = list(parameter.values())
        y = [np.mean(curve[start:end]) for curve in extin_coefs.values()]
        poly_values = np.polyfit(x, y, 1)
        print(f"k = {poly_values[0]:.4f} N + {poly_values[1]:.4f}")
        p = np.poly1d(poly_values)
        r2 = r2_score(y, p(x))
        print("R2 = ", r2)
        
        param_values = np.array([parameter[k] for k in nrj_crop])
        ax.plot([param_values.min()-1,param_values.max()+1],[0.7,0.7], color='orange', linewidth=6, label='STICS')
        for k,curve in extin_coefs.items():
            color = cmap(norm(insertion_angle[k]))
            ax.scatter(parameter[k], np.mean(curve[start:end]), color=color, s=200)
        ax.set_ylim([0.585,1.0])
        ax.set_yticks([])
        ax.set_xticks([12,16,20,24,28])
        ax.set_xlabel(r"$N$", fontsize=24)
 
    ax.tick_params(axis='both', labelsize=30)

    # Save figure
    plt.savefig(f"../figures/plot_extin_coef_f_{parameter_name}.png", dpi=1000)

    plt.show()


def plot_extinction_coef_lai_parameter(extinP_stics, extin_coefs, lai, nrj_crop, archi, parameter, parameter_name, insertion_angle, archicrop_color, stics_color):

    start = 20
    end = len(lai)-1 #34
    lai_trunc = lai[start:end]


    fig, ax = plt.subplots(figsize=(12, 6))
    # lai = [la*density/10000 for la in leaf_area_plant]
    # Normalize parameter values to [0, 1] for colormap
    param_values = np.array([parameter[k] for k in nrj_crop])
    norm = plt.Normalize(param_values.min(), param_values.max())
    cmap = plt.get_cmap('summer')
    # cmap = (mpl.colors.ListedColormap(['green', 'orange', '#BF864D']))
    # bounds = [param_values.min(), param_values.min()+(param_values.max()-param_values.min())/3, param_values.min()+2*(param_values.max()-param_values.min())/3, param_values.max()]
    # norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
    
    ax.plot(lai[:end], [extinP_stics]*len(lai[:end]), color=stics_color, linewidth=6, label="STICS")

    regressions = {}

    # Plot each curve with color from colormap
    for k,curve in extin_coefs.items():
        if parameter[k] in regressions:
            regressions[parameter[k]].append(curve[start:end])
        else:
            regressions[parameter[k]] = [curve[start:end]]
        color = cmap(norm(parameter[k]))
        ax.plot(lai[start:end], curve[start:end], color=color, alpha=0.7)

    for parameter_value, curves in regressions.items():
        mean_series = np.array(curves).mean(axis=0)

        model = LinearRegression().fit(np.array(lai_trunc).reshape(-1,1), mean_series)
        slope = model.coef_[0]
        intercept = model.intercept_
        # equation = f"extinction coef = {slope:.3f}*LAI + {intercept:.3f}"
        fit_line = model.predict(np.array(lai_trunc).reshape(-1,1))
       
        # plt.plot([lai_trunc[0],lai_trunc[0]], [0,1.2], alpha=0.3, color="green")
        # plt.plot(lai_trunc, mean_series, color="black", linewidth=2, label="Mean series")
        color = cmap(norm(parameter_value))
        # ax.plot(lai_trunc, fit_line, color='black', linewidth=3, linestyle="-")
        # ax.plot(lai_trunc, fit_line, color=color, linewidth=2, linestyle="--", label=f"Linear regression for {parameter_name}={parameter_value}")

    # ax.set_xticks(np.arange(0, len(dates)+1, (len(dates)+1)/8))
    ax.set_xlabel("LAI", fontsize=35) 
    ax.set_ylabel(r"$k$", fontsize=35)
    ax.tick_params(axis='both', labelsize=30)
    # ax.set_title(f"Extinction coefficient as a function of {parameter_name}", fontsize=16, fontname="Times New Roman")

    # Add colorbar legend
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax, orientation='vertical', pad=0.02)
    cbar.set_label(r"$\theta_{leaf}$", fontsize=30)
    cbar.set_ticks([param_values.min(), param_values.max()])
    cbar.set_ticklabels([f"{param_values.min()}", f"{param_values.max()}"], fontsize=30)


    # ax.legend(fontsize=14)

    # Save figure
    plt.savefig(f"../figures/plot_extin_coef_LAI_{parameter_name}.png", dpi=1000)

    plt.show()



'''
from pathlib import Path
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

csv_path = "D:/Downloads/data_Oriane.csv"
df = pd.read_csv(csv_path, sep=";")

agg = (
    df.groupby(["idPoint", "year"], as_index=False)
      .agg(rain=("rain", "sum"),
           tmoy=("tmoy", "mean"))
)

outdir = Path("../figures/plots_combines")
outdir.mkdir(exist_ok=True)

created = []
for idp, g in agg.groupby("idPoint"):
    g = g.sort_values("year")

    fig, ax1 = plt.subplots(figsize=(14, 7))

    ax1.plot(g["year"], g["rain"], linewidth=3, color="skyblue", label="Cumulative rainfall")
    ax1.set_xlabel("Year", fontsize=30)
    ax1.set_ylabel("Cumulative rainfall (mm)", fontsize=30)
    ax1.tick_params(labelsize=30)

    ax1.set_ylim([g['rain'].min()*0.85, g['rain'].max()*1.12])

    # Plot vertical line for min and max rainfall
    min_rain_year = g.loc[g["rain"].idxmin(), "year"]
    max_rain_year = g.loc[g["rain"].idxmax(), "year"]
    ax1.plot([min_rain_year]*2, [ax1.get_ylim()[0], min(g["rain"])], linewidth=3, color="powderblue", linestyle="--", label="Driest year")
    ax1.plot([max_rain_year]*2, [ax1.get_ylim()[0], max(g["rain"])], linewidth=3, color="steelblue", linestyle="--", label="Wettest year")
    # Add text annotations for min and max rainfall
    # ax1.text(min_rain_year, ax1.get_ylim()[1]*0.9, f"Min: {g['rain'].min():.1f} mm", color="lightskyblue", fontsize=20, ha="center")
    # ax1.text(max_rain_year, ax1.get_ylim()[1]*1.0, f"Max: {g['rain'].max():.1f} mm", color="deepskyblue", fontsize=20, ha="center")
    ax1.set_yticks([g['rain'].min(), g['rain'].max()], labels=[f"{g['rain'].min():.0f}", f"{g['rain'].max():.0f}"])

    ax2 = ax1.twinx()
    ax2.plot(g["year"], g["tmoy"], linewidth=3, color="salmon", label="Average temperature")
    ax2.set_ylabel("Average temperature (°C)", fontsize=30)
    ax2.set_yticks([g['tmoy'].min(), g['tmoy'].max()], labels=[f"{g['tmoy'].min():.1f}", f"{g['tmoy'].max():.1f}"])
    ax2.tick_params(labelsize=30)

    ax2.set_ylim([g['tmoy'].min()*0.98, g['tmoy'].max()*1.02])

    # Plot vertical line for min and max temperature
    min_temp_year = g.loc[g["tmoy"].idxmin(), "year"]
    max_temp_year = g.loc[g["tmoy"].idxmax(), "year"]
    ax2.plot([min_temp_year]*2, [ax2.get_ylim()[0], min(g["tmoy"])], linewidth=3, color="peachpuff", linestyle="--", label="Coldest year")
    ax2.plot([max_temp_year]*2, [ax2.get_ylim()[0], max(g["tmoy"])], linewidth=3, color="indianred", linestyle="--", label="Warmest year")
    # Add text annotations for min and max temperature
    # ax2.text(min_temp_year, ax2.get_ylim()[1]*0.94, f"Min: {g['tmoy'].min():.1f}°C", color="lightsalmon", fontsize=20, ha="center")
    # ax2.text(max_temp_year, ax2.get_ylim()[1]*1.0, f"Max: {g['tmoy'].max():.1f}°C", color="darksalmon", fontsize=20, ha="center")

    # Add legends
    lines_1, labels_1 = ax1.get_legend_handles_labels()
    lines_2, labels_2 = ax2.get_legend_handles_labels()
    ax1.legend(lines_1, labels_1, fontsize=20, loc="upper left")
    ax2.legend(lines_2, labels_2, fontsize=20, loc="upper right")

    fig.tight_layout()

    outfile = outdir / f"climat_{idp}.png"
    plt.savefig(outfile)
    plt.close(fig)
    created.append(str(outfile))

print("\n".join(created))


'''