from __future__ import annotations

import xml.etree.ElementTree as ET

import numpy as np
import pandas as pd
from pathlib import Path
from datetime import datetime

from .sky_sources import meteo_day


def read_xml_file(file_xml, params):
    """
    Parses an XML file and retrieves the values of the specified parameters.

    :param file_xml: Path to the XML file.
    :param params: List of parameter names to extract.
    :return: Dictionary with parameter names as keys and extracted values.
    """
    tree = ET.parse(file_xml)
    root = tree.getroot()
    
    result = {}
    
    # Search for all 'param' and 'colonne' elements in the XML
    for elem in root.findall(".//param") + root.findall(".//colonne"):
        param_name = elem.get("nom")  # Get the name attribute
        if param_name in params:
            result[param_name] = float(elem.text.strip()) if elem.text else None
    
    return result


def read_sti_file(file_sti, conv_unit=100, end=-1):
    """Reads a STICS mod_s*.sti output file and builds a dictionary.
    
    :param file: str, input file of STICS outputs :
        - somupvtsem : cumulative thermal time (°C.day)
        - laimax : canopy max LAI (m2/m2)
        - laisen(n) : senescent LAI (m2/m2)
        - hauteur : canopy height (m)
        - raint : PAR intercepted (actually, PAR absorbed) by canopy (MJ/m2)
        - trg(n) : global radiation (MJ/m2)
    :return: dict of dicts, for each time step, a dict of values from STICS outputs, converted to be used in ArchiCrop :
        - "Thermal time" (float): thermal time (in °C.day).
        - "Plant leaf area" (float): plant leaf area at a given thermal time (in cm²).
        - "Leaf area increment" (float): leaf area increment at a given thermal time (in cm²).
        - "Plant height" (float): plant height at a given thermal time (in cm).
        - "Height increment" (float): height increment at a given thermal time (in cm).
        - "Absorbed PAR" (float): absorbed PAR at a given thermal time (in MJ/m²)"""
    
    data_dict = {}
    non_zero_height_encountered = False

    with open(file_sti) as file:  # noqa: PTH123
        # Read the header line to get column names
        header = file.readline().strip().split(";")
        # Strip whitespace from column names
        stripped_header = [col.strip() for col in header if col != 'pla']

        # Find indices for date columns
        ian_idx = header.index("ian")
        mo_idx = header.index("mo")
        jo_idx = header.index("jo")

        # Initialize empty lists for each selected column in the dictionary
        data_dict = {col.strip(): [] for col in stripped_header}
        date_list = []

        # Read the rest of the lines (data rows)
        for line in file:
            values = line.strip().split(";")
            if 'pla' in header:
                values = values[:4] + values[5:]
            # Convert the values to floats
            row = {col.strip(): float(value) for col, value in zip(stripped_header, values)}
            if row["hauteur"] > 0.0 and row["laisen(n)"] < row["laimax"]:

                non_zero_height_encountered = True
                
                # Extract date values
                year = int(values[ian_idx])
                month = int(values[mo_idx])
                day = int(values[jo_idx])
                date_str = f"{year:04d}-{month:02d}-{day:02d}"
                date_list.append(date_str)

                for col in stripped_header:
                    data_dict[col.strip()].append(row[col.strip()])

            elif not non_zero_height_encountered:
                prev_tt = float(row["somupvtsem"])

            if non_zero_height_encountered and (row["hauteur"] <= 0.0 or row["laisen(n)"] >= row["laimax"]):
                break

    # start = 21 # 23
    # end = 80

    # Density
    density = data_dict["densite"][-1] # density = 20 plants/m2 = 0.002 plants/cm2

    # Thermal time
    thermal_time = [float(i) - prev_tt for i in data_dict["somupvtsem"]][:end]
    thermal_time_incr = [thermal_time[0]] + [thermal_time[i+1]-thermal_time[i] for i in range(len(thermal_time[1:]))]

    # Green LAI
    plant_leaf_area = [conv_unit**2*float(i)/density for i in data_dict["laimax"]][:end] # from m2/m2 to cm2/plant
    leaf_area_incr = [plant_leaf_area[0]] + [plant_leaf_area[i+1]-plant_leaf_area[i] for i in range(len(plant_leaf_area[1:]))]

    # Senescent LAI
    sen_leaf_area = [conv_unit**2*float(i)/density for i in data_dict["laisen(n)"]][:end] # from m2/m2 to cm2/plant
    sen_leaf_area_incr = [sen_leaf_area[0]] + [sen_leaf_area[i+1]-sen_leaf_area[i] for i in range(len(sen_leaf_area[1:]))]

    # Height
    height = [float(i)*conv_unit for i in data_dict["hauteur"]][:end] # from m to cm
    height_incr = [height[0]] + [height[i+1]-height[i] for i in range(len(height[1:]))]

    # Phenology
    emergence = data_dict["ilevs"][-1] - data_dict["jul"][0] # from pseudo julian day (from the beginning of the year) to day from begining of the simulation
    end_juv = data_dict["iamfs"][-1] - data_dict["jul"][0]
    max_lai = data_dict["ilaxs"][-1] - data_dict["jul"][0]

    # Incident and absorbed PAR
    par_rg_ratio = 0.95*0.48
    par_inc = [par_rg_ratio*float(rg) for rg in data_dict["trg(n)"]][:end]
    par_abs = [float(abs)/inc for abs, inc in zip(data_dict["raint"], par_inc)]#[:end] # to % of light intercepted, in MJ/m^2

    return {
        i+1: {"Date": date_list[i],
            "Thermal time": round(thermal_time[i],4),
            "Thermal time increment": round(thermal_time_incr[i],4),
            "Phenology": 'germination' if i+1 < emergence else 'juvenile' if emergence <= i+1 < end_juv else 'exponential' if end_juv <= i+1 < max_lai else 'repro',
            "Plant leaf area": round(plant_leaf_area[i],4), 
            "Leaf area increment": round(leaf_area_incr[i],4), 
            "Plant senescent leaf area": round(sen_leaf_area[i],4),
            "Senescent leaf area increment": round(sen_leaf_area_incr[i],4),
            "Plant height": round(height[i],4), 
            "Height increment": round(height_incr[i],4), 
            "Incident PAR": round(par_inc[i],4),
            "Absorbed PAR": round(par_abs[i],4)}
        for i in range(len(thermal_time))
    }, density


def read_csv_file_IC(file_csv, conv_unit=100):
    """Reads a STICS .csv output file of several USMs and builds a dictionary.
    
    :param file: str, input file of STICS outputs :
        - somupvtsem : cumulative thermal time (°C.day)
        - laimax : canopy max LAI (m2/m2)
        - laisen(n) : senescent LAI (m2/m2)
        - hauteur : canopy height (m)
        - raint : PAR intercepted (actually, PAR absorbed) by canopy (MJ/m2)
        - trg(n) : global radiation (MJ/m2)
    :return: dict of dicts, for each time step, a dict of values from STICS outputs, converted to be used in ArchiCrop :
        - "Thermal time" (float): thermal time (in °C.day).
        - "Plant leaf area" (float): plant leaf area at a given thermal time (in cm²).
        - "Leaf area increment" (float): leaf area increment at a given thermal time (in cm²).
        - "Plant height" (float): plant height at a given thermal time (in cm).
        - "Height increment" (float): height increment at a given thermal time (in cm).
        - "Absorbed PAR" (float): absorbed PAR at a given thermal time (in MJ/m²)"""
    
    df = pd.read_csv(file_csv)

    variables = [
        "Date",
        "somupvtsem",
        "laimax",
        "laisen_n",
        "hauteur",
        "raint",
        "trg_n",
        "ilevs",
        "iamfs",
        "ilaxs",
        "densite",
        "demande",
    ]

    nested_dict = {}

    for _, row in df.iterrows():

        situation = row["situation"]
        algorithm = row["algorithm"]
        plant = row["Plant"]
        row["Date"] = pd.to_datetime(row["Date"])
        # jul = [dt.date().timetuple().tm_yday for dt in row["Date"]]
        jul = row["Date"].dayofyear

        nested_dict.setdefault(situation, {}) \
                .setdefault(algorithm, {}) \
                .setdefault(plant, {}) [jul] = {
                        var: row[var]
                        for var in variables
                }

    output_dict = {}
    for a,sit in nested_dict.items():
        output_dict[a] = {}
        for b,algo in sit.items():
            output_dict[a][b] = {}
            for c,plt in algo.items():
                output_dict[a][b][c] = {}
                non_zero_height_encountered = False
                thermal_time_prev = 0.0
                leaf_area_prev = 0.0
                sen_leaf_area_prev = 0.0
                height_prev = 0.0
                k_emerg = 732
                k_harv = 732
                for k,day in plt.items():

                    if day["hauteur"] > 0.0 and day["laisen_n"] < day["laimax"]:
                        non_zero_height_encountered = True

                    elif not non_zero_height_encountered:
                        tt_emerg = float(day["somupvtsem"])
                        k_emerg = k

                    if non_zero_height_encountered and (day["hauteur"] <= 0.0 or day["laisen_n"] >= day["laimax"]):
                        k_harv = k
                        break

                    if k >= k_emerg and k < k_harv:
                        # Thermal time
                        thermal_time = float(day["somupvtsem"]) - tt_emerg 
                        thermal_time_incr = thermal_time - thermal_time_prev
                        thermal_time_prev = thermal_time

                        # Green LAI
                        plant_leaf_area = conv_unit**2*float(day["laimax"])/day["densite"] # from m2/m2 to cm2/plant
                        leaf_area_incr = plant_leaf_area - leaf_area_prev
                        leaf_area_prev = plant_leaf_area

                        # Senescent LAI
                        sen_leaf_area = conv_unit**2*float(day["laisen_n"])/day["densite"] # from m2/m2 to cm2/plant
                        sen_leaf_area_incr = sen_leaf_area - sen_leaf_area_prev
                        sen_leaf_area_prev = sen_leaf_area

                        # Height
                        height = float(day["hauteur"])*conv_unit # from m to cm
                        height_incr = height - height_prev
                        height_prev = height

                        # Phenology
                        emergence = 0 # from pseudo julian day (from the beginning of the year) to day from begining of the simulation
                        end_juv = day["iamfs"] - day["ilevs"]
                        max_lai = day["ilaxs"] - day["ilevs"]

                        # Incident and absorbed PAR
                        par_rg_ratio = 0.95*0.48
                        par_inc = par_rg_ratio*float(day["trg_n"])
                        par_abs = float(day["raint"])/par_inc # to % of light intercepted, in MJ/m^2

                        # N demand
                        n_demand = float(day["demande"])

                        # Density
                        density = day["densite"]

                        new_key = k-k_emerg+1
                        output_dict[a][b][c][new_key] = {"Date": day["Date"].date().strftime("%Y-%m-%d"),
                            "Thermal time": round(thermal_time,4),
                            "Thermal time increment": round(thermal_time_incr,4),
                            "Phenology": 'germination' if new_key < emergence else 'juvenile' if emergence <= new_key < end_juv else 'exponential' if end_juv <= new_key < max_lai else 'repro',
                            "Plant leaf area": round(plant_leaf_area,4), 
                            "Leaf area increment": round(leaf_area_incr,4), 
                            "Plant senescent leaf area": round(sen_leaf_area,4),
                            "Senescent leaf area increment": round(sen_leaf_area_incr,4),
                            "Plant height": round(height,4), 
                            "Height increment": round(height_incr,4), 
                            "Incident PAR": round(par_inc,4),
                            "Absorbed PAR": round(par_abs,4),
                            "N demand": round(n_demand,4),
                            "Density": round(density,4),
                            "Emergence": day["ilevs"]
                            }


    return output_dict


def read_doe_intercrop(file_csv):
    """Reads a .csv DOE file of several USMs and builds a dictionary."""
    
    df = pd.read_csv(file_csv)

    variables = [
        "species_principal",
        "species_secondary",
        "design",
        "row_orientation",
        "interrow_distance_principal",
        "interrow_distance_secondary",
        "n_rows_principal",
        "n_rows_secondary",
        "intrarow_distance"
    ]

    nested_dict = {}
    
    i = 1
    for _, row in df.iterrows():

        nested_dict[f"usm_{i}"] = {
                    var: row[var]
                    for var in variables
                }
        i+=1
    
    return nested_dict


def get_stics_management_params(file_tec_xml):
    """Retrieve STICS management parameters from an XML file."""
    params_tec = ['interrang']
    return read_xml_file(file_tec_xml, params_tec)

def get_stics_senescence_params(file_plt_xml):
    """Retrieve STICS senescence parameters from an XML file."""
    params_sen = ['durvieF', 'ratiodurvieI', 'coefb']
    return read_xml_file(file_plt_xml, params_sen)

def get_stics_dynamics(stics_output_file, end=-1):
    """Retrieve STICS growth and senescence dynamics from a STICS output file."""
    return read_sti_file(stics_output_file, end=end)

def get_stics_data(file_tec_xml, file_plt_xml, stics_output_file, end=-1):
    """Retrieve STICS management and senescence parameters, and growth dynamics."""
    tec_stics = get_stics_management_params(file_tec_xml)
    interrow = tec_stics['interrang']
    
    stics_output_data, density = get_stics_dynamics(stics_output_file, end=end)
    
    sen_stics = get_stics_senescence_params(file_plt_xml)
    lifespan = sen_stics['durvieF']
    lifespan_early = sen_stics['ratiodurvieI'] * lifespan
    
    return density, stics_output_data, lifespan, lifespan_early, interrow


def get_stics_data_IC(file_tec_xml, file_plt_xml, d_outputs, usm, algo, plant):
    """Retrieve STICS management and senescence parameters, and growth dynamics."""
    # tec_stics = get_stics_management_params(file_tec_xml)
    # interrow = tec_stics['interrang']
    interrow = None
    
    stics_output_data = d_outputs[usm][algo][plant]
    density = stics_output_data[len(stics_output_data)]["Density"]
    
    sen_stics = get_stics_senescence_params(file_plt_xml)
    lifespan = sen_stics['durvieF']
    lifespan_early = sen_stics['ratiodurvieI'] * lifespan
    
    return density, stics_output_data, lifespan, lifespan_early, interrow


def stics_weather_3d(filename, daily_dynamics):
    """Load the weather data from a file and filter it based on the first and last dates of plant growth."""
    df = meteo_day(filename)  # noqa: PD901

    # Get the first and last dates from daily_dynamics
    first_date = list(daily_dynamics.values())[0]["Date"]  # noqa: RUF015
    last_date = list(daily_dynamics.values())[-1]["Date"]

    # Use these dates to filter DataFrame
    return df[(df.daydate >= pd.to_datetime(first_date)) & (df.daydate <= pd.to_datetime(last_date))]

def stics_weather_3d_bis(filename, dates):
    """Load the weather data from a file and filter it based on the first and last dates of plant growth."""
    df = meteo_day(filename)  # noqa: PD901

    # Get the first and last dates from daily_dynamics
    first_date = dates[0]  # noqa: RUF015
    last_date = dates[-1]

    # Use these dates to filter DataFrame
    return df[(df.daydate >= pd.to_datetime(first_date)) & (df.daydate <= pd.to_datetime(last_date))]


def get_pheno(daily_dynamics: dict):
    thermal_time = [value["Thermal time"] for value in daily_dynamics.values()]

    for key, value in daily_dynamics.items():
        if value["Phenology"] == 'juvenile':
            next_key = key + 1
            if next_key in daily_dynamics and daily_dynamics[next_key]["Phenology"] == 'exponential':
                end_juv = thermal_time[key-1] 

        elif value["Phenology"] == 'exponential':
            next_key = key + 1
            if next_key in daily_dynamics and daily_dynamics[next_key]["Phenology"] == 'repro':
                end_veg = thermal_time[key-1] 
                index_end_veg = key - 1
                break
    
    return index_end_veg, end_veg, end_juv


def stics_output(tec_file, plant_file, stics_output_file):
    """Extract relevant data from STICS output files."""

    density, daily_dynamics, lifespan, lifespan_early, interrow = get_stics_data(
        file_tec_xml=tec_file,  # Path to the STICS management XML file
        file_plt_xml=plant_file,  # Path to the STICS plant XML file
        stics_output_file=stics_output_file  # Path to the STICS output file
    )

    thermal_time = [value["Thermal time"] for value in daily_dynamics.values()]
    leaf_area_plant = [value["Plant leaf area"] for value in daily_dynamics.values()]
    sen_leaf_area_plant = [value["Plant senescent leaf area"] for value in daily_dynamics.values()]
    height_canopy = [value["Plant height"] for value in daily_dynamics.values()]

    index_end_veg, end_veg, end_juv = get_pheno(daily_dynamics)

    return density, daily_dynamics, lifespan, lifespan_early, thermal_time, leaf_area_plant, sen_leaf_area_plant, height_canopy, end_juv, end_veg, index_end_veg


def read_climate_csv(file, sep=";"):
    return pd.read_csv(file, sep=sep)

def export_csv_tab(df, file_out=None, header=True):
    """Exporte en fichier tabulé."""
    df.to_csv(
        file_out,
        sep="\t",
        index=False,
        header=header
    )


def reorder_columns(df, new_order):
    """
    new_order = liste des colonnes dans l'ordre souhaité.
    """
    df_new = df.loc[:,new_order]
    return df_new


def remove_columns(df, columns):
    """
    columns = liste des colonnes à supprimer.
    """
    return df.drop(columns=columns, inplace=True)


def add_column(df, column_name, default_value):
    """
    Ajoute une colonne avec une valeur par défaut.
    """
    df[column_name] = default_value
    return df


def export_without_header(df, file_out, sep="\t"):
    """
    Exporte sans les noms de colonnes.
    """
    df.to_csv(
        file_out,
        sep=sep,
        index=False,
        header=False
    )


def split_climate_csv(df, folder_out):
    """
    Génère un fichier par couple (year, idPoint).
    Exemple :
        Bamako.1981
        Bamako.1982
        ...
    """
    dossier = Path(folder_out)
    dossier.mkdir(parents=True, exist_ok=True)

    for (idpoint, year), groupe in df.groupby(["idPoint", "year"]):
        file_name = dossier / f"{idpoint}.{year}"

        groupe.to_csv(
            file_name,
            sep="\t",
            index=False,
            header=False
        )


'''
df = read_climate_csv("D:/Downloads/data_Oriane.csv")
columns = ["w_date","tmoy","Tdewmin","Tdewmax","original_lat","original_lon"]
remove_columns(df, columns)
add_column(df, "Penman PET", -999.9)
add_column(df, "CO2", -999.9)
new_order = ["idPoint","year","Nmonth","NdayM","DOY","tmin","tmax","srad","Penman PET","rain","wind","Surfpress","CO2"]
df = reorder_columns(df, new_order)
folder_out = "../data/usms_STICS/v10/"
split_climate_csv(df, folder_out)
'''