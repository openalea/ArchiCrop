from __future__ import annotations

archi_maize = {
    "nb_phy": 22, # number of phytomers on the main stem
    "nb_short_phy": 5,
    "short_phy_height": 3,
    
    # Stem
    "height": 250,
    "stem_q": 1, # parameter for ligule height geometric distribution along axis
    "diam_base": 5.0, # stem base diameter cm
    "diam_top": 1.5, # stem top diameter cm

    # Leaf area distribution along the stem 
    "leaf_area": 12000,
    "rmax": [0.5,0.9], # relative position of largest leaf on the stem
    "skew": [0.0005,0.1], # skewness for leaf area distribution along axis 

    # blade area
    "wl": 0.13, # leaf blade width-to-length ratio 
    "klig": 0.6, # parameter for leaf blade shape
    "swmax": 0.55, # relative position of maximal blade width
    "f1": 0.64, # parameter for leaf blade shape
    "f2": 0.92, # parameter for leaf blade shape

    # blade curvature
    "insertion_angle": 30, # leaf blade insertion angle
    "scurv": 0.8, #  relative position of inflexion point
    "curvature": 120, # leaf blade insertion-to-tip angle
    "phyllotactic_angle": 180, # phyllotactic angle
    "phyllotactic_deviation": 20, # half-deviation to phyllotactic angle

    # Development
    "phyllochron": [20,55], # phyllochron, i.e. stem element appearance rate
    # "plastochron": [30,60], # plastochron, i.e. leaf blade appearance rate
    "stem_duration": 1.6,
    "leaf_duration": 1.6,

    # Tillering
    "nb_tillers": 0, # number of tillers
    "tiller_delay": 1, # delay, as factor of phyllochron, between the appearance of a phytomer and the appearance of its tiller
    "tiller_angle": 30, # degree
    "reduction_factor": 1, # reduction factor between tillers of consecutive order
    "gravitropism_coefficient": 0.05, # 0.12
}


archi_wheat = {
    "nb_phy": 10, # number of phytomers on the main stem : [10,40] (Ndiaye et al., 2021; Lafarge and Tardieu, 2002; Clerget, 2008; Ganeme et al., 2022)
    "nb_short_phy": 5,
    "short_phy_height": 3,

    # Stem
    "height": 90,
    "stem_q": 1, # parameter for ligule height distribution along axis : [1.1] (Kaitaniemi et al., 1999) 
    "diam_base": 0.8, # stem base diameter : [2.2] (Ndiaye et al., 2021)
    "diam_top": 0.3, # stem top diameter: [1.2] (Ndiaye et al., 2021)

    # Leaf area distribution along the stem
    "leaf_area": 1500,
    "rmax": 0.83, # parameter for leaf area distribution along axis : [0.6,0.8] (Kaitaniemi et al., 1999; Welcker et al., )
    "skew": 0.001, # parameter for leaf area distribution along axis : [0.0005,0.1] (Kaitaniemi et al., 1999; Welcker et al., )
    
    # blade area 
    "wl": 0.079, # leaf blade width-to-length ratio : 
    "klig": 0.6, # parameter for leaf blade shape
    "swmax": 0.55, # parameter for leaf blade shape
    "f1": 0.64, # parameter for leaf blade shape
    "f2": 0.92, # parameter for leaf blade shape

    # Leaf blade position in space
    "insertion_angle": 40, # leaf blade insertion angle : [10,50] (Truong et al., 2015; Kaitaniemi et al., 1999)
    "scurv": 0.7, # leaf blade relative inflexion point : [0.6, 0.8] ()
    "curvature": 140, # leaf blade insertion-to-tip angle : [0,130] [45, 135] (Kaitaniemi et al., 1999)
    "phyllotactic_angle": 180, # phyllotactic angle : [180] (Davis et al., 2024)
    "phyllotactic_deviation": 90, # half-deviation to phyllotactic angle : [0,90] (Davis et al., 2024)

    # Development
    "phyllochron": 40, # [30,70], # phyllochron, i.e. stem element appearance rate : [40,65 then x1.6-2.5] (Clerget, 2008)
    # "plastochron": [40,65], # plastochron, i.e. leaf blade appearance rate : [34,46 then 80-93] (Rami Kumar et al., 2009)
    "stem_duration": 1.6,
    "leaf_duration": 1.6,

    # Tillering
    "nb_tillers": 6, # number of tillers : [0,6] (Lafarge et al., 2002)
    "tiller_angle": 5,
    "tiller_delay": 1, # delay, as factor of phyllochron, between the appearance of a phytomer and the appearance of its tiller : [] ()
    "reduction_factor": 0.8, # reduction factor between tillers of consecutive order : [0.8,1] ()
    "gravitropism_coefficient": 0.12
    }


archi_rice = {
    "nb_phy": 14, # number of phytomers on the main stem : [10,40] (Ndiaye et al., 2021; Lafarge and Tardieu, 2002; Clerget, 2008; Ganeme et al., 2022)
    "nb_short_phy": 5,
    "short_phy_height": 3,

    # Stem
    "height": 110,
    "stem_q": 1, # parameter for ligule height distribution along axis : [1.1] (Kaitaniemi et al., 1999) 
    "diam_base": 1.0, # stem base diameter : [2.2] (Ndiaye et al., 2021)
    "diam_top": 0.5, # stem top diameter: [1.2] (Ndiaye et al., 2021)

    # Leaf area distribution along the stem
    "leaf_area": 3000,
    "rmax": 0.78, # parameter for leaf area distribution along axis : [0.6,0.8] (Kaitaniemi et al., 1999; Welcker et al., )
    "skew": 0.0001, # parameter for leaf area distribution along axis : [0.0005,0.1] (Kaitaniemi et al., 1999; Welcker et al., )
    
    # blade area 
    "wl": 0.03, # leaf blade width-to-length ratio : [0.1, 0.12] ()
    "klig": 0.6, # parameter for leaf blade shape
    "swmax": 0.55, # parameter for leaf blade shape
    "f1": 0.64, # parameter for leaf blade shape
    "f2": 0.92, # parameter for leaf blade shape

    # Leaf blade position in space
    "insertion_angle": 20, # leaf blade insertion angle : [10,50] (Truong et al., 2015; Kaitaniemi et al., 1999)
    "scurv": 0.6, # leaf blade relative inflexion point : [0.6, 0.8] ()
    "curvature": 20, # leaf blade insertion-to-tip angle : [0,130] [45, 135] (Kaitaniemi et al., 1999)
    "phyllotactic_angle": 180, # phyllotactic angle : [180] (Davis et al., 2024)
    "phyllotactic_deviation": 90, # half-deviation to phyllotactic angle : [0,90] (Davis et al., 2024)

    # Development
    "phyllochron": 60, # [30,70], # phyllochron, i.e. stem element appearance rate : [40,65 then x1.6-2.5] (Clerget, 2008)
    # "plastochron": [40,65], # plastochron, i.e. leaf blade appearance rate : [34,46 then 80-93] (Rami Kumar et al., 2009)
    "stem_duration": 1.6,
    "leaf_duration": 1.6,

    # Tillering
    "nb_tillers": 19, # number of tillers : [0,6] (Lafarge et al., 2002)
    "tiller_angle": 15,
    "tiller_delay": 1, # delay, as factor of phyllochron, between the appearance of a phytomer and the appearance of its tiller : [] ()
    "reduction_factor": 0.93, # reduction factor between tillers of consecutive order : [0.8,1] ()
    "gravitropism_coefficient": 0.12
}


archi_sorghum = {
    "nb_phy": 18, # number of phytomers on the main stem : [10,40] (Ndiaye et al., 2021; Lafarge and Tardieu, 2002; Clerget, 2008; Ganeme et al., 2022)
    "nb_short_phy": 5,
    "short_phy_height": 3,

    # Stem
    "height": 200,
    "stem_q": 1, # parameter for ligule height distribution along axis : [1.1] (Kaitaniemi et al., 1999) 
    "diam_base": 4.0, # stem base diameter : [2.2] (Ndiaye et al., 2021)
    "diam_top": 1.5, # stem top diameter: [1.2] (Ndiaye et al., 2021)

    # Leaf area distribution along the stem
    "leaf_area": 15000,
    "rmax": [0.55,0.9], # parameter for leaf area distribution along axis : [0.6,0.8] (Kaitaniemi et al., 1999; Welcker et al., )
    "skew": [-10,0], # parameter for leaf area distribution along axis : [0.0005,0.1] (Kaitaniemi et al., 1999; Welcker et al., )
    
    # blade area 
    "wl": 0.12, # leaf blade width-to-length ratio : [0.1, 0.12] ()
    "klig": 0.6, # parameter for leaf blade shape
    "swmax": 0.55, # parameter for leaf blade shape
    "f1": 0.64, # parameter for leaf blade shape
    "f2": 0.92, # parameter for leaf blade shape

    # Leaf blade position in space
    "insertion_angle": 60, # leaf blade insertion angle : [10,50] (Truong et al., 2015; Kaitaniemi et al., 1999)
    "scurv": 0.75, # leaf blade relative inflexion point : [0.6, 0.8] ()
    "curvature": 90, # leaf blade insertion-to-tip angle : [0,130] [45, 135] (Kaitaniemi et al., 1999)
    "phyllotactic_angle": 180, # phyllotactic angle : [180] (Davis et al., 2024)
    "phyllotactic_deviation": 20, # half-deviation to phyllotactic angle : [0,90] (Davis et al., 2024)

    # Development
    "phyllochron": [30,80], # phyllochron, i.e. stem element appearance rate : [40,65 then x1.6-2.5] (Clerget, 2008)
    # "plastochron": [40,65], # plastochron, i.e. leaf blade appearance rate : [34,46 then 80-93] (Rami Kumar et al., 2009)
    "stem_duration": 1.6,
    "leaf_duration": 1.6,

    # Tillering
    "nb_tillers": 0,#2, # number of tillers : [0,6] (Lafarge et al., 2002)
    "tiller_angle": 20,
    "tiller_delay": 1, # delay, as factor of phyllochron, between the appearance of a phytomer and the appearance of its tiller : [] ()
    "reduction_factor": 0.9, # reduction factor between tillers of consecutive order : [0.8,1] ()
    "gravitropism_coefficient": 0.12
}


archi_sorghum_angles = {
    "nb_phy": [12,16,20,24,28], # [10,30], # number of phytomers on the main stem : [10,40] (Ndiaye et al., 2021; Lafarge and Tardieu, 2002; Clerget, 2008; Ganeme et al., 2022)
    "nb_short_phy": 4,
    "short_phy_height": 3,

    # Stem
    # "height": [50,400],
    "stem_q": 1, # parameter for ligule height distribution along axis : [1.1] (Kaitaniemi et al., 1999) 
    "diam_base": 4.0, # stem base diameter : [2.2] (Ndiaye et al., 2021)
    "diam_top": 1.5, # stem top diameter: [1.2] (Ndiaye et al., 2021)

    # Leaf area distribution along the stem
    # "leaf_area": [1000,20000],
    "rmax": [0.55,0.9], # parameter for leaf area distribution along axis : [0.6,0.8] (Kaitaniemi et al., 1999; Welcker et al., )
    "skew": [0.0005,0.1], # 0.0001 # parameter for leaf area distribution along axis : [0.0005,0.1] (Kaitaniemi et al., 1999; Welcker et al., )
    
    # blade area 
    "wl": 0.12, # leaf blade width-to-length ratio : [0.1, 0.12] ()
    "klig": 0.6, # parameter for leaf blade shape
    "swmax": 0.55, # parameter for leaf blade shape
    "f1": 0.64, # parameter for leaf blade shape
    "f2": 0.92, # parameter for leaf blade shape

    # Leaf blade position in space
    "insertion_angle": [10,20,30,40,50,60,70,80,90], # [15,30,45,60,75,90], # leaf blade insertion angle : [10,80] (Truong et al., 2015; Kaitaniemi et al., 1999)
    "scurv": 0.5, # leaf blade relative inflexion point : [0.5, 0.9] ()
    "curvature": 60, # leaf blade insertion-to-tip angle : [0, 90] (Kaitaniemi et al., 1999)
    "phyllotactic_angle": 180, # phyllotactic angle : [90;180] (Davis et al., 2024)
    "phyllotactic_deviation": 30, # half-deviation to phyllotactic angle : [0,90] (Davis et al., 2024)

    # Development
    "phyllochron": [30,80], #51 # phyllochron, i.e. stem element appearance rate : [30,70] [40,65 then x1.6-2.5] (Clerget, 2008)
    # "plastochron": 54, # plastochron, i.e. leaf blade appearance rate : [34,46 then 80-93] (Rami Kumar et al., 2009)
    "stem_duration": 1.6,
    "leaf_duration": 1.6,

    # Tillering
    "nb_tillers": 0, # number of tillers : [0,6] (Lafarge et al., 2002)
    "tiller_angle": 20,
    "tiller_delay": 1, # delay, as factor of phyllochron, between the appearance of a phytomer and the appearance of its tiller : [] ()
    "reduction_factor": 1, # reduction factor between tillers of consecutive order : [] ()

    # "inter_row": [0.4,0.6]
}


archi_sorghum_PMA = {
    "nb_phy": 16, # number of phytomers on the main stem : [10,40] (Ndiaye et al., 2021; Lafarge and Tardieu, 2002; Clerget, 2008; Ganeme et al., 2022)
    "nb_short_phy": 5,
    "short_phy_height": 3,

    # Stem
    "height": 200,
    "stem_q": 1.02, # parameter for ligule height distribution along axis : [1.1] (Kaitaniemi et al., 1999) 
    "diam_base": 4.0, # stem base diameter : [2.2] (Ndiaye et al., 2021)
    "diam_top": 1.5, # stem top diameter: [1.2] (Ndiaye et al., 2021)

    # Leaf area distribution along the stem
    "leaf_area": 4000,
    "rmax": 0.9, # parameter for leaf area distribution along axis : [0.6,0.8] (Kaitaniemi et al., 1999; Welcker et al., )
    "skew": 0.001, # parameter for leaf area distribution along axis : [0.0005,0.1] (Kaitaniemi et al., 1999; Welcker et al., )
    
    # blade area 
    "wl": 0.12, # leaf blade width-to-length ratio : [0.1, 0.12] ()
    "klig": 0.6, # parameter for leaf blade shape
    "swmax": 0.55, # parameter for leaf blade shape
    "f1": 0.64, # parameter for leaf blade shape
    "f2": 0.92, # parameter for leaf blade shape

    # Leaf blade position in space
    "insertion_angle": 60, # leaf blade insertion angle : [10,50] (Truong et al., 2015; Kaitaniemi et al., 1999)
    "scurv": 0.7, # leaf blade relative inflexion point : [0.6, 0.8] ()
    "curvature": 90, # leaf blade insertion-to-tip angle : [45, 135] (Kaitaniemi et al., 1999)
    "phyllotactic_angle": 180, # phyllotactic angle : [180] (Davis et al., 2024)
    "phyllotactic_deviation": 20, # half-deviation to phyllotactic angle : [0,90] (Davis et al., 2024)

    # Development
    "phyllochron": 51, # [30,70], # phyllochron, i.e. stem element appearance rate : [40,65 then x1.6-2.5] (Clerget, 2008)
    # "plastochron": [40,65], # plastochron, i.e. leaf blade appearance rate : [34,46 then 80-93] (Rami Kumar et al., 2009)
    "stem_duration": 1.6,
    "leaf_duration": 1.6,

    # Tillering
    "nb_tillers": 0, # number of tillers : [0,6] (Lafarge et al., 2002)
    "tiller_angle": 20,
    "tiller_delay": 1, # delay, as factor of phyllochron, between the appearance of a phytomer and the appearance of its tiller : [] ()
    "reduction_factor": 0.8 # reduction factor between tillers of consecutive order : [0.8,1] ()
}



archi_cowpea = {
    "nb_phy": 30, # number of phytomers on the main stem
    "nb_short_phy": 2,
    "short_phy_height": 3,
    
    # Stem
    "height": 100,
    "stem_q": 1, # parameter for ligule height geometric distribution along axis
    "diam_base": 1.5, # stem base diameter cm
    "diam_top": 0.5, # stem top diameter cm

    # Leaf area distribution along the stem 
    "leaf_area": 12000,
    "rmax": [0.1,1.0], # relative position of largest leaf on the stem
    "skew": [-10,0], # skewness for leaf area distribution along axis 

    # blade area
    "wl": 0.8, # leaf blade width-to-length ratio 
    "klig": 0.05, # parameter for leaf blade shape
    "swmax": 0.75, # relative position of maximal blade width
    "f1": 0.5, # parameter for leaf blade shape
    "f2": 0.9, # parameter for leaf blade shape

    # blade curvature
    "insertion_angle": 80, # leaf blade insertion angle
    "scurv": 0.5, #  relative position of inflexion point
    "curvature": 20, # leaf blade insertion-to-tip angle
    "phyllotactic_angle": 137.5, # phyllotactic angle
    "phyllotactic_deviation": 45, # half-deviation to phyllotactic angle

    # Development
    "phyllochron": [10,50], # phyllochron, i.e. stem element appearance rate
    # "plastochron": [30,60], # plastochron, i.e. leaf blade appearance rate
    "stem_duration": 3.0,
    "leaf_duration": 3.0,

    # Tillering
    "nb_tillers": 0, # number of tillers
    "tiller_delay": 1, # delay, as factor of phyllochron, between the appearance of a phytomer and the appearance of its tiller
    "tiller_angle": 30, # degree
    "reduction_factor": 1, # reduction factor between tillers of consecutive order
    "gravitropism_coefficient": 0.05, # 0.12
}