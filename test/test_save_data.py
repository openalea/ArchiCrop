from __future__ import annotations

import numpy as np
import pandas as pd
import xarray as xr

dates = pd.to_datetime(["2026-01-21", "2026-01-22","2026-01-23","2026-01-24","2026-01-25"])

params = {
    0: {"a" : 1, "b" : 2, "c" : 3, "d" : 4, "e" : 5},
    2: {"a" : 6, "b" : 7, "c" : 8, "d" : 9, "e" : 10},
    4: {"a" : 11, "b" : 12, "c" : 13, "d" : 14, "e" : 15}
}

ds_archi = (
    pd.DataFrame.from_dict(params, orient="index")
    .to_xarray()
    .rename({"index": "id"})
)

daily_dyn = {
    "tps_ther": [20, 21, 22, 23, 24],
    "lai": [100, 101, 102, 103, 104]
}

pot_la = {
    0: [1, 2, 3, 4, 5],
    2: [2, 3, 4, 5, 6],
    4: [3, 4, 5, 6, 7]
}

pot_h = {
    0: [10, 11, 12, 13, 14],
    2: [11, 12, 13, 14, 15],
    4: [12, 13, 14, 15, 16]
}

nrj_per_leaf = {
    0: [[10], [11], [12,12], [13,13], [14,14]],
    2: [[11], [12], [13,13], [14,14], [15,15]],
    4: [[12], [13], [14,14], [15,15], [16,16]]
}

ids = list(nrj_per_leaf.keys())
n_id = len(ids)
n_time = len(dates)
n_leaf = len(nrj_per_leaf[ids[0]][-1])  # length of inner list

arr_nrj = np.zeros((n_id, n_time, n_leaf), dtype=float)

# for k,v in nrj_per_leaf.items():
#     for t, lst in enumerate(nrj_per_leaf[k]):
#         arr_nrj[k, t, :len(lst)] = lst

for i, plant_id in enumerate(ids):
    plant = nrj_per_leaf[plant_id]
    for t, lst in enumerate(plant):
        arr_nrj[i, t, :len(lst)] = lst

df_pot_la = pd.DataFrame.from_dict(pot_la, orient="index", columns=dates)
df_pot_h = pd.DataFrame.from_dict(pot_h, orient="index", columns=dates)

ds_res = xr.Dataset(
    data_vars={
        "STICS_tps_ther": (("time"), daily_dyn["tps_ther"]),
        "STICS_lai": (("time"), daily_dyn["lai"]),
        "pot_la": (("id", "time"), df_pot_la),
        "pot_h": (("id", "time"), df_pot_h),
        "nrj_per_leaf": (("id", "time", "leaf"), arr_nrj)
    },
    coords={
        "id": df_pot_la.index,
        "time": dates,
        "leaf": np.arange(n_leaf)
    }
)

ds1 = xr.merge([ds_res, ds_archi])

ds2 = xr.merge([ds_res, ds_archi])

# print(dict(zip(ds.id.values, ds.pot_la.values)))
# print(dict(zip(ds.id.values, ds.a.values)))

def concat_ds(ds_list):
    # ds_new = xr.Dataset()
    for i,ds in enumerate(ds_list):
        if i > 0:
            ds_new = ds.assign_coords(id=ds.id + offset)
            ds_tot = xr.concat([ds_tot, ds_new], dim="id")
        else:
            ds_tot = ds
        offset = int(ds.id.max()) + 2
    return ds_tot

ds_new = concat_ds([ds1,ds2])

# print(ds_new)
