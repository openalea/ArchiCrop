from __future__ import annotations

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

df_pot_la = pd.DataFrame.from_dict(pot_la, orient="index", columns=dates)
df_pot_h = pd.DataFrame.from_dict(pot_h, orient="index", columns=dates)

ds = xr.Dataset(
    data_vars={
        "STICS_tps_ther": (["time"], daily_dyn["tps_ther"]),
        "STICS_lai": (["time"], daily_dyn["lai"]),
        "pot_la": (["id", "time"], df_pot_la),
        "pot_h": (["id", "time"], df_pot_h)
    },
    coords={
        "id": df_pot_la.index,
        # "id": range(len(pot_la)),
        "time": dates
    }
)

ds = xr.merge([ds, ds_archi])

# print(dict(zip(ds.id.values, ds.pot_la.values)))
# print(dict(zip(ds.id.values, ds.pot_la.values)))