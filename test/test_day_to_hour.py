from __future__ import annotations

from openalea.archicrop.sky_sources import meteo_day
from openalea.astk.sky_irradiance import sky_irradiance
from openalea.astk.sky_sources import caribu_light_sources, sky_sources


def test_day_to_hour():

    fn = 'climsorj.meteo'
    df = meteo_day()  
    location ={
    'longitude': 3.87,
    'latitude': 45,
    'altitude': 56,
    'timezone': 'Europe/Paris'}

    for row in df.itertuples():
        irr = sky_irradiance(daydate=row.daydate, day_ghi=row.rad, **location)
        sun, sky = sky_sources(sky_type='blended', sky_irradiance=irr, scale='global')
        lights = caribu_light_sources(sun, sky)
        # then caribu with caribuscene(scene,light=lights,...)
        print(lights)  


