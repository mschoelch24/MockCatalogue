import time
import numpy as np
import pandas as pd
from dustmaps.config import config
config.reset()
from dustmaps.bayestar import BayestarQuery
from dustmaps.decaps import DECaPSQueryLite
import dustmaps.decaps
from pathlib import Path
from astropy.coordinates import SkyCoord
from astropy import units as u
from cat_functions import compext_parallel, extinction_calc_wrapper, extinction_lm

# Lallement+ Marshall extinction map:
# read extmap_Lallement22_Marshall06.dat to a list of lists
datContent = [i.strip().split() for i in open('./extmap_Lallement22_Marshall06.dat').readlines()]
df2 = pd.DataFrame(datContent, dtype = float)
df2.rename(columns={0: 'l', 1: 'b'}, inplace=True)
df2.rename(columns={x:y for x,y in zip(df2.columns,range(0,(len(df2.columns))))})

allow_cache_file = True
unique_hash = "replace"
cache_file = Path(f".cache/extinction_grid_{unique_hash}.parquet")

if cache_file.exists() and allow_cache_file:
    points = pd.read_parquet(cache_file)
else:
    print(
        "No cache file found for these settings, or cache file use disallowed; "
        "looking up extinction values for this grid."
    )
    dustmaps.bayestar.fetch()
    bayestar = BayestarQuery(version='bayestar2019')
    print("--------------------------------")
    print("Downloading Zucker+25 dust map.")
    print("WARNING: this will take at least 8 GB of disk space!")
    print("--------------------------------")
    dustmaps.decaps.fetch(mean_only=True)
    decaps = DECaPSQueryLite(mean_only=True)

def extinction_gz(l, b, d):
    l = np.asarray(l, dtype=float)
    b = np.asarray(b, dtype=float)
    d = np.asarray(d, dtype=float)

    gal_coords = SkyCoord(l * u.degree,
                          b * u.degree,
                          d * u.kpc,
                          frame="galactic")

    mask_north = gal_coords.icrs.dec.deg > -30
    mask_south = gal_coords.icrs.dec.deg <= -30

    print("Performing Bayestar query...")
    #Av = Rv * E(B-V), where Rv = 3.1, and Bayestar to E(B-V) multiplication factor 0.86
    extinction_bayestar = bayestar.query(gal_coords[mask_north], mode="mean") * 3.1 * 0.86

    print("Performing DECaPS query...")
    extinction_decaps = decaps.query(gal_coords[mask_south]) * 3.1

    extinction = np.empty((len(gal_coords)))
    extinction[mask_north] = extinction_bayestar
    extinction[mask_south] = extinction_decaps

    mask_nans = np.isnan(extinction)
    if np.any(mask_nans):
        print("Performing supplemental Lallement & Marshall query...")
        dummy_ids = np.arange(np.sum(mask_nans))
        extinction[mask_nans] = compext_parallel(l[mask_nans],b[mask_nans],d[mask_nans],dummy_ids, 1)[0]
    return extinction
