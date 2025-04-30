import pandas as pd
import numpy as np
import glob
import sys
import time
import os
from cat_functions import *

def main():
    input_file = sys.argv[1]
    simname,_ = os.path.splitext(input_file)
    
    t0 = time.time()
    
    # coordinate.py output dataframe:
    df1 = pd.read_pickle(simname + '_coords.pkl',compression='zip') 

    # adding all the observables from df fragments:
    fragmt_list = sorted(glob.glob('*observ_out_pt*.pkl'))

    df2 = pd.read_pickle(fragmt_list[0],compression='zip')
    for fragmt in fragmt_list[1:]:
        dfn = pd.read_pickle(fragmt)
        df2 = pd.concat([df2, dfn], ignore_index=True)

    if len(df2) == len(df1):
        df1 = df1.merge(df2, how = 'inner', on='source_id')
    else:
        "Dataframes cannot be merged due to different lengths."
        pass

    df1.to_pickle(simname + '_out.pkl',compression='zip')
    print("Final dataframe contains columns", list(df1), "and has length", len(df1))

    # drawing ra,dec,parallax,pmra, and pmdec from Gaussian distribution with standard deviation of the respective uncertainties
    np.random.seed(42)
    ra = np.random.normal(loc=np.array(df1['ra']), scale=np.array(df1['ra_error']/3.6e6)) #converting uncertainty from mas to degrees
    dec = np.random.normal(loc=np.array(df1['dec']), scale=np.array(df1['dec_error']/3.6e6))
    parallax = np.random.normal(loc=np.array(df1['parallax']), scale=np.array(df1['plx_error'])) #parallax and parallax uncertainty in mas
    pmra = np.random.normal(loc=df1['pmra'],scale=df1['pmra_error'])
    pmdec = np.random.normal(loc=df1['pmdec'], scale=df1['pmdec_error'])
    radial_velocity = np.array(df1['radial_velocity'])

    # removing unphysical negative distances
    #parallax = np.where(parallax >= 0, parallax, np.nan) 
    
    # converting back to cartesian coordinates
    x, y, z, vx, vy, vz = equatorial2cartesian(ra, dec, 1/parallax, pmra, pmdec, radial_velocity)

    dferrors = pd.DataFrame()
    dferrors['x'] = x
    dferrors['y'] = y
    dferrors['z'] = z
    dferrors['vx'] = vx
    dferrors['vy'] = vy
    dferrors['vz'] = vz
    dferrors.to_pickle(simname + '_mock.pkl',compression='zip')
    print("Mock dataframe contains columns", list(dferrors), "and has length", len(dferrors))

    tf = time.time()
    t_total = tf - t0
    print("Total time compile.py:", '%.2f' % (t_total/60) ,"min")

if __name__ == "__main__":
    main()


