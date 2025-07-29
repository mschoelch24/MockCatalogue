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
    df = pd.read_pickle(simname + '_coords.pkl',compression='zip') 

    # adding all the observables from df fragments:
    fragmt_list = sorted(glob.glob(simname + '_observ_out_pt*.pkl'))

    df2 = pd.read_pickle(fragmt_list[0],compression='zip')
    for fragmt in fragmt_list[1:]:
        dfn = pd.read_pickle(fragmt)
        df2 = pd.concat([df2, dfn], ignore_index=True)

    if len(df2) == len(df):
        df = df.merge(df2, how = 'inner', on='source_id')
    else:
        "Dataframes cannot be merged due to different lengths."
        pass

    # drawing ra,dec,parallax,pmra, and pmdec from Gaussian distribution with standard deviation of the respective uncertainties
    np.random.seed(42)
    ra = np.random.normal(loc=np.array(df['ra']), scale=np.array(df['ra_error']/3.6e6)) #converting uncertainty from mas to degrees
    dec = np.random.normal(loc=np.array(df['dec']), scale=np.array(df['dec_error']/3.6e6))
    parallax = np.random.normal(loc=np.array(df['parallax']), scale=np.array(df['plx_error'])) #parallax and parallax uncertainty in mas
    pmra = np.random.normal(loc=df['pmra'],scale=df['pmra_error'])
    pmdec = np.random.normal(loc=df['pmdec'], scale=df['pmdec_error'])
    radial_velocity = np.array(df['radial_velocity'])

    # update dataframe to include uncertainties
    df['ra'] = ra
    df['dec'] = dec
    df['parallax'] = parallax
    df['pmra'] = pmra
    df['pmdec'] = pmdec
    df['radial_velocity'] = radial_velocity

    # compute parallax error from RGB uncertainties:
    rel_uncert_samples = np.load("kde_rel_uncert_samples.npy")
    sampled_rel_uncert = np.random.choice(rel_uncert_samples, size=len(df), replace=True)
    sampled_error = sampled_rel_uncert * df['parallax']
    df['plx_error_RGB'] = sampled_error
    
    plx_samples = np.random.normal(loc=df['parallax'],scale=sampled_error)
    df['parallax_RGB'] = plx_samples

    # compute distance using Weiler+25
    df['R']= 1. / (df['parallax'] + df['plx_error'] * Weiler_C(df['parallax']/df['plx_error'],0.5) )
    df['R_RGB']= 1. / (df['parallax_RGB'] + df['plx_error_RGB'] * Weiler_C(df['parallax_RGB']/df['plx_error_RGB'],0.5) )
    distance = np.array(df['R_RGB'])
    
    # converting back to cartesian coordinates
    x, y, z, vx, vy, vz = equatorial2cartesian(ra, dec, distance, pmra, pmdec, radial_velocity)

    df['X'] = x
    df['Y'] = y
    df['Z'] = z
    df['Vx'] = vx
    df['Vy'] = vy
    df['Vz'] = vz

    df.to_pickle(simname + '_mock.pkl',compression='zip')
    print("Mock dataframe contains columns", list(df), "and has length", len(df))

    tf = time.time()
    t_total = tf - t0
    print("Total time compile.py:", '%.2f' % (t_total/60) ,"min")

if __name__ == "__main__":
    main()


