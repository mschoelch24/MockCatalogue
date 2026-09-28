import pandas as pd
import numpy as np
import glob
import sys
import time
import os
from cat_functions import *

def main():
    input_file = sys.argv[1]
    tracer = sys.argv[2]
    rls = sys.argv[3]
    oldrls = sys.argv[4]
    simname,_ = os.path.splitext(input_file)
    
    t0 = time.time()
    
    # existing mock dataframe:
    df = pd.read_pickle(simname + '_mock_' + tracer + '_' + oldrls + '.pkl',compression='zip') 
    print(list(df), len(df))

    # computing uncertainties for additional release:
    print("****Starting uncertainties calculation****")
    plx_error, ra_error, dec_error, pmra_error, pmdec_error, radial_velocity_error = uncertainties(df['G'], df['bp_rp'], rls)
    df['plx_error'] = plx_error/1000 #converting to mas
    df['ra_error'] = ra_error/1000
    df['dec_error'] = dec_error/1000
    df['pmra_error'] = pmra_error/1000 #mas/yr
    df['pmdec_error'] = pmdec_error/1000
    df['radial_velocity_error'] = radial_velocity_error #in km/s

    # drawing ra, dec, parallax, pmra, and pmdec from Gaussian distribution with standard deviation of the respective uncertainties
    np.random.seed(42)
    ra = np.random.normal(loc=np.array(df['ra']), scale=np.array(df['ra_error']/3.6e6)) #converting uncertainty from mas to degrees
    dec = np.random.normal(loc=np.array(df['dec']), scale=np.array(df['dec_error']/3.6e6))
    parallax = np.random.normal(loc=np.array(df['parallax']), scale=np.array(df['plx_error'])) #parallax and parallax uncertainty in mas
    pmra = np.random.normal(loc=df['pmra'],scale=df['pmra_error'])
    pmdec = np.random.normal(loc=df['pmdec'], scale=df['pmdec_error'])
    radial_velocity = np.random.normal(loc=df['radial_velocity'], scale=df['radial_velocity_error'])
    
    # compute parallax error from RGB uncertainties:
    rel_uncert_samples = np.load("kde_rel_uncert_samples.npy")
    sampled_rel_uncert = np.random.choice(rel_uncert_samples, size=len(df), replace=True)
    sampled_error = sampled_rel_uncert * df['parallax'] 
    #this parallax is crazy large. Need to preserve unedited plx in previous files to have something more reasonable!
    df['Plx_error_RGB'] = sampled_error
    plx_samples = np.random.normal(loc=np.array(df['parallax']), scale=sampled_error)
    df['Parallax_RGB'] = plx_samples

    # update dataframe to include uncertainties
    df['Ra'] = ra
    df['Dec'] = dec
    df['Parallax'] = parallax
    df['Pmra'] = pmra
    df['Pmdec'] = pmdec
    df['Radial_velocity'] = radial_velocity

    # compute distance using Weiler+25
    df['R']= 1. / (df['parallax'] + df['plx_error'] * Weiler_C(df['parallax']/df['plx_error'],0.5) )
    df['R_RGB']= 1. / (df['parallax_RGB'] + df['plx_error_RGB'] * Weiler_C(df['parallax_RGB']/df['plx_error_RGB'],0.5) )
    distance = np.array(df['R_RGB'])

    # converting back to cartesian coordinates
    x, y, z, vx, vy, vz = equatorial2cartesian(ra, dec, distance, pmra, pmdec, np.array(radial_velocity))

    df['X'] = x
    df['Y'] = y
    df['Z'] = z
    df['Vx'] = vx
    df['Vy'] = vy
    df['Vz'] = vz

    df.to_pickle(simname+ '_mock_'+ tracer+ '_'+ rls+'.pkl',compression='zip')
    print("Mock dataframe contains columns", list(df), "and has length", len(df))
    df.to_csv(simname+ '_mock_'+ tracer+ '_'+ rls+'.csv')
 
    tf = time.time()
    t_total = tf - t0
    print("Total time compile.py:", '%.2f' % (t_total/60) ,"min")

if __name__ == "__main__":
    main()