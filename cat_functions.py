import os.path
import pandas as pd
from astropy.io import fits
import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord
import astropy.coordinates as coord
from multiprocessing import Pool
import gc
from scipy.interpolate import RegularGridInterpolator
from pygaia.errors.astrometric import parallax_uncertainty, proper_motion_uncertainty
from scipy.interpolate import interp1d


def open_fits(file):
    with fits.open(file) as data:
        return pd.DataFrame(data[1].data)
def open_pickle(file):
    return pd.read_pickle(file, compression='zip')

READER_MAP = {
    '.csv': pd.read_csv,
    '.pkl': open_pickle,
    '.fits': open_fits  # add more?
}

def read_sim(file):
    filename, ext = os.path.splitext(file)
    try:
        reader = READER_MAP[ext]
    except KeyError:
        raise ValueError(f'Unsupported filetype: {ext}')
    return reader(file)

def equat_heliocen(x,y,z,vx,vy,vz):
  gc = SkyCoord(x*u.kpc, y*u.kpc, z*u.kpc, v_x = vx*(u.km/u.s), v_y = vy*(u.km/u.s), v_z = vz*(u.km/u.s), frame=coord.Galactocentric, z_sun = 2 *u.pc)
  hc = gc.transform_to(coord.ICRS)
  return hc.ra.degree, hc.dec.degree, hc.distance.kpc, hc.pm_ra_cosdec.value, hc.pm_dec.value, hc.radial_velocity.value #*u.s/u.km

def gal_heliocen(x,y,z,vx,vy,vz):
  gc = SkyCoord(x*u.kpc, y*u.kpc, z*u.kpc, v_x = vx*(u.km/u.s), v_y = vy*(u.km/u.s), v_z = vz*(u.km/u.s), frame=coord.Galactocentric, z_sun = 2 *u.pc)
  hc = gc.transform_to('galactic')
  return hc.l.degree, hc.b.degree, hc.distance.kpc, hc.pm_l_cosb.value, hc.pm_b.value, hc.radial_velocity.value #*u.s/u.km

def cart_galactocen(ra, dec, distance, pmra, pmdec, vr):
  hc = coord.SkyCoord(ra*u.degree, dec*u.degree, distance*u.kpc, pm_ra_cosdec = pmra *u.mas/u.yr, pm_dec = pmdec *u.mas/u.yr, radial_velocity = vr*u.km/u.s, frame='icrs')
  gc = hc.transform_to(coord.Galactocentric) #(galcen_distance=1*u.kpc))
  return gc.x.value, gc.y.value, gc.z.value, gc.v_x.value, gc.v_y.value, gc.v_z.value

def rotationz(x,y,z,vx,vy,vz, theta = 0):
    """
    Rotating the galaxy frame positions and velocities around the z axis by angle theta (in degrees).
    """
    x2 = x * np.cos(np.radians(theta)) - y * np.sin(np.radians(theta))
    y2 = x * np.sin(np.radians(theta)) + y * np.cos(np.radians(theta))
    z2 = z

    vx2 = vx * np.cos(np.radians(theta)) - vy * np.sin(np.radians(theta))
    vy2 = vx * np.sin(np.radians(theta)) + vy * np.cos(np.radians(theta))
    vz2 = vz
    return x2, y2, z2, vx2, vy2, vz2


# Lallement+ Marshall extinction map:
# read extmap_Lallement22_Marshall06.dat to a list of lists
datContent = [i.strip().split() for i in open('./extmap_Lallement22_Marshall06.dat').readlines()]
df2 = pd.DataFrame(datContent, dtype = float)
df2.rename(columns={0: 'l', 1: 'b'}, inplace=True)
df2.rename(columns={x:y for x,y in zip(df2.columns,range(0,(len(df2.columns))))})

def extinction_calc(l, b, d, id_n):
    try:
        # index where source l, b match grid latitude and longitude
        l_index = np.floor(l)
        b_index = np.floor(b)
        l_index_up = l_index +1
        b_index_up = b_index +1

        if (l_index > 360.0):
            l_index_up = 0.0
        if (l_index < 0.0):
            l_index_up = 360.0
        if (abs(b_index) >= 89.0):
            b_index = b_index_up

        l_lim_lo = np.array(df2.index[(df2['l']== int(l_index))])
        l_lim_hi = np.array(df2.index[(df2['l']== int(l_index_up))])

        b_lim_lo = np.array(df2.index[(df2['b']== int(b_index))])
        b_lim_hi = np.array(df2.index[(df2['b']== int(b_index_up))])
        
        row_index1 = np.intersect1d(l_lim_lo, b_lim_lo)[0]
        row_index2 = np.intersect1d(l_lim_lo, b_lim_hi)[0]
        row_index3 = np.intersect1d(l_lim_hi, b_lim_lo)[0]
        row_index4 = np.intersect1d(l_lim_hi, b_lim_hi)[0]

        # make difference array from distance row - distance of source, excluding extinction values (odd indices):
        ext_row = np.asarray(df2.iloc[[row_index1]])[0]
        column_indices = np.linspace(0, len(ext_row)-1, len(ext_row))
        ext_row[~(column_indices %2==0)] = np.nan # making all extinction values nans
        ext_row[0]= np.nan
        difference_array = np.absolute(ext_row - d)

        # find column index = index of minimum element of difference array... NOT NAN!!
        col_index = column_indices[np.nanargmin(difference_array)] # -2 to reproduce the java code
        if col_index >= 266:
            col_index = 264
        col_index2 = col_index +2

        if ext_row[int(col_index)] > d:
                col_index = col_index - 2
                col_index2 = col_index2 -2

        A000 = df2.iloc[int(row_index1), int(col_index +1)]
        A001 = df2.iloc[int(row_index1), int(col_index2 +1)]
        A010 = df2.iloc[int(row_index2), int(col_index +1)]
        A011 = df2.iloc[int(row_index2), int(col_index2 +1)]
        A100 = df2.iloc[int(row_index3), int(col_index +1)]
        A101 = df2.iloc[int(row_index3), int(col_index2 +1)]
        A110 = df2.iloc[int(row_index4), int(col_index +1)]
        A111 = df2.iloc[int(row_index4), int(col_index2 +1)]

        l_vals = np.array([l_index, l_index+1])
        b_vals = np.array([b_index, b_index+1])
        d0 = df2.iloc[int(row_index1),int(col_index)]
        d1 = df2.iloc[int(row_index1),int(col_index2)]
        d_vals = np.array([d0,d1])

        extinction_vals = np.array([[[A000,A001],[A010,A011]],[[A100,A101],[A110,A111]]])
        rgi = RegularGridInterpolator((l_vals, b_vals, d_vals), extinction_vals, bounds_error=False, fill_value = None)
        ext_d = rgi(np.array([l, b, d]))[0]
    except Exception as e:
        print("Failed on star with id ", id_n)
        print(e)
        ext_d = np.nan
    return ext_d, id_n

def compext_parallel(l_input, b_input, distance, sid, ncores):
    inputs = list(zip(l_input, b_input, distance, sid))
    with Pool(processes=ncores) as pool:
        results = [pool.apply_async(extinction_calc, i) for i in inputs]
        results = [r.get() for r in results]
        sorted_results = sorted(results, key=lambda a_entry: a_entry[1])
        sorted_results = [item[0] for item in sorted_results]
    return sorted_results

def magnitude(d, Av):
    """
    Calculating the G magnitude of stars from apparent K magnitude and color.
    Input:
        d - heliocentric distance of the star in parsec
        Av - extinction
    Output: 
        G - in mag
    """
    K = -1.62 + 5 * np.log10(d*10e3) - 5 + (0.114* Av)
    color = (0.282 - 0.114) * Av + 0.55
    G = K - 0.286 + 4.023 * color - 0.35 * color **2 + 0.021 * color ** 3
    return G

def uncertainties(G, rls = 'dr3'):
    """
    Importing uncertainty factors from the PyGaia package (https://github.com/agabrown/PyGaia)
    Input:
        G - G magnitude of star (array or scalar)
    Optional: 
        rls - Gaia data release (default: 'dr3')
    Output: 
        Uncertainties (in micro-arcsec) in parallax, ra, dec, proper motion ra and proper motion dec
    
    *Note that GaiaNIR requires input file (e.g. 'GmagSig_K5IIIAv0_M5.csv', where M5 stands for the medium mission over a 5-yr baseline).
    Parallax uncertainties are computed from a RC spectrum, using GaiaNIR simulation of Hobbs et al. (in prep). *
    """
    _t_factor = {"dr3": 1.0, "dr4": 0.749, "dr5": 0.527}
    pos_alpha_factor = {"dr3": 0.8, "dr4": 0.8, "dr5": 0.8, "NIR":0.8}
    pos_delta_factor = {"dr3": 0.7, "dr4": 0.7, "dr5": 0.7, "NIR":0.7}
    pm_alpha_factor = {"dr3": 1.03, "dr4": 0.58, "dr5": 0.29, "NIR":0.29}
    pm_delta_factor = {"dr3": 0.89, "dr4": 0.50, "dr5": 0.25, "NIR":0.25}
    gatefloor = np.power(10.0, 0.4 * (13.0 - 15.0))
    z = np.maximum(gatefloor, np.power(10.0, 0.4 * (G - 15.0)))
    if (rls == 'NIR'):
        nir = pd.read_csv('GmagSig_K5IIIAv0_M5.csv')
        interp_func = interp1d(nir['Gmag'], nir['sigma'], kind='linear', fill_value="extrapolate")
        plx_unc = interp_func(G) #in micro-as
    else:
        plx_unc = np.sqrt(40 + 800 * z + 30 * z * z) * _t_factor[rls]

    ra_unc = plx_unc * pos_alpha_factor[rls]
    dec_unc = plx_unc * pos_delta_factor[rls]
    pmra_unc = plx_unc * pm_alpha_factor[rls]
    pmdec_unc = plx_unc * pm_delta_factor[rls]
    return plx_unc, ra_unc, dec_unc, pmra_unc, pmdec_unc
