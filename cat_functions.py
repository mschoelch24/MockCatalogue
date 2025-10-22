import os.path
import pandas as pd
from astropy.io import fits
import numpy as np
import math
import astropy.units as u
from astropy.coordinates import SkyCoord
import astropy.coordinates as coord
from multiprocessing import Pool, Manager, Value, Lock
import gc
from pygaia.errors.astrometric import parallax_uncertainty, proper_motion_uncertainty
from scipy.interpolate import RegularGridInterpolator
from scipy.interpolate import interp1d
from functools import partial
import time
import scipy.special as scp

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

# Coordinate transformations:
def cartesian2equatorial(x,y,z,vx,vy,vz):
  gc = SkyCoord(x*u.kpc, y*u.kpc, z*u.kpc, v_x = vx*(u.km/u.s), v_y = vy*(u.km/u.s), v_z = vz*(u.km/u.s), frame=coord.Galactocentric, z_sun = 2 *u.pc)
  hc = gc.transform_to(coord.ICRS)
  return hc.ra.degree, hc.dec.degree, hc.distance.kpc, hc.pm_ra_cosdec.value, hc.pm_dec.value, hc.radial_velocity.value #*u.mas/u.yr

def cartesian2galactic(x,y,z,vx,vy,vz):
  gc = SkyCoord(x*u.kpc, y*u.kpc, z*u.kpc, v_x = vx*(u.km/u.s), v_y = vy*(u.km/u.s), v_z = vz*(u.km/u.s), frame=coord.Galactocentric, z_sun = 2 *u.pc)
  hc = gc.transform_to('galactic')
  return hc.l.degree, hc.b.degree, hc.distance.kpc, hc.pm_l_cosb.value, hc.pm_b.value, hc.radial_velocity.value #*u.s/u.km

def equatorial2cartesian(ra, dec, distance, pmra, pmdec, vr):
  mask = (dec > -90) & (dec < 90) & (ra >= 0) & (ra < 360)
  hc = coord.SkyCoord(ra[mask]*u.degree, dec[mask]*u.degree, distance[mask]*u.kpc,
                      pm_ra_cosdec = pmra[mask] *u.mas/u.yr, pm_dec = pmdec[mask] *u.mas/u.yr, radial_velocity = vr[mask]*u.km/u.s, frame='icrs')
  gc = hc.transform_to(coord.Galactocentric) #(galcen_distance=1*u.kpc))
  components = (gc.x, gc.y, gc.z, gc.v_x, gc.v_y, gc.v_z)
  output = []
  for comp in components:
    arr = np.full_like(ra, np.nan)
    arr[mask] = comp.value
    output.append(arr)
  return tuple(output)

def equatorial2galactic(ra, dec, distance, pmra, pmdec, vr):
  hc = coord.SkyCoord(ra*u.degree, dec*u.degree, distance*u.kpc, pm_ra_cosdec = pmra *u.mas/u.yr, pm_dec = pmdec *u.mas/u.yr, radial_velocity = vr*u.km/u.s, frame='icrs')
  gh = hc.transform_to('galactic') #(galcen_distance=1*u.kpc))
  return gh.l.degree, gh.b.degree, gh.distance.kpc, gh.pm_l_cosb.value, gh.pm_b.value, gh.radial_velocity.value

def galactic2equatorial(l, b, distance, pm_l_cosb, pm_b, vr):
    gh = coord.SkyCoord(l*u.degree, b*u.degree, distance*u.kpc, pm_l_cosb=pm_l_cosb*u.mas/u.yr, pm_b=pm_b*u.mas/u.yr, radial_velocity=vr*u.km/u.s, frame='galactic')
    hc = gh.transform_to('icrs')
    return hc.ra.degree, hc.dec.degree, hc.distance.kpc, hc.pm_ra_cosdec.value, hc.pm_dec.value, hc.radial_velocity.value

def galactic2cartesian(x_l, x_b, x_distance, x_pm_l_cosb, x_pm_b, x_vr):
    hc = SkyCoord(l=x_l*u.degree, b=x_b*u.degree, distance=x_distance*u.kpc, pm_l_cosb=x_pm_l_cosb*u.mas/u.yr, pm_b=x_pm_b*u.mas/u.yr, radial_velocity=x_vr*u.km/u.s, frame='galactic')
    gc = hc.transform_to(coord.Galactocentric(z_sun=2*u.pc))  # Reverse transformation
    return gc.x.to(u.kpc).value, gc.y.to(u.kpc).value, gc.z.to(u.kpc).value, gc.v_x.to(u.km/u.s).value, gc.v_y.to(u.km/u.s).value, gc.v_z.to(u.km/u.s).value

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

        # find column index = index of minimum element of difference array...
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
        
        # extrapolation inflates extinction - we set boundary case:
        d_vals_lim = np.array([14.9, 15])
        extinction_vals_lim =np.array([[[A001,A001],[A011,A011]],[[A101,A101],[A111,A111]]])
        if d > 15:
            rgi = RegularGridInterpolator((l_vals, b_vals, d_vals_lim), extinction_vals_lim, bounds_error=False, fill_value = None)
        ext_d = rgi(np.array([l, b, d]))[0]
    except Exception as e:
        print("Failed on star with id ", id_n)
        print(e)
        ext_d = np.nan
    return ext_d, id_n

def extinction_calc_wrapper(args, counter, lock, start_time, total):
    result = extinction_calc(*args)
    with lock:
        counter.value += 1
        processed = counter.value
        elapsed = time.time() - start_time.value
        elapsed_str = f"{elapsed / 60:.1f} min" if elapsed >= 60 else f"{elapsed:.1f} sec"
        if processed in {1000, 10000} or processed % 100000 == 0 or processed == total:
            if processed < 300000:
                eta_str = "calculating..."
            else:
                avg_time = elapsed / processed
                eta = avg_time * (total - processed)
                eta_str = f"{eta / 60:.1f} min" if eta >= 60 else f"{eta:.1f} sec"
            print(f"{processed}/{total} stars processed - Elapsed time: {elapsed_str}, Estimated time to finish: {eta_str}")
    return result

def compext_parallel(l_input, b_input, distance, sid, ncores):
    inputs = list(zip(l_input, b_input, distance, sid))
    total = len(inputs)
    with Manager() as manager:
        counter = manager.Value('i', 0)
        lock = manager.Lock()
        start_time = manager.Value('d', time.time())
        
        with Pool(processes=ncores) as pool:
            func = partial(extinction_calc_wrapper, counter=counter, lock=lock,
                           start_time=start_time, total=total)
            results = pool.map(func, inputs)
        sorted_results = sorted(results, key=lambda a_entry: a_entry[1])
        sorted_results = [item[0] for item in sorted_results]
    return sorted_results

def magnitude(d, Av):
    """
    Calculating the G magnitude of stars from apparent K magnitude and color. 
    We assign Red Clump stellar parameters, i.e. absolute magnitude Mk = -1.62 and intrinsic colour (J-K)0 = 0.55, and use the extinction relations Ak = 0.114*Av and Aj = 0.282*Av.
    Input:
        d - heliocentric distance of the star in parsec
        Av - extinction in V band
    Output: 
        G - in mag
    """
    K = -1.62 + 5 * np.log10(d*1e3) - 5 + (0.114* Av)
    color = (0.282 - 0.114) * Av + 0.55
    G = K - 0.286 + 4.023 * color - 0.35 * color **2 + 0.021 * color ** 3
    return G

def magnitude_RGB(d, Av):
    """
    Calculating the G magnitude from RGB stellar parameters, i.e. drawing observed color and absolute magnitude from RGB distribution ('kde_color_mag_samples.npz') and converting to apparent G magnitude using the star's distance and the extinction conversion Ag/Av (Ramos et al. 2020 and J.M.Carassco, private communication).
    Input: 
        d - heliocentric distance of the star in parsec
        Av - extinction in V band
    Output: 
        G - in mag
    """
    color_mag_samples = np.load("kde_color_mag_samples.npz")
    sampled_colors = np.random.choice(color_mag_samples["color"], size=len(d), replace=True)
    sampled_magnitudes = np.random.choice(color_mag_samples["magnitude"], size=len(d), replace=True)
    
    Ag = 0.98703 * Av - 0.14452 * sampled_colors * Av + 0.01126 * sampled_colors**2 * Av + 0.03313 * Av**2 - 0.00389 * sampled_colors * Av**2
    G = sampled_magnitudes + 5 * np.log10(d*1e3) - 5 + Ag
    return G

def uncertainties(G, rls = 'dr3'):
    """
    Importing uncertainty factors from PyGaia (https://github.com/agabrown/PyGaia)
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
    if (rls[:3] == 'NIR'):
        nir = pd.read_csv('GmagSig_K5IIIAv0_'+rls[4:]+'.csv')
        interp_func = interp1d(nir['Gmag'], nir['sigma'], kind='linear', fill_value="extrapolate")
        plx_unc = interp_func(G) #in micro-as
    else:
        plx_unc = np.sqrt(40 + 800 * z + 30 * z * z) * _t_factor[rls]

    ra_unc = plx_unc * pos_alpha_factor[rls[:3]]
    dec_unc = plx_unc * pos_delta_factor[rls[:3]]
    pmra_unc = plx_unc * pm_alpha_factor[rls[:3]]
    pmdec_unc = plx_unc * pm_delta_factor[rls[:3]]
    return plx_unc, ra_unc, dec_unc, pmra_unc, pmdec_unc

# conversion of parallax to distance, Weiler+25 (https://arxiv.org/abs/2505.16588)
def Weiler_C(x,p):
    return math.sqrt(2.) * scp.erfinv( 1. - p * ( 1. + scp.erf( x / math.sqrt(2.)) ) )
