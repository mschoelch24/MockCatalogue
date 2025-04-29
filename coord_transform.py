import time
import os.path
from cat_functions import *
import sys

def main():
    input_file = sys.argv[1]
    t0 = time.time()

    simname,_ = os.path.splitext(input_file)
    df = read_sim(input_file)
    print("Opening simulation",simname, "with columns", list(df),"and length", len(df))

    x = np.array(df['x']) #change if input columns are labeled differently, or use df.iloc[:,0] etc.
    y = np.array(df['y'])
    z = np.array(df['z'])

    vx = np.array(df['vx'])
    vy = np.array(df['vy'])
    vz = np.array(df['vz'])

    source_id = range(0,len(df)) # Assigning a source id for sorting.
    df['source_id'] = source_id

    t1 = time.time()
    t_import = t1 - t0
    print("Import time: ", '%.2f' % (t_import/60) ,"min")

    # adding rotation:
    #x,y,z,v_x,v_y,v_z = rotationz(x,y,z,v_x,v_y,v_z, 30)
    
    print("****Starting coordinate transformation****")

    # Coordinate transformation: cartesian galactocentric to equatorial heliocentric:
    ra, dec, distance, pmra, pmdec, radial_velocity = cartesian2equatorial(x,y,z,vx,vy,vz)

    df['ra'] = ra
    df['dec'] = dec
    df['parallax'] = 1/distance
    df['pmra'] = pmra
    df['pmdec'] = pmdec
    df['radial_velocity'] = radial_velocity

    # Coordinate transformation: cartesian galactocentric to galactic heliocentric:
    l_input, b_input, distance, pml, pmb, vr = cartesian2galactic(x,y,z,vx,vy,vz)

    df['l'] = l_input
    df['b'] = b_input
    df['d'] = distance

    print("Coordinate transformations complete, save df.")
    df.to_pickle(simname + '_coords.pkl',compression='zip')

    """
    del x, y, z, vx, vy, vz, ra, dec, parallax, pmra, pmdec, radial_velocity, pml, pmb
    gc.collect()
    """
    t4 = time.time()
    t_total = t4 - t0
    print("Total time coord_transform.py:", '%.2f' % (t_total/60) ,"min")

if __name__ == "__main__":
    main()


