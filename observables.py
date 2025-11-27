from cat_functions import *
import time
import sys

def main():
    input_file = sys.argv[1]
    ncores = sys.argv[2]
    nruns = int(sys.argv[3])
    tracer = sys.argv[4]
    rls = sys.argv[5]
    extinction = sys.argv[6]
    n = int(sys.argv[7])

    simname,_ = os.path.splitext(input_file)
    input_df = simname + '_coords.pkl'

    t0 = time.time()
    df = pd.read_pickle(input_df,compression='zip')
    
    # set limits of sources to process in this run:
    s_tot = len(df)
    s_i = n * (s_tot//nruns) - (s_tot//nruns)
    s_f = n * (s_tot//nruns) - 1

    # apply limits to the dataframe:
    df = pd.read_pickle(input_df,compression='zip')[s_i:s_f+1]
    print("Processing stars", s_i, "to", s_f, "in run", n , "/", nruns)

    source_id = range(s_i, s_f+1) # Assigning a source id for sorting.

    # New data frame to save only the observables.
    dfobs = pd.DataFrame()
    dfobs['source_id'] = source_id

    t1 = time.time()
    print("****Starting extinction calculation****")

    if extinction == 'gz':
        print("Using Green-Zucker extinction map.")
        from extinction_functions import extinction_gz
        extarray = extinction_gz(df['l'],df['b'],df['d'])#, df['source_id'], int(ncores))
    else:
        print("Using default Lallement-Marshall extinction map.")
        extarray = compext_parallel(df['l'],df['b'],df['d'],df['source_id'], int(ncores))
    dfobs['Av'] = np.array(extarray)

    t2 = time.time()
    t_ext = t2 - t1
    print("Extinction computation finished in", '%.2f' % (t_ext/60) ,"min")

    
    print("****Starting G mag calculation****")
    if tracer == 'RGB':
        G, bp_rp = magnitude_RGB(np.array(df['d']), np.array(dfobs['Av']))
        dfobs['bp_rp'] = bp_rp
        #print("G mag min, max, median:", np.min(G), np.max(G), np.nanmedian(G))
    else: 
        G = magnitude(np.array(df['d']), np.array(dfobs['Av']))
    dfobs['G'] = G

    t3 = time.time()
    t_Gmag = t3 - t2
    print("G mag computation finished in", '%.2f' % (t_Gmag/60) ,"min")

    
    print("****Starting uncertainties calculation****")
    plx_error, ra_error, dec_error, pmra_error, pmdec_error = uncertainties(G, rls)
    dfobs['plx_error'] = plx_error/1000 #converting to mas
    dfobs['ra_error'] = ra_error/1000
    dfobs['dec_error'] = dec_error/1000
    dfobs['pmra_error'] = pmra_error/1000 #mas/yr
    dfobs['pmdec_error'] = pmdec_error/1000

    dfobs.to_pickle(simname + '_observ_out_pt'+ str(n) +'.pkl',compression='zip')
    print("Observables dataframe no.", n ,"contains columns", list(dfobs), "and has length", len(dfobs))

    t4 = time.time()
    t_total = t4 - t0
    print("Total time observables.py {}/{}".format(n, nruns),":", '%.2f' % (t_total/60) ,"min")

if __name__ == "__main__":
    main()
