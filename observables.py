from cat_functions import *
import time
import sys

def main():
    input_file = sys.argv[1]
    ncores = sys.argv[2]
    nruns = int(sys.argv[3])
    n = int(sys.argv[4])

    simname,_ = os.path.splitext(input_file)
    input_df = simname + '_coords.csv'

    t0 = time.time()
    df = pd.read_csv(input_df, index_col=0)
    
    # set limits of sources to process in this run:
    s_tot = len(df)
    s_i = n * (s_tot//nruns) - (s_tot//nruns)
    s_f = n * (s_tot//nruns)

    # apply limits to the dataframe:
    df = pd.read_csv(input_df, index_col=0)[s_i:s_f]
    print("Processing stars", s_i, "to", s_f, "in run", n , "/", nruns)

    source_id = range(s_i, s_f) # Assigning a source id for sorting.

    # New data frame to save only the observables.
    dfobs = pd.DataFrame()
    dfobs['source_id'] = source_id

    t1 = time.time()
    print("****Starting extinction calculation****")
    extarray = compext_parallel(df['l'],df['b'],df['d'],df['source_id'], int(ncores))
    dfobs['Av_lm'] = np.array(extarray)

    t2 = time.time()
    t_ext = t2 - t1
    print("Extinction computation finished in", '%.2f' % (t_ext/60) ,"min")

    
    print("****Starting G mag calculation****")
    G = magnitude(np.array(df['d']), np.array(dfobs['Av_lm']))
    #print("G mag min, max, median:", np.min(G), np.max(G), np.nanmedian(G))
    dfobs['G'] = G

    t3 = time.time()
    t_Gmag = t3 - t1
    print("G mag computation finished in", '%.2f' % (t_ext/60) ,"min")

    
    print("****Starting uncertainties calculation****")
    plx_error, pmra_error, pmdec_error = uncertainties(G)
    dfobs['plx_error'] = plx_error
    dfobs['pmra_error'] = pmra_error
    dfobs['pmdec_error'] = pmdec_error

    dfobs.to_csv(simname + '_observ_out_pt'+ str(n) +'.csv')
    print("Observables dataframe no.", n ,"contains columns", list(dfobs), "and has length", len(dfobs))

    t4 = time.time()
    t_total = t4 - t0
    print("Total time observables.py {}/{}".format(n, nruns),":", '%.2f' % (t_total/60) ,"min")

if __name__ == "__main__":
    main()
