from cat_functions import *
import time
import sys

def main():
    input_file = sys.argv[1]
    ncores = sys.argv[2]

    simname,_ = os.path.splitext(input_file)
    input_df = simname + '_coords.csv'

    t0 = time.time()
    df = pd.read_csv(input_df, index_col=0)
    print("Opening simulation",simname, "with columns", list(df),"and length", len(df))
    
    source_id = range(1, len(df) + 1) # Assigning a source id for sorting.
    df['source_id'] = source_id

    t1 = time.time()
    t_import = t1 - t0
    print("Import time: ", '%.2f' % (t_import/60) ,"min")

    
    print("****Starting extinction calculation****")
    extarray = compext_parallel(df['l'],df['b'],df['d'],df['source_id'], int(ncores))
    df['Av_lm'] = np.array(extarray)

    t2 = time.time()
    t_ext = t2 - t1
    print("Extinction computation finished in", '%.2f' % (t_ext/60) ,"min")

    
    print("****Starting G mag calculation****")
    G = magnitude(df['d'], df['Av_lm'])
    print("G mag min, max, median:", np.min(G), np.max(G), np.nanmedian(G))
    df['G'] = G

    t3 = time.time()
    t_Gmag = t3 - t1
    print("G mag computation finished in", '%.2f' % (t_ext/60) ,"min")

    
    print("****Starting uncertainties calculation****")
    plx_unc, pmra_unc, pmdec_unc = uncertainties(G)
    df['plx_unc'] = plx_unc
    df['pmra_unc'] = pmra_unc
    df['pmdec_unc'] = pmdec_unc

    df.to_csv(simname + '_out.csv') 
    print("Final dataframe contains columns", list(df), "and has length", len(df))

    t4 = time.time()
    t_total = t4 - t0
    print("Total time observables:", '%.2f' % (t_total/60) ,"min")

if __name__ == "__main__":
    main()
