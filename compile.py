import pandas as pd
import glob
import sys
import time
import os

def main():
    input_file = sys.argv[1]
    simname,_ = os.path.splitext(input_file)
    
    t0 = time.time()
    
    # coordinate.py output dataframe:
    df1 = pd.read_csv(simname + '_coords.csv', index_col=0) 

    # adding all the observables from df fragments:
    fragmt_list = sorted(glob.glob('*observ_out_pt*.csv'))

    df2 = pd.read_csv(fragmt_list[0], index_col=0)
    for fragmt in fragmt_list[1:]:
        dfn = pd.read_csv(fragmt, index_col=0)
        df2 = pd.concat([df2, dfn], ignore_index=True)

    if len(df2) == len(df1):
        df1 = df1.merge(df2, how = 'inner', on='source_id')
    else:
        "Dataframes cannot be merged due to different lengths."
        pass

    df1.to_csv(simname + '_out.csv')
    print("Final dataframe contains columns", list(df1), "and has length", len(df1))

    tf = time.time()
    t_total = tf - t0
    print("Total time compile.py:", '%.2f' % (t_total/60) ,"min")

if __name__ == "__main__":
    main()


