#!/bin/bash

input_file='Bar_ome45_short_carte.fits' #example input file, can be full path, e.g.
#input_file='../2Project/snapshots/PMcrs0a0.9000_wo_LMC_DM.pkl'

ncores='5' # input number of cores here
nruns='10' # input number of runs of observables.py

echo "input_file: $input_file"
echo "ncores: $ncores"

start=$(date +%s)

# running the coordinate transformations
python3 coord_transform.py "$input_file"

# adding other observables such as the G magnitude and uncertainties, and the Marshall & Lallement extinction model
N=$nruns
# looping this file for n times, include n as input so files can be labelled pt n
for ((i=1; i<=N; i++))
do
        echo "Executing iteration $i of observables.py"
        python3 observables.py "$input_file" "$ncores" "$nruns" $i
done

# compiling the separate dataframes from each iteration of observables.py and merging with the dataframe from coord_transform.py 
python3 compile.py "$input_file"

end=$(date +%s)

seconds=$(echo "$end - $start" | bc)
awk -v t=$seconds 'BEGIN{t=int(t*1000); printf "%d:%02d:%02d\n", t/3600000, t/60000%60, t/1000%60}'
