#!/bin/bash

input_file='Bar_ome45_short_carte.fits' #example input file
tracer='RC' #stellar tracer ('RC' or 'RGB')
rls='dr3' #Gaia data release for uncertainties, for GaiaNIR: 'NIR_M5', 'NIR_M10', 'NIR_L5', or 'NIR_L10'

ncores='5' # input number of cores here
nruns='1' # input number of runs of observables.py

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
        python3 observables.py "$input_file" "$ncores" "$nruns" "$tracer" "$rls" $i
done

# compiling the separate dataframes from each iteration of observables.py and merging with the dataframe from coord_transform.py 
python3 compile.py "$input_file" "$tracer" "$rls"

end=$(date +%s)

seconds=$(echo "$end - $start" | bc)
awk -v t=$seconds 'BEGIN{t=int(t*1000); printf "%d:%02d:%02d\n", t/3600000, t/60000%60, t/1000%60}'

# removing all intermediate files
rm "${input_file%.*}"_coords.pkl
rm "${input_file%.*}_observ_out_pt"*.pkl
