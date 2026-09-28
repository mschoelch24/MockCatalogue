#!/bin/bash

input_file='Bar_ome45_short_carte.fits' #example input file
tracer='RC' #stellar tracer ('RC' or 'RGB')
rls='dr3' #Gaia data release for uncertainties, for GaiaNIR: 'NIR_M5', 'NIR_M10', 'NIR_L5', or 'NIR_L10'
extinction='lm' #extinction model: 'gz'- Green+19 & Zucker+25 (lm to fill), 'lm'- Lallement+22 & Marshall+06

ncores='5' # input number of cores here

echo "input_file: $input_file"
echo "ncores: $ncores"
echo "rls: $rls"

start=$(date +%s)

output_pattern="${input_file%.*}_mock_${tracer}_*.pkl"
files=( $output_pattern )

# Check if another output file for this simulation already exists:
if [ ${#files[@]} -gt 0 ] && [ -f "${files[0]}" ]; then
    oldfile="${files[0]}"
    oldrls="${oldfile##*_mock_${tracer}_}"
    oldrls="${oldrls%.pkl}"

    echo "Found existing file: $oldfile"
    python3 uncertainties.py "$input_file" "$tracer" "$rls" "$oldrls"

    exit 0
fi

echo "No existing file found."
echo "Running full observable calculation..."

# running the coordinate transformations
python3 coord_transform.py "$input_file"

# adding other observables such as the G magnitude and uncertainties, and the Marshall & Lallement extinction model
python3 observables.py "$input_file" "$ncores" "$tracer" "$rls" "$extinction"

# compiling the separate dataframes from each iteration of observables.py and merging with the dataframe from coord_transform.py 
python3 compile.py "$input_file" "$tracer" "$rls"

end=$(date +%s)

seconds=$(echo "$end - $start" | bc)
awk -v t=$seconds 'BEGIN{t=int(t*1000); printf "%d:%02d:%02d\n", t/3600000, t/60000%60, t/1000%60}'

# removing all intermediate files
rm "${input_file%.*}"_coords.pkl
rm "${input_file%.*}"_observ.pkl
