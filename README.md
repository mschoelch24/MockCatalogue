# MockCatalogue
Python tool for creating mock catalogues of Gaia data.\
Generates an output of data including the Gaia observables 'source_id', 'ra', 'dec', 'parallax', 'parallax_error', 'pmra', 'pmdec', 'pmra_error', 'pmdec_error', 'radial_velocity', 'l', 'b', 'd', 'G', and 'Av', from an input simulation with parameters 'x', 'y', 'z', 'vx', 'vy', and 'vz'. Also outputs a file '_mock.pkl' with cartesian coordinates 'x', 'y', 'z', 'vx', 'vy', and 'vz' AFTER applying the Gaia uncertainties.  
Catalogues are adjustable for different Gaia error models and GaiaNIR assuming red clump stars (different absorption models and stellar tracers to be added).

## Files
- `exe_mockcat.sh`: Script to run `coord_transform.py`, `observables.py`, and `compile.py` in succession. Requires an "input_file" and "ncores" edit. The keyword 'rls' for release can also be added to apply different Gaia data release uncertainties. To split a dataset into parts, "nruns" may be edited. Current file format options are *.fits*, *.pkl*, and *.csv*.
- `cat_functions.py`: All functions used in the catalogue tool.
- `coord_transform.py`: Reads input simulation, converts cartesian coordinates to equatorial heliocentric and galactic heliocentric. Outputs dataframe with sim input coordinates, as well as 'source_id', 'ra', 'dec', 'parallax', 'pmra', 'pmdec', 'radial_velocity', 'l', 'b', 'd'.
- `observables.py`: Inputs 'l', 'b', 'd' from `coord_transform.py` to compute the extinction 'Av_lm' (using the Lallement 2022 & Marshall 2006 map), Gaia G magnitude 'G', and GaiaDR3 uncertainties in parallax and proper motion: 'parallax_error', 'pmra_error', and 'pmdec_error', all in microarcseconds.
- `extmap_Lallement22_Marshall2006.dat`: Grid of extinction values by distance, l, and b, from Marshall et al. (2006) and  Lallement et al. (2022). Exceeds upload limit in GitHub, access [here](https://drive.google.com/file/d/17iK-KLArShT6dU-B81VFPNgdRq2UQ8gM/view?usp=share_link) through Google Drive.
- `compile.py`: Compiles the separate dataframes created during the repeated iteration over `observables.py` to split the dataset into smaller parts.
- Requires input simulation, for example [this Milky Way simulation](https://drive.google.com/file/d/1aRysr_2sHH2sOgYvQcIqF5UIKf26Qyok/view?usp=share_link).

## Usage
After downloading all files, make `exe_mockcat.sh` executable by running
```bash
chmod +x exe_mockcat.sh
```
In `exe_mockcat.sh`, edit the "input_file" with the name of the simulation file. Add "ncores" depending on the cores available, and change "nruns" to split the full source sample into fragments to be run in `observables.py`. If the input file does not have the ordered columns 'x', 'y', 'z', 'vx', 'vy', and 'vz', edit the column identifiers in `coord_transform.py`.\
Finally, run `exe_mockcat.sh`:
```bash
./exe_mockcat.sh
```

## Dependencies
- Python >= 3.7
- Numpy
- Pandas
- [PyGaia](https://github.com/agabrown/PyGaia)

## License
...
