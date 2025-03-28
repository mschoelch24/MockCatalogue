# MockCatalogue
Python tool for creating mock catalogues of Gaia data.\
Generates an output of data including the Gaia observables 'source_id', 'ra', 'dec', 'parallax', 'parallax_error', 'pmra', 'pmdec', 'pmra_error', 'pmdec_error', 'radial_velocity', 'l', 'b', 'd', 'G', and 'Av', from an input simulation with parameters 'x', 'y', 'z', 'vx', 'vy', and 'vz'. Also outputs a mock catalogue with cartesian coordinates 'x', 'y', 'z', 'vx', 'vy', and 'vz' after applying Gaia uncertainties.  
Catalogues are adjustable for different Gaia error models, including GaiaNIR, assuming red clump (RC) stars (different absorption models and stellar tracers to be added).

## Files
- `exe_mockcat.sh`: Script to run `coord_transform.py`, `observables.py`, and `compile.py` in succession. Requires an 'input_file' edit. Current input file format options are *.fits*, *.pkl*, and *.csv*. The keyword 'rls' can be changed to apply different Gaia error models (dr3, dr4, dr5, NIR). To use a different number of cores or split a dataset into parts, 'ncores' and 'nruns' may be edited. 
- `cat_functions.py`: All functions used in the catalogue tool. Includes all coordinate transformations and the extinction, magnitude, and uncertainties computations (using performance models from [PyGaia](https://github.com/agabrown/PyGaia) and Hobbs et al. (in prep)).
- `coord_transform.py`: Reads input simulation, converts cartesian coordinates to equatorial heliocentric and galactic heliocentric. Outputs dataframe with simulation input coordinates, as well as 'source_id', 'ra', 'dec', 'parallax', 'pmra', 'pmdec', 'radial_velocity', 'l', 'b', 'd'.
- `observables.py`: Inputs 'l', 'b', 'd' from `coord_transform.py` to compute the extinction 'Av_lm' (using the Lallement 2022 & Marshall 2006 map), Gaia G magnitude 'G', and Gaia uncertainties in position and proper motion: 'parallax_error', 'ra_error', 'dec_error', 'pmra_error', and 'pmdec_error'.
- `extmap_Lallement22_Marshall2006.dat`: Grid of extinction values by distance, l, and b, from Marshall et al. (2006) and  Lallement et al. (2022). Exceeds upload limit in GitHub, access [here](https://drive.google.com/file/d/17iK-KLArShT6dU-B81VFPNgdRq2UQ8gM/view?usp=share_link) through Google Drive.
- `GmagSig_K5IIIAv0_*.csv`: Parallax uncertainties as a function of G magnitude for GaiaNIR (M-median mission, L-large mission, 5 or 10 year baseline) using simulator by Hobbs et al. (in prep) and assuming a RC stellar spectrum. 
- `compile.py`: Compiles the separate dataframes created during the repeated iteration over `observables.py`. Outputs both `_out.pkl`, which includes all Gaia observable parameters, and `_mock.pkl`, which provides the cartesian coordinates after applying Gaia uncertainties.
- Requires input simulation, for example [this Milky Way simulation](https://drive.google.com/file/d/1aRysr_2sHH2sOgYvQcIqF5UIKf26Qyok/view?usp=share_link).

## Usage
After downloading all files, make `exe_mockcat.sh` executable by running
```bash
chmod +x exe_mockcat.sh
```
In `exe_mockcat.sh`, edit 'input_file' with the name of your simulation file. If the input file does not have ordered columns 'x', 'y', 'z', 'vx', 'vy', and 'vz', edit the column identifiers in `coord_transform.py`.\
Lastly, run `exe_mockcat.sh`:
```bash
./exe_mockcat.sh
```

## Dependencies
- Python >= 3.7
- Numpy
- Pandas
- Astropy
- Scipy

## License
...
