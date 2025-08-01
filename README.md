# MockCatalogue
Python tool for creating mock catalogues of Gaia data.\
Uses an input simulation with parameters `x`, `y`, `z`, `vx`, `vy`, and `vz` to generate an output of data, including the Gaia observables `source_id`, `ra`, `dec`,`ra_error`, `dec_error`, `parallax`, `parallax_error`, `pmra`, `pmdec`, `pmra_error`, `pmdec_error`, `radial_velocity`, `l`, `b`, `d`, `G`, as well as the extinction `Av` and the new mock cartesian coordinates `X`, `Y`, `Z`, `Vx`, `Vy`, and `Vz`, after Gaia uncertainties have been applied. 
Catalogues are adjustable for different Gaia error models, including GaiaNIR, assuming red clump (RC) stars (different absorption models and stellar tracers to be added in the future).

## Files
- `exe_mockcat.sh`: Script to run `coord_transform.py`, `observables.py`, and `compile.py` in succession. Requires an 'input_file' edit. Current input file format options are *.fits*, *.pkl*, and *.csv*. The keyword 'rls' can be changed to apply different Gaia error models (dr3, dr4, dr5, NIR). To use a different number of cores or split a dataset into parts, 'ncores' and 'nruns' may be edited. 
- `cat_functions.py`: All functions used in the catalogue tool. Includes all coordinate transformations and the extinction, magnitude, and uncertainties computations (using performance models from [PyGaia](https://github.com/agabrown/PyGaia) and Hobbs et al. (in prep)). Also includes the Weiler_C function from [Weiler et al. (2025)](https://arxiv.org/abs/2505.16588) for parallax to distance conversion.
- `coord_transform.py`: Reads input simulation, converts cartesian coordinates to equatorial heliocentric and galactic heliocentric. Outputs dataframe with simulation input coordinates, as well as `source_id`, `ra`, `dec`, `parallax`, `parallax_error`, `pmra`, `pmdec`, `pmra_error`, `pmdec_error`, `radial_velocity`, `l`, `b`, and `d`.
- `observables.py`: Inputs `l`, `b`, `d` from `coord_transform.py` to compute the extinction `Av_lm` (using the Lallement 2022 & Marshall 2006 grid), Gaia G magnitude `G`, and Gaia uncertainties in position and proper motion: `parallax_error`, `ra_error`, `dec_error`, `pmra_error`, and `pmdec_error`. Note that all uncertainties are converted to mas and mas/yr to match Gaia uncertainties.
- `extmap_Lallement22_Marshall2006.dat`: Grid of extinction values by distance, l, and b, from Marshall et al. (2006) and  Lallement et al. (2022). Also used in the Gaia Object Generator (GOG, Antiche et al. (2014)). Exceeds upload limit in GitHub, access [here](https://drive.google.com/file/d/17iK-KLArShT6dU-B81VFPNgdRq2UQ8gM/view?usp=share_link) through Google Drive.
- `GmagSig_K5IIIAv0_*.csv`: Parallax uncertainties as a function of G magnitude for GaiaNIR (M-median mission, L-large mission, 5 or 10 year baseline) using simulator by Hobbs et al. (in prep) and assuming a RC stellar spectrum. 
- `compile.py`: Compiles the separate dataframes created during the repeated iteration over `observables.py`. Also applies the uncertainties by drawing from a Gaussian, and computes distances using Weiler et al. (2025). Outputs file: `_mock.pkl`, which includes the input coordinates, all Gaia observable parameters, and the new `X`, `Y`, `Z`, `Vx`, `Vy`, and `Vz` cartesian coordinates after Gaia uncertainties have been applied and heliocentric coordinates have been converted back to the cartesian coordinate frame.
- Requires input simulation, for example [this Milky Way simulation](https://drive.google.com/file/d/1aRysr_2sHH2sOgYvQcIqF5UIKf26Qyok/view?usp=share_link).

## Usage
After downloading all files, make `exe_mockcat.sh` executable by running
```bash
chmod +x exe_mockcat.sh
```
In `exe_mockcat.sh`, edit 'input_file' with the name of your simulation file. If the input file does not have ordered columns `x`, `y`, `z`, `vx`, `vy`, and `vz`, edit the column identifiers in `coord_transform.py`.\
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
