# MockCatalogue
Python tool for creating mock catalogues of Gaia data.\
Uses an input simulation with parameters `x`, `y`, `z`, `vx`, `vy`, and `vz` to generate an output of data, including the Gaia observables `source_id`, `ra`, `dec`,`ra_error`, `dec_error`, `parallax`, `parallax_error`, `pmra`, `pmdec`, `pmra_error`, `pmdec_error`, `radial_velocity`, `radial_velocity_error`, `l`, `b`, `d`, `G`, `bp_rp`, as well as the extinction `Av` and the new mock cartesian coordinates `X`, `Y`, `Z`, `Vx`, `Vy`, and `Vz`, after Gaia uncertainties have been applied. 
Catalogues are adjustable for different Gaia error models, including GaiaNIR, assuming red clump (RC) or red giant branch (RGB) stars (different absorption models and stellar tracers to be added in the future).

## Files
- `exe_mockcat.sh`: Script to run `coord_transform.py`, `observables.py`, and `compile.py` in succession. Requires an 'input_file' edit. Current input file format options are *.fits*, *.pkl*, and *.csv*. The keyword 'tracer' can be 'RC' or 'RGB', depending on the desired tracer stellar population and 'rls' can be changed to apply different Gaia error models (dr3, dr4, dr5, NIR). To use a different number of cores, 'ncores' and may be edited. If there is an existing version of a mock catalogue with a different error model than desired, the code will skip the computation of the observables (to avoid recomputing in particular the extinction `Av`, which is very computationally expensive) and run `uncertainties.py` instead. The output file will still include all observables and the mock catalogue with the new error model applied.
- Code requires input simulation, for example [this Milky Way simulation](https://drive.google.com/file/d/1aRysr_2sHH2sOgYvQcIqF5UIKf26Qyok/view?usp=share_link).
- `cat_functions.py`: All functions used in the catalogue tool. Includes all coordinate transformations and the extinction, magnitude, and uncertainties computations (using performance models from [PyGaia](https://github.com/agabrown/PyGaia) and Hobbs et al. (in prep)). For the RC stellar tracer, apparent magnitudes are calculated using absolute magnitude Mk = -1.62 and intrinsic color (J-K)0 = 0.55. For RGB stars, absolute magnitudes and colors are drawn from a KDE `kde_color_mag_samples.npz` to compute apparent magnitudes. Line-of-sight velocity uncertainties are computed using the [Gaia spectroscopic performance model](https://www.cosmos.esa.int/web/gaia/science-performance#spectroscopic%20performance) with G3 giant star parameters to represent the red giant branch (see Fig. E.1 in Katz et al. 2023).
G_RVS magnitudes of each particle are estimated using the Sartoretti et al. (2023) color-color relation, and the availability of line-of-sight velocities in the mock catalogues is limited by G_RVS ≤ 14 mag and 𝜎𝑉𝑙𝑜𝑠 ≤ 20 km s−1 for Gaia DR3 and G_RVS ≤ 16 mag for Gaia DR4 and beyond. Parallaxes are converted to distance with the Weiler_C function from [Weiler et al. (2025)](https://arxiv.org/abs/2505.16588).
- `coord_transform.py`: Reads input simulation, converts cartesian coordinates to equatorial heliocentric and galactic heliocentric. Outputs dataframe with simulation input coordinates, as well as `source_id`, `ra`, `dec`, `parallax`, `parallax_error`, `pmra`, `pmdec`, `pmra_error`, `pmdec_error`, `radial_velocity`, `radial_velocity_error`, `l`, `b`, and `d`.
- `observables.py`: Inputs `l`, `b`, `d` from `coord_transform.py` to compute the extinction `Av_lm` (using the Lallement 2022 & Marshall 2006 grid), Gaia magnitude `G` and color `bp_rp`, and Gaia uncertainties in position, proper motion, and line-of-sight velocity: `parallax_error`, `ra_error`, `dec_error`, `pmra_error`, `pmdec_error`, and `radial_velocity_error`. Note that all uncertainties are converted to mas, mas/yr, and km/s to match Gaia uncertainties.
- `extmap_Lallement22_Marshall2006.dat`: Grid of extinction values by distance, l, and b, from Marshall et al. (2006) and  Lallement et al. (2022). Also used in the Gaia Object Generator (GOG, Antiche et al. (2014)). Exceeds upload limit in GitHub, access [here](https://drive.google.com/file/d/17iK-KLArShT6dU-B81VFPNgdRq2UQ8gM/view?usp=share_link) through Google Drive.
- `GmagSig_K5IIIAv0_*.csv`: Parallax uncertainties as a function of G magnitude for GaiaNIR (M-median mission, L-large mission, 5 or 10 year baseline) using simulator by Hobbs et al. (in prep) and assuming a RC stellar spectrum. 
- `compile.py`: Compiles the separate dataframes created from the coordinate transformation (`coord_transform.py`) and addition of observables (`observables.py`). Also applies the uncertainties by drawing from a Gaussian centered on the true value, where the standard deviation is the respective true uncertainty. Distances computed using Weiler et al. (2025). Outputs file: `_mock_tracer_rls.pkl`, which includes the input coordinates, all Gaia observable parameters, the observables with uncertainties applied (capitalized parameters, e.g. `Ra`, `Dec`), and the new `X`, `Y`, `Z`, `Vx`, `Vy`, and `Vz` cartesian coordinates after Gaia uncertainties have been applied and heliocentric coordinates have been converted back to the cartesian coordinate frame.
- `uncertainties.py`: Applies the uncertainties of the specified error model and creates output mock catalogue (`_mock_tracer_rls.pkl`) if a previously computed mock catalogue was already available.
- `sky_projection.py`: Functions for projecting simulation into sky at LMC position and orientation. Contact [marie.schoelch@fqa.ub.edu](mailto:marie.schoelch@fqa.ub.edu) for supplemental code to create LMC mock catalogues.

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
