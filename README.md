# MockCatalogue
Python tool for creating mock catalogues of Gaia data.

## Files
- exe_mockcat.sh: Script to run coord_transform.py and observables.py in succession. Requires an "input_file" and "ncores" edit. Current file format options are .fits, .pkl, and .csv.
- cat_functions.py: All functions used in the catalogue tool (todo: more documentation).
- coord_transform.py: Reads input simulation, converts cartesian coordinates to equatorial heliocentric (todo: optional inclusion in final df) and galactic heliocentric. Outputs dataframe with sim input coords, as well as 'source_id', 'ra', 'dec', 'parallax', 'pmra', 'pmdec', 'vr', 'l', 'b', 'd'.
- observables.py: Inputs 'l', 'b', 'd' from coord_transform.py output to compute the extinction 'Av_lm' (using the Lallement&Marshall map), Gaia G magnitude 'Gmag', and GaiaDR3 uncertainties in parallax and proper motion 'plx_unc', 'pmra_unc', 'pmdec_unc' in microarcseconds.
- combi_lallement22_Marshall2006.dat: Grid of extinction values by distance, l, and b, from Marshall 2006 and Lallement 2022. Exceeds upload limit in GitHub, access through Google Drive (https://drive.google.com/file/d/17iK-KLArShT6dU-B81VFPNgdRq2UQ8gM/view?usp=share_link)
- Requires input simulation, for example (https://drive.google.com/file/d/1aRysr_2sHH2sOgYvQcIqF5UIKf26Qyok/view?usp=share_link).

## Usage

Edit the "input_file" and "ncores" in the exe_mockcat.sh script and run.

## License
...
