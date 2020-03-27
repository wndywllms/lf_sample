# lf_sample

Tools to compute LF for generic radio-optical catalogues:

* `examples` directory contains 
`use_lf_sample.py` : with some example code using lf_sample

An `sample` instance needs to be initilised minimally with an identifier and catalogue. The catalogue is an astropy Table with fields `z`, `power`, `opt_lum`. The example contains a catalogue utility function 


`sample`  Routines:

* `calc_zmin_zmax` - calculates zmin and zmax needed for LF (will save this to files and read from a saved file if it exists)

* `sub_z_sample` - select a subsample on redshift (will automatically update zmin/zmax)

* `sub_sample_by_field` - select a subsample on anyother field present in the catalogue record array

* `get_LF_f_areal` - call to  `get_LF_f_areal` in `LF_util`


`util`  Routines:

* `get_LF_f_areal` - compute LF for given power bins, power, zmin, zmax, fcor, areal, area
