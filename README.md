# radio_lf

Tools to compute LF for generic radio-optical catalogues:

* `examples` directory contains 
`use_lf_sample.py` : with some example code using lf_sample

A `radio_lf.sample` instance needs to be initilised minimally with an identifier and catalogue. The catalogue is an astropy Table with fields `z`, `power`, `opt_lum`. The example contains a catalogue utility function 


`sample`  routines:

* `calc_zmin_zmax` - calculates zmin and zmax needed for LF (will save this to files and read from a saved file if it exists)

* `sub_z_sample` - select a subsample on redshift (will automatically update zmin/zmax)

* `sub_sample_by_field` - select a subsample on anyother field present in the catalogue record array

* `get_LF_f_areal` - call to  `get_LF_f_areal` in `LF_util`


`util`  contains various utilities including:

* class `rmsmapz` - to calculate volums from an rms object (fits image or histogram)

* `RadioPower`, `RadioFlux`, `OpticalLuminosity`, `OpticalFlux`

* `get_LF_f_areal` - compute LF for given power bins, power, zmin, zmax, fcor, areal, area


`model` contains functions that return literature and model LFs


# to install
run `python setup.py install --user`
