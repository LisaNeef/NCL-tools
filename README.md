# NCL-tools
A few odd scripts for the NCAR Command Language (NCL) 

### `tropopause_pressure.ncl`  

This script takes every netcdf file in a directory and computes the pressure of the lapse-rate tropopause using the `trop_wmo.ncl` function, and produces a set of files where the `h1` in the filename is replaced by `tropopause`. 
It's designed to work with WACCM `h1` history fiels. 
Currently, all these output files hold is an array of tropopause pressures.
