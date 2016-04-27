# NCL-tools
A few odd scripts for the NCAR Command Language (NCL) 

### Script `tropopause_pressure.ncl`  

This script takes every netcdf file in a directory and computes the pressure of the lapse-rate tropopause using the `trop_wmo.ncl` function, and produces a set of files where the `h1` in the filename is replaced by `tropopause`. 
It's designed to work with WACCM `h1` history fiels. 
Currently, all these output files hold is an array of tropopause pressures.

### Script `tropopause_pressure_DART.ncl`  

This script does the same thing as `tropopause_pressure.ncl`, but for the DART state-space output files (i.e. the ones with names like `Prior_Diag_STUFF.nc`. 
DART output files are pretty similar to WACCM/CAM files, but have some quirks, e.g.:  
+ variable dimensions are (time,copy,lat,lon,lev) whereas in WACCM/CAM it's (time, lev,lat,lon) - so the vertical level dimension is in  a different place, and there are additional copies for ensemble members, their mean and spread, etc. 
+ the reference pressure `P0` for computing pressure on hybrid model levels only has one value, but is in an array of as many vertical levels as the model. (Not sure if this is the case for all DART models -- but it is in DART-WACCM output from DART version kodiak). 

It doesn't make sense to compute tropopause height for copies like ensemble std or inflation factors, and in fact will yield crazy values that make it hard to view the file in ncview. T
Therefore (and so save computation time), this script only does the computation for the ensemble mean and individual ensemble members. 
