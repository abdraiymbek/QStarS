Direction files were written with ndarray.tofile (always C order) but read by Fortran as (NT,ND) column-major: the applied perturbations were the transposes of the intended ones. Superseded.
