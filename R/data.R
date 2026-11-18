# distance-contact lookup at 10kb-resolution
'small_lookup'

# compartment lookup table
'compartment_lookup'

# Hi-C model kernels (scripts/modeltraining.R; GM12878 chr1-22 + X, kept 10 kb bins): distance_spline = log point density
# rho(d) per haploid copy at unit depth (cubic spline in log d to 1e8, 1/d beyond; integrated over each pixel by the C++ kernel);
# diag_spline = log diagonal count vs log bin area (linear)
'splineobj'

# inter-chromosomal contact density per haploid copy pair per bp^2 at unit (haploid) depth, genome-wide on kept bins
# (scripts/modeltraining.R)
'interchrom_density'

# depth calibration of estimate.depthratio(): off-diagonal kept x kept 10 kb contacts (10 kb self-pixels excluded) per kept copy-bp per unit
# depth on GM12878, at the training depth (valid-pair coverage of kept bins) (scripts/modeltraining.R)
'depth_calibration'
