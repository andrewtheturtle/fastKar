# distance-contact lookup at 100kb-resolution
'medium_lookup'

# distance-contact lookup at 10kb-resolution
'small_lookup'

# distance-contact lookup at 1mb-resolution
'big_lookup'

# compartment lookup table
'comparment_lookup'

# Hi-C model kernels (scripts/modeltraining.R; GM12878 chr1-22 + X, kept 10 kb bins): distance_spline = log point density
# rho(d) per haploid copy at unit depth (cubic spline in log d to 1e8, 1/d beyond; integrated over each pixel by the C++ kernel);
# diag_spline = log diagonal count vs log bin area (linear)
'splineobj'

# inter-chromosomal contact density per haploid copy pair per bp^2 at unit (haploid) depth, genome-wide on kept bins
# (scripts/modeltraining.R)
'interchrom_density'