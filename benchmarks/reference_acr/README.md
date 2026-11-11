Reference kinetic-module and ACR/ACRR stage, run with the adapted R code, on the six models whose
concordance is sound or adjudicated. Inputs were the completed exhaustive .mat files with metabolite
names patched in (patch_mat_names.m) and exact full-coupling matrices from coupling_exact/.

Overview_*.csv columns of interest: [4] modules, [5] singletons, [6] largest module, [14] ACR count,
[15] ACRR count -- note a value of 1 in [15] means exactly one pair. MetDouble_*.csv lists the ACRR
pairs; a single pair is written as a two-row column. Every pair listed is in Kimonu's `acrr` by name.
