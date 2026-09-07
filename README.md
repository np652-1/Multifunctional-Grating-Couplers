# Multifunctional-Grating-Couplers

## Supplementary: perturbative many-mode design and mode counts

`SM_perturbative_optimization.m` is the first-order perturbation-theory solver
behind the mode-count panels. It builds the free-space-to-guided scattering
matrix of a weakly patterned slab analytically, optimises the design's Fourier
coefficients to maximise the rank of that matrix, and counts the channels whose
coupling efficiency clears a threshold.

Everything is in that one file: a task switch at the top, subroutines below.

    SM_perturbative_optimization_driver('figures')   rebuild the panels from the
                                                     data in this repository
    SM_perturbative_optimization_driver('check')     validation tests
    SM_perturbative_optimization_driver('etch')      efficiency vs patterned-layer depth
    SM_perturbative_optimization_driver('quick')     recompute the NA = 0.1 cells
    SM_perturbative_optimization_driver('full')      recompute all thirty cells

`SM_perturbative_optimization_data_hetch0.075/` holds one `.mat` per swept cell,
each with the full coupling-efficiency spectrum, so the figures can be rebuilt
in seconds without recomputing. Cells 1-10 are NA = 0.1, 11-20 are NA = 0.5,
21-30 are NA = 1.0, each with `L/lambda` = 25 to 100 in ten steps.

The panels are written to `paper_panels/` at the exact pixel size of the images
they replace in the manuscript, with transparent backgrounds.

### Cost

The number of free-space modes grows as `(k0*NA*L)^2`, from 39 at NA = 0.1,
`L = 25*lambda` to about 63000 at NA = 1.0, `L = 100*lambda`. The NA = 0.1 cells
run in seconds; the largest NA = 1.0 cells take hours and roughly 100 GB. Cells
are independent, so on a cluster run them as separate tasks calling
`SM_perturbative_optimization('cell', k)`.

### Validation

`('check')` verifies the slab mode against its dispersion relation and unit
power normalisation, the slab reflection against an independent transfer-matrix
evaluation, the analytic figure-of-merit gradient against per-parameter central
finite differences, convergence in the guided-azimuth sampling, and the rank
bound. It also reports the first-Born validity diagnostic: the largest total
power any single free-space mode couples into guided channels, which must stay
well below one for first-order perturbation theory to apply.
