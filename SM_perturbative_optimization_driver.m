function SM_perturbative_optimization_driver(mode)
%SM_PERTURBATIVE_OPTIMIZATION_DRIVER  Reproduce the perturbative supplementary figure.
%
%   SM_perturbative_optimization_driver('figures')  rebuild the panels from the
%       cell data shipped with the repository.  Seconds.
%
%   SM_perturbative_optimization_driver('check')    run the validation tests:
%       slab mode, slab reflection against an independent transfer matrix,
%       analytic gradient against central finite differences, sampling
%       convergence, rank bound.  About two minutes.
%
%   SM_perturbative_optimization_driver('etch')     coupling efficiency versus
%       patterned-layer depth, the scan behind the choice of h_etch.  Minutes.
%
%   SM_perturbative_optimization_driver('quick')    recompute the ten NA = 0.1
%       cells from scratch and rebuild the panels.  Minutes.
%
%   SM_perturbative_optimization_driver('full')     recompute all thirty cells.
%       Hours to a day: cost grows as (NA*L)^2 in the number of free-space
%       modes, so the NA = 1.0 cells at L = 100*lambda need roughly 100 GB of
%       memory and several hours each.  Cells are independent, so on a cluster
%       run them as thirty separate tasks calling
%       SM_perturbative_optimization('cell', k) for k = 1..30.
%
%   Cells 1-10 are NA = 0.1, 11-20 are NA = 0.5, 21-30 are NA = 1.0, each with
%   L/lambda = 25 to 100 in ten steps.  Recomputing overwrites the shipped .mat
%   for that cell with an identical result.

if nargin < 1, mode = 'figures'; end

switch mode
    case 'figures'
        SM_perturbative_optimization('figures');
        SM_perturbative_optimization('paper');
    case 'check'
        SM_perturbative_optimization('tests');
    case 'etch'
        SM_perturbative_optimization('etch');
    case 'quick'
        for k = 1:10
            SM_perturbative_optimization('cell', k);
        end
        SM_perturbative_optimization('figures');
        SM_perturbative_optimization('paper');
    case 'full'
        for k = 1:30
            SM_perturbative_optimization('cell', k);
        end
        SM_perturbative_optimization('figures');
        SM_perturbative_optimization('paper');
    otherwise
        error('unknown mode ''%s''; see help SM_perturbative_optimization_driver', mode);
end
end
