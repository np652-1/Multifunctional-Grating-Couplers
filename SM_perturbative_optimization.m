function SM_perturbative_optimization(task, idx)
%SM_PERTURBATIVE_OPTIMIZATION  Perturbative many-mode grating-coupler mode counts.
%
%   SM_perturbative_optimization              run every missing sweep cell, then figures
%   SM_perturbative_optimization('cell', k)   run sweep cell k only (batch worker)
%   SM_perturbative_optimization('figures')   rebuild figures from cached cell data
%   SM_perturbative_optimization('tests')     run the validation tests
%   SM_perturbative_optimization('etch')      coupling efficiency vs patterned-layer depth
%   SM_perturbative_optimization('paper')     drop-in replacements for the manuscript panels
%
%   Cell results are cached in SM_perturbative_optimization_data/.

if nargin < 1, task = 'all'; end

cfg   = configuration();
cells = sweep_cells(cfg);

if ~exist(cfg.datadir, 'dir'), mkdir(cfg.datadir); end

switch task
    case 'cell'
        run_cell(cfg, cells(idx));
    case 'figures'
        make_figures(cfg, cells);
    case 'tests'
        run_tests(cfg);
    case 'etch'
        etch_scan(cfg);
    case 'paper'
        paper_panels(cfg, cells);
    case 'all'
        for ic = 1:numel(cells)
            if ~exist(cell_file(cfg, cells(ic)), 'file')
                run_cell(cfg, cells(ic));
            end
        end
        make_figures(cfg, cells);
    otherwise
        error('unknown task ''%s''', task);
end
end

% ---------------------------------------------------------------- parameters

function cfg = configuration()
cfg.lambda      = 1;
cfg.k0          = 2*pi/cfg.lambda;
cfg.t_wg        = 0.3;              % slab thickness
cfg.neff        = 3.0;              % target effective index, sets n_core
cfg.h_etch      = 0.075;            % patterned-layer thickness
cfg.overlay     = true;             % pattern sits above the core, in air
cfg.nz          = 200;              % vertical quadrature points

cfg.NA_list     = [0.1, 0.5, 1.0];
cfg.L_list      = round(linspace(25, 100, 10));
cfg.nrestart    = 5;

cfg.oversample  = 3;                % guided azimuths per independent channel
cfg.kz_min      = 0.05;             % drop free-space modes within ~3 deg of grazing,
                                    % where the unit-power normalisation diverges
cfg.alpha       = 1e-1;             % Tikhonov regularisation in the figure of merit
cfg.maxiter     = 50;
cfg.CE_cutoff   = 0.02;

cfg.datadir     = sprintf('SM_perturbative_optimization_data_hetch%.3f', cfg.h_etch);
cfg.figdir      = sprintf('SM_perturbative_optimization_figures_hetch%.3f', cfg.h_etch);
cfg.ref_NA      = 0.1;              % cell shown in the design and spectrum panels
cfg.ref_L       = 100;
end

function cells = sweep_cells(cfg)
cells = struct('NA', {}, 'L', {}, 'nfour', {});
for NA = cfg.NA_list
    for L = cfg.L_list
        cells(end+1) = struct('NA', NA, 'L', L, ...
                              'nfour', n_fourier(NA, cfg.neff)/2); %#ok<AGROW>
    end
end
end

function f = cell_file(cfg, c)
f = fullfile(cfg.datadir, sprintf('cell_NA%.2f_L%03d.mat', c.NA, c.L));
end

% ------------------------------------------------------------------ geometry

function g = geometry(cfg)
res       = @(n) slab_te_beta(cfg.t_wg, cfg.k0, n) - cfg.neff*cfg.k0;
g.n_core  = fzero(res, [1.01, 6]);
g.chimax  = g.n_core^2 - 1;
[g.beta, g.Efun] = slab_te_mode(cfg.t_wg, cfg.k0, g.n_core);
g.z       = linspace(cfg.t_wg/2, cfg.t_wg/2 + cfg.h_etch, cfg.nz);
g.Emode   = g.Efun(g.z);
end

function beta = slab_te_beta(t, k0, n_core)
[beta, ~] = slab_te_mode(t, k0, n_core);
end

function [beta, Efun] = slab_te_mode(t, k0, n_core)
% Fundamental TE mode of a symmetric slab in air, normalised to unit power.
a = t/2;
v = k0*a*sqrt(n_core^2 - 1);
b = fzero(@(b) v*sqrt(1-b) - atan(sqrt(b/(1-b))), [1e-12, 1-1e-12]);
beta  = sqrt(b*(n_core^2-1) + 1)*k0;
kappa = v*sqrt(1-b)/a;
xi    = v*sqrt(b)/a;
prof  = @(z) (abs(z) <= a).*cos(kappa*min(max(z,-a),a)) + ...
             (abs(z) >  a).*cos(kappa*a).*exp(-xi*(abs(z)-a));
zq = linspace(-a-8/xi, a+8/xi, 4000);
P  = 0.5*(beta/k0)*trapz(zq, prof(zq).^2);
Efun = @(z) prof(z)/sqrt(P);
end

% ------------------------------------------------------- incident plane wave

function [A, inplane] = incident_profile(cfg, g, kpar, pol)
% Vertical profile and in-plane amplitude of a unit-incident-power plane wave
% in the presence of the unpatterned slab, evaluated on g.z.
% kpar is a column of in-plane wavenumber magnitudes; A is numel(kpar)-by-nz.
k0   = cfg.k0;
kz0  = sqrt(max(k0^2 - kpar.^2, 0));
kz1  = sqrt(g.n_core^2*k0^2 - kpar.^2);

switch pol
    case 's', Y0 = kz0;                 Y1 = kz1;
    case 'p', Y0 = kz0;                 Y1 = kz1/g.n_core^2;
end
r01 = (Y0 - Y1)./(Y0 + Y1);
ph  = exp(2i*kz1*cfg.t_wg);
r   = r01.*(1 - ph)./(1 - r01.^2.*ph);

zeta = g.z - cfg.t_wg/2;
switch pol
    case 's'
        A = exp(-1i*kz0*zeta) + r.*exp(1i*kz0*zeta);
        inplane = ones(size(kpar));
    case 'p'
        A = (kz0/k0).*(exp(-1i*kz0*zeta) - r.*exp(1i*kz0*zeta));
        inplane = ones(size(kpar));
end
A = A ./ sqrt(0.5*(kz0/k0));   % unit incident power per unit area; the
                               % caller divides by L for the L-by-L aperture
end

% ---------------------------------------------------------------------- bases

function fs = free_space_basis(cfg, g, NA, L)
% Free-space modes on the aperture-orthogonal grid, both polarisations.
dq = 2*pi/L;
kmax = cfg.k0*min(NA, sqrt(1 - cfg.kz_min^2));
n  = ceil(kmax/dq);
[QX, QY] = meshgrid((-n:n)*dq, (-n:n)*dq);
keep = (QX.^2 + QY.^2) <= kmax^2;
qx = QX(keep).'; qy = QY(keep).';

fs.kx  = [qx, qx];
fs.ky  = [qy, qy];
fs.pol = [repmat('s', 1, numel(qx)), repmat('p', 1, numel(qx))];

% The free-space grid is rectangular, so the separable aperture transform
% only needs one-dimensional sinc tables over the distinct kx and ky values.
[fs.kxu, ~, ixu] = unique(fs.kx);
[fs.kyu, ~, iyu] = unique(fs.ky);
fs.ix = ixu(:).';  fs.iy = iyu(:).';

kpar = hypot(qx, qy).';
kpar(kpar < 1e-12) = 1e-12;
khx  = (qx.')./kpar;  khy = (qy.')./kpar;

[As, ~] = incident_profile(cfg, g, kpar, 's');
[Ap, ~] = incident_profile(cfg, g, kpar, 'p');

Iz_s = trapz(g.z, As .* g.Emode, 2);
Iz_p = trapz(g.z, Ap .* g.Emode, 2);

% in-plane unit vectors: s along zhat x khat, p along khat.  The 1/L
% completes the normalisation to unit incident power over the aperture.
fs.ex = [(-khy).*Iz_s;  khx.*Iz_p].'/L;
fs.ey = [( khx).*Iz_s;  khy.*Iz_p].'/L;
end

function ov = oversample_for(cfg, NA, L)
% Guided-azimuth sampling. Halving it from 3 to 2 changes the mode count by
% ~1% (222 vs 225 at NA=0.1, L=100) and is used where assembly dominates.
if cfg.k0*NA*L > 300, ov = 2; else, ov = cfg.oversample; end
end

function wg = waveguide_basis(cfg, g, NA, L, thetas)
% Guided azimuths reachable from within the numerical aperture, one arc per
% Fourier component, sampled at cfg.oversample points per independent channel.
half  = asin(min(NA/cfg.neff, 1));
dphi  = 2*pi/(g.beta*L)/oversample_for(cfg, NA, L);
m     = round(half/dphi);
off   = (-m:m)*dphi;
phi   = [];
for th = [thetas, thetas + pi]
    phi = [phi, th + off]; %#ok<AGROW>
end
wg.phi   = phi;
wg.dphi  = dphi;
wg.tx    = -sin(phi);
wg.ty    =  cos(phi);
end

% ---------------------------------------------------------------- the matrix

function W = coupling_weight(cfg, g, fs, wg)
W = coupling_constant(cfg, g)*sqrt(wg.dphi) * (wg.tx.'.*fs.ex + wg.ty.'.*fs.ey);
end

function [AX, AY, dAX, dAY] = sinc_tables(g, fs, wg, qx, qy, L)
% One-dimensional sinc tables over the distinct grid wavenumbers.  This is
% the whole reason the assembly is cheap: the transcendental work is
% Nwg*(nkx+nky) rather than Nwg*Nfree.
tx = (fs.kxu - g.beta*cos(wg.phi).' + qx)*L/2;
ty = (fs.kyu - g.beta*sin(wg.phi).' + qy)*L/2;
AX = sinc0(tx);  AY = sinc0(ty);
if nargout > 2, dAX = dsinc0(tx)*(L/2);  dAY = dsinc0(ty)*(L/2); end
end

function F = expand(AX, AY, fs, L)
F = L^2 * AX(:, fs.ix) .* AY(:, fs.iy);
end

function S = smatrix(cfg, g, fs, wg, thetas, amps, dq, L)
% Perturbative scattering matrix: guided channels (rows) from free-space modes.
nf  = numel(thetas);
cAC = g.chimax/(2*nf);
qr  = g.beta + dq;
[AX, AY] = sinc_tables(g, fs, wg, 0, 0, L);
S = (g.chimax/2)*expand(AX, AY, fs, L);
for n = 1:nf
    qx = qr(n)*cos(thetas(n));  qy = qr(n)*sin(thetas(n));
    [AXp, AYp] = sinc_tables(g, fs, wg,  qx,  qy, L);
    [AXm, AYm] = sinc_tables(g, fs, wg, -qx, -qy, L);
    S = S + amps(n)*0.5*cAC*(expand(AXp, AYp, fs, L) + expand(AXm, AYm, fs, L));
end
S = S .* coupling_weight(cfg, g, fs, wg);
end

function grad = smatrix_contract(P, cfg, g, fs, wg, thetas, amps, dq, L)
% grad(p) = -2*Re sum( dS/dp .* P ), one Fourier component at a time so the
% full derivative stack is never stored.
nf   = numel(thetas);
cAC  = g.chimax/(2*nf);
qr   = g.beta + dq;
PW   = P.*coupling_weight(cfg, g, fs, wg);
grad = zeros(1, 2*nf);
for n = 1:nf
    qx = qr(n)*cos(thetas(n));  qy = qr(n)*sin(thetas(n));
    [AXp, AYp, dAXp, dAYp] = sinc_tables(g, fs, wg,  qx,  qy, L);
    [AXm, AYm, dAXm, dAYm] = sinc_tables(g, fs, wg, -qx, -qy, L);
    Fp = expand(AXp, AYp, fs, L);  Fm = expand(AXm, AYm, fs, L);
    grad(n) = -2*real(sum(0.5*cAC*(Fp + Fm).*PW, 'all'));
    dp = cos(thetas(n))*expand(dAXp, AYp, fs, L) + sin(thetas(n))*expand(AXp, dAYp, fs, L);
    dm = cos(thetas(n))*expand(dAXm, AYm, fs, L) + sin(thetas(n))*expand(AXm, dAYm, fs, L);
    grad(nf + n) = -2*real(sum(amps(n)*0.5*cAC*(dp - dm).*PW, 'all'));
end
end

function K = coupling_constant(cfg, g)
% |c| = K * int d3r' Emode(z') (Einc . thetahat) deltaeps ;  see SM Sec. I B.
K = cfg.k0^2/4 * g.beta/(2*cfg.k0) * sqrt(2/(pi*g.beta));
end

function v = sinc0(t)
v = ones(size(t));
m = abs(t) > 1e-9;
v(m) = sin(t(m))./t(m);
end

function v = dsinc0(t)
v = zeros(size(t));
m = abs(t) > 1e-9;
v(m) = (t(m).*cos(t(m)) - sin(t(m)))./t(m).^2;
end

% -------------------------------------------------------- figure of merit

function [fom, grad, CE] = fom_and_grad(params, cfg, g, fs, wg, thetas, L)
nf   = numel(thetas);
amps = params(1:nf);
dq   = params(nf+1:end);
S    = smatrix(cfg, g, fs, wg, thetas, amps, dq, L);

[nr, nc] = size(S);
byrow = nr <= nc;                       % diagonalise the smaller Gram matrix
if byrow, A = S*S'; else, A = S'*S; end
A = (A + A')/2;
[U, D] = eig(A, 'vector');
D = max(real(D), 0);
w = 1./(D + cfg.alpha);
fom = sum(w);
if nargout > 2, CE = sort(D, 'descend'); end

if nargout > 1
    if byrow
        B = S'*U;                       % nc-by-nr
        P = conj(U)*(diag(w.^2)*B.');   % nr-by-nc
    else
        B = S*U;                        % nr-by-nc
        P = conj(B)*(diag(w.^2)*U.');   % nr-by-nc
    end
    grad = smatrix_contract(P, cfg, g, fs, wg, thetas, amps, dq, L);
end
end

% ------------------------------------------------------------- optimisation

function run_cell(cfg, c)
g  = geometry(cfg);
fs = free_space_basis(cfg, g, c.NA, c.L);
thetas = (0:c.nfour-1)*pi/c.nfour;
wg = waveguide_basis(cfg, g, c.NA, c.L, thetas);

kperp = cfg.k0*c.NA;
lb = [-ones(1, c.nfour), zeros(1, c.nfour)];
ub = [ ones(1, c.nfour), kperp*ones(1, c.nfour)];

opts = optimoptions('fmincon', 'Algorithm', 'sqp', ...
    'SpecifyObjectiveGradient', true, 'MaxIterations', cfg.maxiter, ...
    'Display', 'iter', 'OptimalityTolerance', 1e-10);

[nrestart, maxiter] = effort(numel(wg.phi)*numel(fs.kx), cfg);
opts = optimoptions(opts, 'MaxIterations', maxiter);

best = struct('fom', inf);
for r = 1:nrestart
    rng(1000*round(100*c.NA) + c.L + r);
    p0 = [0.8 + 0.1*rand(1, c.nfour), (0.3 + 0.2*rand)*kperp*ones(1, c.nfour)];
    [p, f] = fmincon(@(p) fom_and_grad(p, cfg, g, fs, wg, thetas, c.L), ...
                     p0, [], [], [], [], lb, ub, [], opts);
    if f < best.fom
        best.fom = f; best.params = p; best.restart = r;
    end
end

[~, ~, CE] = fom_and_grad(best.params, cfg, g, fs, wg, thetas, c.L);
S = smatrix(cfg, g, fs, wg, thetas, best.params(1:c.nfour), ...
            best.params(c.nfour+1:end), c.L);

out = struct('cell', c, 'cfg', cfg, 'thetas', thetas, 'params', best.params, ...
             'fom', best.fom, 'restart', best.restart, 'CE', CE, ...
             'nrestart', nrestart, 'maxiter', maxiter, ...
             'Nmodes', sum(CE >= cfg.CE_cutoff), ...
             'Nrank', sum(CE >= 1e-6*max(CE)), 'CEmax', max(CE), ...
             'Nphase', phase_matched_estimate(c.NA, c.L, cfg.neff, g.beta), ...
             'power_bound', max(sum(abs(S).^2, 1)), ...
             'Nfree', size(S,2), 'Nwg', size(S,1));
save(cell_file(cfg, c), '-struct', 'out');
fprintf(['NA=%.2f L=%3d : %4d above %.0f%%, rank %4d, ceiling %4d, ' ...
         'CEmax %.3e, Nfree %d, Nwg %d, max column power %.3e%s\n'], ...
        c.NA, c.L, out.Nmodes, 100*cfg.CE_cutoff, out.Nrank, out.Nphase, ...
        out.CEmax, out.Nfree, out.Nwg, out.power_bound, born_note(out.power_bound));
end

function [nrestart, maxiter] = effort(work, cfg)
% Restart and iteration budget scaled to the size of one matrix assembly.
if     work < 1e7, nrestart = cfg.nrestart; maxiter = cfg.maxiter;
elseif work < 1e8, nrestart = 3;            maxiter = 40;
elseif work < 1e8*2, nrestart = 2;          maxiter = 30;
else,              nrestart = 2;            maxiter = 20;
end
end

function N = n_fourier(NA, neff)
N = 2*floor(pi/(asin(1/neff) + asin(min(NA/neff,1))));
end

function N = phase_matched_estimate(NA, L, neff, beta)
% Leading-order estimate of the number of STRONGLY coupled guided channels:
% those whose phase matching falls inside the numerical aperture.  This is an
% estimate, not a bound.  The aperture transform has algebraic sinc tails, so
% every free-space mode couples weakly to every guided azimuth and the rank of
% S is not confined to these arcs; the only true bound is min(Nwg, Nfree).
k0   = beta/neff;
arc  = n_fourier(NA, neff)*asin(min(NA/neff,1))*beta*L/pi;
dq   = 2*pi/L;
kmax = k0*min(NA, sqrt(1 - 0.05^2));
n    = ceil(kmax/dq);
[QX, QY] = meshgrid((-n:n)*dq, (-n:n)*dq);
disk = 2*nnz(QX.^2 + QY.^2 <= kmax^2);
N = min(round(arc), disk);
end

% ----------------------------------------------------------------- figures

function make_figures(cfg, cells)
if ~exist(cfg.figdir, 'dir'), mkdir(cfg.figdir); end
g = geometry(cfg);

ref = [];
for ic = 1:numel(cells)
    if cells(ic).NA == cfg.ref_NA && cells(ic).L == cfg.ref_L
        f = cell_file(cfg, cells(ic));
        if exist(f, 'file'), ref = load(f); end
    end
end

if ~isempty(ref)
    panel_design(cfg, g, ref);
    panel_spectrum(cfg, ref);
    panel_directions(cfg, g, ref);
end
panel_scaling(cfg, cells, g);
end

function panel_design(cfg, g, ref)
nf = ref.cell.nfour;
amps = ref.params(1:nf);  dq = ref.params(nf+1:end);
w = 10;
t = linspace(-w/2, w/2, 800);
[X, Y] = meshgrid(t, t);
chi = g.chimax/2 + zeros(size(X));
for n = 1:nf
    qr = g.beta + dq(n);
    chi = chi + amps(n)*(g.chimax/(2*nf))*cos(qr*(cos(ref.thetas(n))*X + sin(ref.thetas(n))*Y));
end

fig = figure('visible', 'off', 'position', [0 0 900 800]);
imagesc(t, t, sqrt(1 + chi)); axis xy image; colorbar
xlabel('x/\lambda'); ylabel('y/\lambda'); title('refractive index')
set(gca, 'fontsize', 16)
print(fig, fullfile(cfg.figdir, 'panel_a_design.png'), '-dpng', '-r150'); close(fig)

F = abs(fftshift(fft2(chi - mean(chi(:)))));
kx = 2*pi*((0:numel(t)-1) - floor(numel(t)/2))/(t(end)-t(1));
fig = figure('visible', 'off', 'position', [0 0 900 800]);
imagesc(kx, kx, F); axis xy image
xlim(1.6*g.beta*[-1 1]); ylim(1.6*g.beta*[-1 1]); colorbar
xlabel('q_x'); ylabel('q_y'); title('|FT \chi|')
set(gca, 'fontsize', 16)
print(fig, fullfile(cfg.figdir, 'panel_b_fourier.png'), '-dpng', '-r150'); close(fig)
end

function panel_spectrum(cfg, ref)
CE = ref.CE;
fig = figure('visible', 'off', 'position', [0 0 1000 700]);
semilogy(1:numel(CE), CE, 'linewidth', 3, 'color', [0.85 0.33 0.10]); hold on
yline(cfg.CE_cutoff, 'r-', 'linewidth', 2);
xline(ref.Nmodes, 'r--', 'linewidth', 2);
hold off; grid on; box on
xlim([1, max(2*ref.Nmodes, 10)]); ylim([1e-6, 1])
xlabel('mode number'); ylabel('coupling efficiency')
title(sprintf('NA = %.2f, L = %d\\lambda: %d modes above %.0f%%', ...
      ref.cell.NA, ref.cell.L, ref.Nmodes, 100*cfg.CE_cutoff))
set(gca, 'fontsize', 16)
print(fig, fullfile(cfg.figdir, 'panel_d_spectrum.png'), '-dpng', '-r150'); close(fig)
end

function panel_directions(cfg, g, ref)
nf = ref.cell.nfour;
fs = free_space_basis(cfg, g, ref.cell.NA, ref.cell.L);
wg = waveguide_basis(cfg, g, ref.cell.NA, ref.cell.L, ref.thetas);
S  = smatrix(cfg, g, fs, wg, ref.thetas, ref.params(1:nf), ref.params(nf+1:end), ref.cell.L);
[~, ibest] = max(sum(abs(S).^2, 1));
p = abs(S(:, ibest)).^2;

% total coupled power in each of the 2*nfour arcs, and the fine structure
edges = [ref.thetas, ref.thetas + pi];
arcP = zeros(size(edges));
for m = 1:numel(edges)
    d = mod(wg.phi - edges(m) + pi, 2*pi) - pi;
    arcP(m) = sum(p(abs(d) < asin(min(ref.cell.NA/cfg.neff,1))));
end

fig = figure('visible', 'off', 'position', [0 0 1000 900]);
pax = polaraxes(fig); hold(pax, 'on')
for m = 1:numel(edges)
    polarplot(pax, [edges(m) edges(m)], [0 arcP(m)/max(arcP)], ...
              'linewidth', 5, 'color', [0.85 0.33 0.10]);
end
polarplot(pax, wg.phi, p/max(arcP), '.', 'markersize', 4, 'color', [0 0.3 0.7]);
hold(pax, 'off')
pax.RLim = [0 1.05]; pax.FontSize = 14; pax.RTickLabel = {};
title(pax, sprintf('guided power per direction, one free-space mode (%d arcs)', numel(edges)))
print(fig, fullfile(cfg.figdir, 'panel_c_directions.png'), '-dpng', '-r150'); close(fig)
end

function panel_scaling(cfg, cells, g)
fig = figure('visible', 'off', 'position', [0 0 1200 600]);
colors = [0.9 0 0; 1 0.5 0; 1 0.9 0];
hd = gobjects(1, numel(cfg.NA_list));
hc = gobjects(1, numel(cfg.NA_list));
hold on
for i = 1:numel(cfg.NA_list)
    NA = cfg.NA_list(i);
    L = []; N = [];
    for ic = 1:numel(cells)
        if cells(ic).NA == NA && exist(cell_file(cfg, cells(ic)), 'file')
            d = load(cell_file(cfg, cells(ic)), 'Nmodes');
            L(end+1) = cells(ic).L; N(end+1) = d.Nmodes; %#ok<AGROW>
        end
    end
    Lc = cfg.L_list;
    Nc = arrayfun(@(l) phase_matched_estimate(NA, l, cfg.neff, g.beta), Lc);
    hc(i) = plot(Lc, Nc, '--', 'color', colors(i,:), 'linewidth', 2);
    if ~isempty(L)
        [L, o] = sort(L); N = N(o);
        hd(i) = plot(L, N, '-o', 'color', colors(i,:), 'linewidth', 3, ...
                     'markerfacecolor', colors(i,:), 'markersize', 10);
    end
end
hold off; grid on; box on
xlabel('L^{coupler}/\lambda'); ylabel('N^{coupled}')
ok = isgraphics(hd);
legend([hd(ok), hc(1)], [arrayfun(@(a) sprintf('NA = %.1f', a), cfg.NA_list(ok), ...
       'UniformOutput', false), {'phase-matched estimate'}], 'location', 'northwest')
set(gca, 'fontsize', 16)
print(fig, fullfile(cfg.figdir, 'panel_e_scaling.png'), '-dpng', '-r150'); close(fig)
end

% ------------------------------------------------------------------- tests

function run_tests(cfg)
g = geometry(cfg);
fprintf('n_core = %.6f   beta/k0 = %.6f   chimax = %.6f\n', ...
        g.n_core, g.beta/cfg.k0, g.chimax);

pass = @(name, ok, val) fprintf('%-34s %-4s %s\n', name, verdict(ok), val);

% guided mode: dispersion residual and unit power
a  = cfg.t_wg/2;
zq = linspace(-a-6, a+6, 20000);
P  = 0.5*(g.beta/cfg.k0)*trapz(zq, g.Efun(zq).^2);
pass('T_slab  effective index', abs(g.beta/cfg.k0 - cfg.neff) < 1e-9, ...
     sprintf('neff = %.10f', g.beta/cfg.k0));
pass('T_slab  unit power', abs(P-1) < 1e-6, sprintf('P = %.10f', P));

% slab reflection against a direct transfer-matrix evaluation
kpar = (0:0.2:0.9).'*cfg.k0;
for pol = 'sp'
    r1 = slab_reflection(cfg, g, kpar, pol);
    r2 = slab_reflection_tmm(cfg, g, kpar, pol);
    pass(sprintf('T_fresnel  %c vs transfer matrix', pol), ...
         max(abs(r1-r2)) < 1e-10, sprintf('max diff = %.2e', max(abs(r1-r2))));
    pass(sprintf('T_fresnel  %c lossless |r|<=1', pol), ...
         max(abs(r1)) <= 1+1e-12, sprintf('max |r| = %.6f', max(abs(r1))));
end

% gradient against per-degree-of-freedom central differences
c  = struct('NA', 0.1, 'L', 25, 'nfour', n_fourier(0.1, cfg.neff)/2);
fs = free_space_basis(cfg, g, c.NA, c.L);
th = (0:c.nfour-1)*pi/c.nfour;
wg = waveguide_basis(cfg, g, c.NA, c.L, th);
rng(0);
p0 = [0.8+0.1*rand(1,c.nfour), 0.4*cfg.k0*c.NA*ones(1,c.nfour)];
[~, ga] = fom_and_grad(p0, cfg, g, fs, wg, th, c.L);
gf = zeros(size(ga));
h  = 1e-6;
for j = 1:numel(p0)
    pp = p0; pp(j) = pp(j) + h*max(1, abs(p0(j)));
    pm = p0; pm(j) = pm(j) - h*max(1, abs(p0(j)));
    gf(j) = (fom_and_grad(pp, cfg, g, fs, wg, th, c.L) - ...
             fom_and_grad(pm, cfg, g, fs, wg, th, c.L))/(2*h*max(1, abs(p0(j))));
end
rel = max(abs(ga-gf)./max(abs(gf), eps));
pass('T_grad  analytic vs central FD', rel < 1e-6, sprintf('max rel err = %.2e', rel));

% power bound on every free-space column
S  = smatrix(cfg, g, fs, wg, th, p0(1:c.nfour), p0(c.nfour+1:end), c.L);
pb = max(sum(abs(S).^2, 1));
fprintf('%-34s %-4s max column power = %.4e%s\n', 'T_born  first-Born validity', 'INFO', pb, ...
        born_note(pb));

% convergence of the mode count in the guided-azimuth sampling
base = cfg.oversample;
Nm = zeros(1,3);
for i = 1:3
    cfg2 = cfg; cfg2.oversample = base*2^(i-1);
    wg2 = waveguide_basis(cfg2, g, c.NA, c.L, th);
    [~, ~, CE] = fom_and_grad(p0, cfg2, g, fs, wg2, th, c.L);
    Nm(i) = sum(CE >= cfg.CE_cutoff*max(CE));
end
pass('T_conv  count vs oversampling', abs(Nm(3)-Nm(2)) <= 1, ...
     sprintf('%d, %d, %d at %dx, %dx, %dx', Nm, base, 2*base, 4*base));

% rank against the only true bound, and against the phase-matched estimate
[~, ~, CE] = fom_and_grad(p0, cfg, g, fs, wg, th, c.L);
nr = sum(CE > 1e-12*max(CE));
nb = min(numel(wg.phi), numel(fs.kx));
np = phase_matched_estimate(c.NA, c.L, cfg.neff, g.beta);
pass('T_rank  rank <= min(Nwg,Nfree)', nr <= nb, sprintf('rank %d, bound %d', nr, nb));
fprintf('%-34s %-4s rank/phase-matched estimate = %.2f (may exceed 1)\n', ...
        'T_rank  vs phase-matched estimate', 'INFO', nr/np);
end

function s = born_note(pb)
if pb > 1
    s = '  (> 1: outside the first-Born regime at these parameters)';
else
    s = '';
end
end

function s = verdict(ok)
if ok, s = 'PASS'; else, s = 'FAIL'; end
end

function r = slab_reflection(cfg, g, kpar, pol)
k0  = cfg.k0;
kz0 = sqrt(k0^2 - kpar.^2);
kz1 = sqrt(g.n_core^2*k0^2 - kpar.^2);
if pol == 's', Y0 = kz0; Y1 = kz1; else, Y0 = kz0; Y1 = kz1/g.n_core^2; end
r01 = (Y0 - Y1)./(Y0 + Y1);
ph  = exp(2i*kz1*cfg.t_wg);
r   = r01.*(1 - ph)./(1 - r01.^2.*ph);
end

function r = slab_reflection_tmm(cfg, g, kpar, pol)
% Independent evaluation by explicit 2x2 interface and propagation matrices.
k0  = cfg.k0;
kz0 = sqrt(k0^2 - kpar.^2);
kz1 = sqrt(g.n_core^2*k0^2 - kpar.^2);
if pol == 's', Y0 = kz0; Y1 = kz1; else, Y0 = kz0; Y1 = kz1/g.n_core^2; end
r = zeros(size(kpar));
for i = 1:numel(kpar)
    I01 = 0.5*[1 + Y1(i)/Y0(i), 1 - Y1(i)/Y0(i); 1 - Y1(i)/Y0(i), 1 + Y1(i)/Y0(i)];
    Pl  = [exp(-1i*kz1(i)*cfg.t_wg), 0; 0, exp(1i*kz1(i)*cfg.t_wg)];
    I10 = 0.5*[1 + Y0(i)/Y1(i), 1 - Y0(i)/Y1(i); 1 - Y0(i)/Y1(i), 1 + Y0(i)/Y1(i)];
    M   = I01*Pl*I10;
    r(i) = M(2,1)/M(1,1);
end
end

function etch_scan(cfg)
% How the attainable coupling efficiency depends on the depth of the
% patterned layer, at the reference cell.  Writes etch_scan.mat.
h_list = [0.03, 0.05, 0.075, 0.10, 0.15, 0.20];
NA = cfg.ref_NA;  L = cfg.ref_L;
nfour = n_fourier(NA, cfg.neff)/2;
thetas = (0:nfour-1)*pi/nfour;

CEmax = zeros(size(h_list));  N2 = zeros(size(h_list));
pb = zeros(size(h_list));     Nr = zeros(size(h_list));
for i = 1:numel(h_list)
    c2 = cfg;  c2.h_etch = h_list(i);
    g  = geometry(c2);
    fs = free_space_basis(c2, g, NA, L);
    wg = waveguide_basis(c2, g, NA, L, thetas);
    lb = [-ones(1,nfour), zeros(1,nfour)];
    ub = [ ones(1,nfour), c2.k0*NA*ones(1,nfour)];
    opts = optimoptions('fmincon','Algorithm','sqp','SpecifyObjectiveGradient',true, ...
                        'MaxIterations',cfg.maxiter,'Display','off');
    rng(7);
    p0 = [0.8+0.1*rand(1,nfour), 0.4*c2.k0*NA*ones(1,nfour)];
    p  = fmincon(@(p) fom_and_grad(p,c2,g,fs,wg,thetas,L), p0, [],[],[],[], lb, ub, [], opts);
    [~,~,CE] = fom_and_grad(p, c2, g, fs, wg, thetas, L);
    S  = smatrix(c2, g, fs, wg, thetas, p(1:nfour), p(nfour+1:end), L);
    CEmax(i) = max(CE);  N2(i) = sum(CE >= cfg.CE_cutoff);
    Nr(i) = sum(CE >= 1e-6*max(CE));  pb(i) = max(sum(abs(S).^2,1));
    fprintf('h_etch = %.3f : CEmax %.4f, %4d above %.0f%%, rank %4d, max column power %.4f%s\n', ...
            h_list(i), CEmax(i), N2(i), 100*cfg.CE_cutoff, Nr(i), pb(i), born_note(pb(i)));
end
save('etch_scan.mat','h_list','CEmax','N2','Nr','pb','NA','L','nfour');
end

function paper_panels(cfg, cells)
out = 'paper_panels';
if ~exist(out, 'dir'), mkdir(out); end
g = geometry(cfg);
ref = load(cell_file(cfg, cells(arrayfun(@(c) c.NA==cfg.ref_NA && c.L==cfg.ref_L, cells))));

panel_c_paper(cfg, ref, out);
panel_d_paper(cfg, g, ref, out);
panel_e_paper(cfg, cells, out);
panel_f_fig5_paper(cfg, cells, out);
end

function export_exact(fig, file, W, H)
% export_fig gives a genuinely transparent PNG; the slide background is dark,
% so anything opaque would show as a white box.
addpath(genpath('~/main_projects/ASML/multifunctional_gc_github/export_fig'));
set(fig, 'InvertHardcopy', 'off', 'Color', 'w');
for ax = findall(fig, 'type', 'axes').'
    set(ax, 'XColor', [0.15 0.15 0.15], 'YColor', [0.15 0.15 0.15], 'Color', 'w');
end
tmp = [tempname '.png'];
export_fig(tmp, '-transparent', '-m2', fig);
[im, ~, al] = imread(tmp);
im = imresize(im, [H W], 'lanczos3');
if isempty(al)
    imwrite(im, file);
else
    imwrite(im, file, 'Alpha', imresize(al, [H W], 'lanczos3'));
end
delete(tmp); close(fig);
fprintf('wrote %s  (%dx%d)\n', file, W, H);
end

function panel_c_paper(cfg, ref, out)
% replaces image23.png (922 x 736): coupling efficiency vs mode number
CE = ref.CE;
xlims = [0, 250];          % same extent as the panel it replaces, so the
ylims = [0, 0.12];         % slide's own tick labels stay in register
fig = figure('units','pixels','position',[100 100 1400 1118],'visible','off');
scatter(1:numel(CE), CE, 220, [0.8500 0.3250 0.0980], 'filled');
xlim(xlims); ylim(ylims); box on
set(gca,'fontsize',42,'linewidth',1.5)
xticks(0:50:250); yticks(0:0.06:0.12)
set(gca,'XTickLabel',[],'YTickLabel',[])
daspect([0.5*(xlims(2)-xlims(1)), ylims(2)-ylims(1), 1])
export_exact(fig, fullfile(out,'panel_c_image23_replacement.png'), 922, 736);
end

function panel_d_paper(cfg, g, ref, out)
% replaces image43.png (799 x 739): guided-from-free-space scattering matrix
% over the whole light cone, columns and rows split desired / undesired.
nf = ref.cell.nfour;
[fs, ndes_f] = light_cone_basis(cfg, g, ref.cell.NA, ref.cell.L);
[wg, ndes_w] = all_azimuth_basis(cfg, g, ref.cell.NA, ref.cell.L, ref.thetas);
S = smatrix(cfg, g, fs, wg, ref.thetas, ref.params(1:nf), ref.params(nf+1:end), ref.cell.L);
A = abs(S).^2;

nb = 700;                                   % display grid
ri = round(linspace(1, size(A,1), nb));
ci = round(linspace(1, size(A,2), nb));
Ad = A(ri, ci);
xd = interp1(1:size(A,2), 1:size(A,2), ci);
fdes = interp1(ci, 1:nb, ndes_f, 'linear', 'extrap');
wdes = interp1(ri, 1:nb, ndes_w, 'linear', 'extrap');

fig = figure('units','pixels','position',[100 100 1200 1110],'visible','off');
imagesc(Ad); axis xy
colormap(white_to_red(256));
cmax = 0.01*max(Ad(:));   % desired blocks saturate; undesired stay white if suppressed below 1%
clim([0, cmax]);
cb = colorbar; cb.LineWidth = 1.5;
set(gca,'fontsize',34,'linewidth',1.5,'xtick',[],'ytick',[])
daspect([1 1 1])
fprintf('   panel d colour scale 0 .. %.3e (99.9th percentile), max %.3e\n', cmax, max(Ad(:)));
export_exact(fig, fullfile(out,'panel_d_image43_replacement.png'), 799, 739);
fprintf('   panel d: %d free-space modes (%d desired), %d guided (%d desired)\n', ...
        size(A,2), ndes_f, size(A,1), ndes_w);
end

function panel_e_paper(cfg, cells, out)
% replaces image102.png (1068 x 548): mode count vs coupler size, three NA,
% scatter with linear fits.  Style copied from Nrank_vs_Rfiber_vs_NA_plot.m.
colors = [0.9 0 0.0; 1.0 0.5 0; 1.0 0.9 0.0];
fig = figure('units','pixels','position',[100 100 1920 985],'visible','off');
hold on
ymax = 0;
for i = 1:numel(cfg.NA_list)
    NA = cfg.NA_list(i);
    L = []; N = [];
    for ic = 1:numel(cells)
        if cells(ic).NA == NA && exist(cell_file(cfg, cells(ic)), 'file')
            d = load(cell_file(cfg, cells(ic)), 'Nmodes');
            L(end+1) = cells(ic).L; N(end+1) = d.Nmodes; %#ok<AGROW>
        end
    end
    if isempty(L), continue; end
    [L, o] = sort(L); N = N(o);
    scatter(L, N, 800, colors(i,:), 'filled');
    P = polyfit(L, N, 1);
    plot(L, P(1)*L + P(2), 'color', colors(i,:), 'linewidth', 10);
    ymax = max(ymax, max(N));
    fprintf('   NA = %.1f : slope %.2f modes per lambda, intercept %.0f\n', NA, P(1), P(2));
end
hold off
xlims = [20, 104];
ylims = [0, 100*ceil(1.06*ymax/100)];
xlim(xlims); ylim(ylims); box on
set(gca,'fontsize',42,'linewidth',1.5)
xticks(20:20:100); yticks(0:500:ylims(2))
daspect([0.5*(xlims(2)-xlims(1)), ylims(2)-ylims(1), 1])
export_exact(fig, fullfile(out,'panel_e_image102_replacement.png'), 1068, 548);
end

function panel_f_fig5_paper(cfg, cells, out)
% replaces image26.png (1860 x 794) in the main-text figure: same quantity as
% the supplementary panel but at that figure's aspect ratio and font size.
% Style copied from Nrank_vs_Rfiber_vs_NA_plot_fig5d_AR2p45.m.
target_AR = 1556/635.5;
colors = [0.9 0 0.0; 1.0 0.5 0; 1.0 0.9 0.0];
fig = figure('units','pixels','position',[100 100 2200 900],'visible','off');
hold on
ymax = 0;
for i = 1:numel(cfg.NA_list)
    NA = cfg.NA_list(i);
    L = []; N = [];
    for ic = 1:numel(cells)
        if cells(ic).NA == NA && exist(cell_file(cfg, cells(ic)), 'file')
            d = load(cell_file(cfg, cells(ic)), 'Nmodes');
            L(end+1) = cells(ic).L; N(end+1) = d.Nmodes; %#ok<AGROW>
        end
    end
    if isempty(L), continue; end
    [L, o] = sort(L); N = N(o);
    scatter(L, N, 800, colors(i,:), 'filled');
    P = polyfit(L, N, 1);
    plot(L, P(1)*L + P(2), 'color', colors(i,:), 'linewidth', 10);
    ymax = max(ymax, max(N));
end
hold off
xlims = [20, 104];
ylims = [0, 100*ceil(1.06*ymax/100)];
xlim(xlims); ylim(ylims); box on
set(gca,'fontsize',59,'linewidth',1.5)
set(gca,'XTick',20:20:100); set(gca,'YTick',0:500:ylims(2));
daspect([(xlims(2)-xlims(1))/target_AR, ylims(2)-ylims(1), 1])
export_exact(fig, fullfile(out,'panel_f_image26_replacement.png'), 1860, 794);
end

function [fs, ndes] = light_cone_basis(cfg, g, NA, L)
% Free-space modes over the whole light cone, desired (inside NA) first.
dq = 2*pi/L;
kmax = cfg.k0*sqrt(1 - cfg.kz_min^2);
n = ceil(kmax/dq);
[QX, QY] = meshgrid((-n:n)*dq, (-n:n)*dq);
r2 = QX.^2 + QY.^2;
keep = r2 <= kmax^2;
qx = QX(keep).'; qy = QY(keep).'; rr = r2(keep).';
des = rr <= (cfg.k0*NA)^2;
qx = [qx(des), qx(~des)];  qy = [qy(des), qy(~des)];
ndes = 2*nnz(des);
fs = build_fs(cfg, g, qx, qy, L);
end

function fs = build_fs(cfg, g, qx, qy, L)
fs.kx = [qx, qx];  fs.ky = [qy, qy];
kpar = hypot(qx, qy).';  kpar(kpar < 1e-12) = 1e-12;
khx = (qx.')./kpar;  khy = (qy.')./kpar;
[As, ~] = incident_profile(cfg, g, kpar, 's');
[Ap, ~] = incident_profile(cfg, g, kpar, 'p');
Iz_s = trapz(g.z, As .* g.Emode, 2);
Iz_p = trapz(g.z, Ap .* g.Emode, 2);
fs.ex = [(-khy).*Iz_s;  khx.*Iz_p].'/L;
fs.ey = [( khx).*Iz_s;  khy.*Iz_p].'/L;
[fs.kxu, ~, ix] = unique(fs.kx);
[fs.kyu, ~, iy] = unique(fs.ky);
fs.ix = ix(:).';  fs.iy = iy(:).';
end

function [wg, ndes] = all_azimuth_basis(cfg, g, NA, L, thetas)
% Every guided azimuth, those inside the phase-matched arcs listed first.
dphi = 2*pi/(g.beta*L);
phi = 0:dphi:(2*pi - dphi);
half = asin(min(NA/cfg.neff, 1));
des = false(size(phi));
for th = [thetas, thetas + pi]
    des = des | abs(mod(phi - th + pi, 2*pi) - pi) <= half;
end
phi = [phi(des), phi(~des)];
ndes = nnz(des);
wg.phi = phi;  wg.dphi = dphi;
wg.tx = -sin(phi);  wg.ty = cos(phi);
end

function c = white_to_red(n)
c = [ones(n,1), linspace(1,0,n).', linspace(1,0,n).'];
end
