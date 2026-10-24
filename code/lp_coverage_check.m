%% Calibrated Monte Carlo: coverage and size-adjusted inference for the LPs
% Simulates artificial samples from the estimated baseline SVAR (median-target
% draw), feeds the true structural shocks to local_projection.m, and records
%  (i)  how often nominal 68/90/95% intervals contain the true response, and
%  (ii) the simulated distribution of |t| = |b - b_true|/se, whose quantiles
%       serve as size-adjusted critical values (a parametric bootstrap-t).
% The bootstrap p-value of the actual medium-run LP estimate is then the share
% of simulated |t| at least as large as the actual |b/se|. Report both in the
% paper. Run main.m first (needs output/baseline.mat).

clear; clc;
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, 'lib'));
cfg = default_config(here);
load(fullfile(cfg.outDir, 'baseline.mat'), 'base');
warning('off', 'MATLAB:rankDeficientMatrix');
rng(cfg.seed + 1);

M = base.M; R = base.R;
d = median_target(R.irf);
B = R.B(:, :, d);
A = R.A(:, :, d);
n = M.n; p = M.p; T = M.T + p;
ird = find(strcmp(M.vars, cfg.roleRD));
w = cfg.lrWindow;
hs = [4 12 24];
nrep = 1000;
lab = [arrayfun(@(h) sprintf('$h=%d$', h), hs, 'UniformOutput', false), {sprintf('avg. $h=%d$--$%d$', w)}];

body = cell(0, 7);
shocks = find(ismember(R.shockNames, {'Credit supply', 'Uncertainty'}));
for s = shocks
    irfTrue = var_irf(B, n, p, M.ndet, A(:, s), w(2));
    truthH = irfTrue(:, ird);
    tru = [truthH(hs+1); mean(truthH(w(1)+1:w(2)+1))];
    absT = zeros(nrep, numel(tru));
    hit = zeros(nrep, numel(tru), 3);
    for r = 1:nrep
        E = randn(T, n);
        Y = zeros(T, n);
        Y(1:p, :) = M.Y(1:p, :);
        for t = p+1:T
            x = 1;
            for l = 1:p, x = [x, Y(t-l, :)]; end %#ok<AGROW>
            Y(t, :) = x*B([1, M.ndet+1:end], :) + E(t, :)*A';
        end
        L = local_projection(Y(:, ird), E(:, s), Y, ones(T, 1), w(2), cfg.lpL, w);
        est = [L.beta(hs+1); L.lr.beta];
        se = [L.se(hs+1); L.lr.se];
        nu = [L.df(hs+1); L.lr.df];
        absT(r, :) = abs(est - tru)./se;
        lv = [0.84 0.95 0.975];
        for j = 1:3
            hit(r, :, j) = absT(r, :)' <= arrayfun(@(v) tq(lv(j), v), nu);
        end
    end
    cov = squeeze(mean(hit, 1));
    cAdj = pctl(absT, 95);

    % Bootstrap p-value of the actual medium-run estimate.
    f = matlab.lang.makeValidName([R.shockNames{s} '_' cfg.roleRD]);
    Lact = base.lp.(f);
    tAct = abs(Lact.lr.beta/Lact.lr.se);
    pBoot = mean(absT(:, end) >= tAct);

    fprintf('\n%s shock -> R&D (nominal 68/90/95%% coverage, size-adjusted 5%% critical value)\n', R.shockNames{s});
    for i = 1:numel(lab)
        fprintf('  %-22s %.2f  %.2f  %.2f   c* = %.2f\n', lab{i}, cov(i, :), cAdj(i));
        pcol = '';
        if i == numel(lab), pcol = sprintf('%.3f', pBoot); end
        body(end+1, :) = {R.shockNames{s}, lab{i}, sprintf('%.2f', cov(i, 1)), sprintf('%.2f', cov(i, 2)), ...
            sprintf('%.2f', cov(i, 3)), sprintf('%.2f', cAdj(i)), pcol}; %#ok<AGROW>
    end
    fprintf('  actual medium-run estimate %.2f (se %.2f): bootstrap p-value %.3f\n', Lact.lr.beta, Lact.lr.se, pBoot);
end
write_tex_table(fullfile(cfg.tabDir, 'tab_lp_coverage.tex'), ...
    {'Shock', 'Horizon', 'Cov. 68\%', 'Cov. 90\%', 'Cov. 95\%', '$c^{*}_{5\%}$', 'Boot. $p$'}, body, 'llccccc');

function q = tq(p, nu)
lo = 0; hi = 50;
for it = 1:60
    mid = (lo + hi)/2;
    if 1 - 0.5*betainc(nu/(nu + mid^2), nu/2, 0.5) < p, lo = mid; else, hi = mid; end
end
q = (lo + hi)/2;
end
