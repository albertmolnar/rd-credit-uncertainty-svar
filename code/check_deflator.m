%% Does deflating R&D by the (sign-restricted) GDP deflator drive the result?
% Re-estimates the baseline with log NOMINAL R&D. If real R&D = nominal - p, the
% restriction p < 0 mechanically pushes measured real R&D up, so the nominal
% response bounds how much of the real response is mechanical.
clear; clc;
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, 'lib'));
cfg = default_config(here);
rng(cfg.seed);
D = build_dataset(cfg);
c = setf(cfg, 'vars', {'y', 'p', 'rdnom', 'r', 'ebp', 'epu'}, 'roleRD', 'rdnom', 'nDraws', 500);
o = run_spec(D, c, false);
ip = find(strcmp(c.vars, 'p'));
ird = find(strcmp(c.vars, 'rdnom'));
for s = 1:2
    nm = o.R.shockNames{s};
    qn = pctl(o.R.irf(:, 9, ird, s), [16 50 84]);
    qr = pctl(o.R.irf(:, 9, ird, s) - o.R.irf(:, 9, ip, s), [16 50 84]);
    qp = pctl(o.R.irf(:, 9, ip, s), [16 50 84]);
    fprintf('%-14s h=8: nominal R&D %.2f [%.2f, %.2f], P(<0)=%.2f | prices %.2f | nominal-p %.2f [%.2f, %.2f]\n', ...
        nm, qn(2), qn(1), qn(3), mean(o.R.irf(:, 9, ird, s) < 0), qp(2), qr(2), qr(1), qr(3));
end
