%% Compact figures and tables for the 6-page conference version (run after main.m)
clear; close all; clc;
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, 'lib'));
cfg = default_config(here);
load(fullfile(cfg.outDir, 'baseline.mat'), 'base');
load(fullfile(cfg.outDir, 'robustness.mat'), 'results', 'specs');
M = base.M; R = base.R;
ird = find(strcmp(M.vars, 'rd'));

% Local projections of R&D, one column wide
fig = figure('Units', 'inches', 'Position', [1 1 3.5 1.7], 'Color', 'w');
nm = {'CreditSupply_rd', 'Uncertainty_rd'};
tt = {'Credit supply shock', 'Uncertainty shock'};
for a = 1:2
    L = base.lp.(nm{a});
    ax = subplot(1, 2, a);
    plot_band(ax, L.h, L.beta - L.c90.*L.se, L.beta - L.c68.*L.se, L.beta, ...
        L.beta + L.c68.*L.se, L.beta + L.c90.*L.se, tt{a});
    set(ax, 'FontSize', 7);
    xlabel(ax, 'Quarters', 'FontSize', 7);
    if a == 1, ylabel(ax, 'Percent', 'FontSize', 7); end
end
exportgraphics(fig, fullfile(cfg.figDir, 'fig_lp_rd.pdf'), 'ContentType', 'vector');
close(fig);

% Historical decomposition, one column wide
d = median_target(R.irf);
[contrib, baseline] = hist_decomp(M, R.B(:, :, d), R.A(:, :, d));
yr = floor(M.qidx/4) + mod(M.qidx, 4)/4;
fig = figure('Units', 'inches', 'Position', [1 1 3.5 2.3], 'Color', 'w');
ax = axes(fig); hold(ax, 'on');
plot(ax, yr, M.Y(:, ird) - baseline(:, ird), 'k-', 'LineWidth', 1.4);
cols = lines(4);
for s = 1:3
    plot(ax, yr, contrib(:, ird, s), '-', 'Color', cols(s, :), 'LineWidth', 1.0);
end
plot(ax, yr, sum(contrib(:, ird, 4:end), 3), ':', 'Color', cols(4, :), 'LineWidth', 1.0);
yline(ax, 0, 'k-');
legend(ax, {'Total', 'Credit supply', 'Uncertainty', 'Aggr. supply', 'Other'}, ...
    'Location', 'southoutside', 'Orientation', 'horizontal', 'NumColumns', 3, 'FontSize', 6, 'Box', 'off');
ylabel(ax, 'Percent', 'FontSize', 7);
set(ax, 'FontSize', 7, 'Box', 'on');
exportgraphics(fig, fullfile(cfg.figDir, 'fig_hd_rd_small.pdf'), 'ContentType', 'vector');
close(fig);

% Short robustness table
short = {'Baseline', 'Baa--10y spread for EBP', 'VIX for EPU (from 1990)', 'R\&D/GDP for log R\&D', ...
    'Sample from 1994', 'Sample to 2026Q2, COVID dummies', 'Narrative restrictions', ...
    'Restrictions at $h=0$ only', 'Restrictions at $h=0,\dots,3$', 'Two lags', 'Six lags', ...
    'Credit: yield unrestricted', 'Uncertainty: EPU sign only', 'Recursive, financial last', ...
    'Recursive, financial first'};
body = cell(0, 5);
for i = 1:numel(results)
    if isempty(results{i}), continue; end
    S = results{i}.summary;
    c = @(f) sprintf('%.2f [%.2f, %.2f]', S.(f).h8(2), S.(f).h8(1), S.(f).h8(3));
    body(end+1, :) = {short{i}, c('CreditSupply'), sprintf('%.2f', S.CreditSupply.pneg8), ...
        c('Uncertainty'), sprintf('%.2f', S.Uncertainty.pneg8)}; %#ok<SAGROW>
end
write_tex_table(fullfile(cfg.tabDir, 'tab_robustness_short.tex'), ...
    {'Specification', 'Credit, $h=8$', '$P(<0)$', 'Uncertainty, $h=8$', '$P(<0)$'}, body, 'lcccc');
% Unit-root table with readable names
U = readtable(fullfile(cfg.outDir, 'unitroot.csv'), 'TextType', 'string');
names = struct('y', 'Real GDP', 'p', 'GDP deflator', 'rd', 'Real R\&D investment', 'r', '1-year yield', ...
    'ebp', 'EBP', 'epu', 'EPU', 'baa', 'Baa--10y spread', 'vix', 'VIX');
body = cell(0, 5);
for v = ["y", "p", "rd", "r", "ebp", "epu", "baa", "vix"]
    t = U(U.variable == v, :);
    if isempty(t), continue; end
    body(end+1, :) = {names.(v), sprintf('%.2f (%.2f)', t.adf_level, t.adf_level_cv5), ...
        sprintf('%.2f (%.2f)', t.adf_diff, t.adf_diff_cv5), sprintf('%.2f (%.2f)', t.kpss_level, t.kpss_level_cv5), ...
        sprintf('%.2f', t.kpss_diff)}; %#ok<SAGROW>
end
write_tex_table(fullfile(cfg.tabDir, 'tab_unitroot_short.tex'), ...
    {'Variable', 'ADF level', 'ADF diff.', 'KPSS level', 'KPSS diff.'}, body, 'lcccc');

% Medium-run (h = 20..28) table: LP, VAR and implied capital stock
q = @(x) sprintf('%.2f [%.2f, %.2f]', pctl(x, 50), pctl(x, 16), pctl(x, 84));
Lc = base.lp.CreditSupply_rd; Lu = base.lp.Uncertainty_rd;
body = {
    'R\&D investment, LP', sprintf('%.2f (%.2f)', Lc.lr.beta, Lc.lr.se), sprintf('%.2f (%.2f)', Lu.lr.beta, Lu.lr.se)
    'R\&D investment, VAR', q(base.lrVAR(:, ird, 1)), q(base.lrVAR(:, ird, 2))
    'R\&D capital stock, VAR', q(base.lrK(:, 1)), q(base.lrK(:, 2))
    };
write_tex_table(fullfile(cfg.tabDir, 'tab_longrun_short.tex'), {'', 'Credit supply', 'Uncertainty'}, body, 'lcc');

for a = 1:2
    L = base.lp.(nm{a});
    [mn, hm] = min(L.beta(1:21));
    fprintf('LP %s: min %.2f at h=%d; beta(h=4,8,12,16,20,24) = %s\n', nm{a}, mn, hm - 1, ...
        sprintf('%.2f ', L.beta([5 9 13 17 21 25])));
end
fprintf('Compact outputs written.\n');
