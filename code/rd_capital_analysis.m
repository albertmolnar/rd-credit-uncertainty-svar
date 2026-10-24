%% R&D capital stock following Sliker (2007) -- construction and diagnostics
% 1. Annual real R&D investment from the quarterly NIPA series.
% 2. Benchmark regression (eq. 4) with the R^2 subsample rule, implied g_I, g_K.
% 3. Perpetual inventory (eq. 1, mid-year convention) for several depreciation
%    rates and for the accelerating-depreciation variant (eqs. 5-9).
% 4. How much the benchmark and delta matter over the VAR sample.
% 5. Annual vs quarterly stocks, and how smooth (interpolated) quarterly R&D is.
% Writes tables to output/tables and figures to output/figures; run after
% import/build (main.m sections 1-2), or on its own.

clear; close all; clc;
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, 'lib'));
cfg = default_config(here);
for d = {cfg.rawDir, cfg.outDir, cfg.figDir, cfg.tabDir}
    if ~exist(d{1}, 'dir'), mkdir(d{1}); end
end
import_local_data(cfg);
D = build_dataset(cfg);

%% 1. Annual real investment (calendar-year sums of quarterly flows)
ok = ~isnan(D.rd);
q = D.qidx(ok);
Iq = exp(D.rd(ok)/100)/4;                  % real quarterly flow
yr = floor(q/4);
ya = unique(yr);
nq = accumarray(yr - ya(1) + 1, 1);
Ia = accumarray(yr - ya(1) + 1, Iq);
full = nq == 4;
ya = ya(full); Ia = Ia(full);
fprintf('Annual real R&D investment %d-%d (%d years)\n', ya(1), ya(end), numel(ya));

%% 2-3. Benchmarks and stocks
deltas = [0.10 0.15 0.20 0.25];
nSliker = 8.809e-4;                         % Sliker (2007) slope of delta_t, all other nonfin. corp. business
variants = {};
for dlt = deltas
    variants(end+1, :) = {sprintf('$\\delta=%.2f$', dlt), dlt*ones(size(ya)), 'r2'}; %#ok<SAGROW>
end
variants(end+1, :) = {'$\delta=0.15$, full-sample $g$', 0.15*ones(size(ya)), 'full'};
variants(end+1, :) = {'$\delta=0.15$, zero benchmark', 0.15*ones(size(ya)), 'zero'};
variants(end+1, :) = {'$\delta_t=0.15+n(t-2000)$', 0.15 + nSliker*(ya - 2000), 'r2'};

in = ya >= 1986 & ya <= 2019;
body = cell(size(variants, 1), 9);
K = zeros(numel(ya), size(variants, 1));
for v = 1:size(variants, 1)
    S = rd_capital_sliker(Ia, ya, variants{v, 2}, 'benchmark', variants{v, 3});
    K(:, v) = S.K;
    gk = mean(diff(log(S.K(in))));
    sh86 = S.benchShare(ya == 1986);
    body(v, :) = {variants{v, 1}, sprintf('%d--%d', S.subsample), sprintf('%.4f (%.4f)', S.m, S.m_se), ...
        sprintf('%.4f', S.R2), sprintf('%.4f', S.gI), sprintf('%.4f', S.gK), sprintf('%.1f', S.K0), ...
        sprintf('%.2e', sh86), sprintf('%.2f', 100*gk)};
    if v == 2, base = S; end
end
write_tex_table(fullfile(cfg.tabDir, 'tab_rdk_sliker.tex'), ...
    {'Variant', 'Subsample', '$m$ (s.e.)', '$R^2$', '$g_I$', '$g_K$', '$K_{1946}$', 'Bench. share 1986', '$\Delta\ln K$ 1986--2019 (\%)'}, ...
    body, 'llccccccc');
disp(cell2table(body));

% Correlation of quarterly-relevant dynamics across variants: growth rates 1986-2019
G = diff(log(K(in, :)));
C = corrcoef(G);
fprintf('Correlation of annual growth of K across variants (1986-2019): min %.3f\n', min(C(:)));

%% 4. Annual vs quarterly stocks (quarterly stock as built in build_dataset)
Kq = exp(D.rdk/100);
qq = D.qidx;
q4 = mod(qq, 4) == 3 & ~isnan(Kq);
[tf, loc] = ismember(floor(qq(q4)/4), ya);
ratio = Kq(q4)./K(loc(tf), 2);
ratio = ratio(tf);
fprintf('Quarterly Q4 stock / annual end-of-year stock: mean %.4f, range [%.4f, %.4f]\n', ...
    mean(ratio), min(ratio), max(ratio));

%% 5. Smoothness of quarterly NIPA R&D (sign of interpolation from annual data)
s = D.qidx >= 1986*4 & D.qidx <= 2019*4 + 3;
gr = diff(D.rd(s));
gy = diff(D.y(s));
yrq = floor(D.qidx(s)/4); yrq = yrq(2:end);
ac = @(x, l) autocorr1(x, l);
within = @(x) within_share(x, yrq);
rows = {
    'Real R\&D investment', sprintf('%.2f', ac(gr, 1)), sprintf('%.2f', ac(gr, 2)), sprintf('%.2f', ac(gr, 4)), sprintf('%.2f', std(gr)), sprintf('%.2f', within(gr))
    'Real GDP',             sprintf('%.2f', ac(gy, 1)), sprintf('%.2f', ac(gy, 2)), sprintf('%.2f', ac(gy, 4)), sprintf('%.2f', std(gy)), sprintf('%.2f', within(gy))
    };
write_tex_table(fullfile(cfg.tabDir, 'tab_rd_smoothness.tex'), ...
    {'Quarterly growth of', 'AC(1)', 'AC(2)', 'AC(4)', 'S.d. (pp)', 'Within-year var. share'}, rows, 'lccccc');
disp(cell2table(rows));

%% Figures
fig = figure('Units', 'inches', 'Position', [1 1 7 2.4], 'Color', 'w');
subplot(1, 3, 1);
plot(ya(3:end), base.r2(3:end), 'k-', 'LineWidth', 1.2); hold on;
xline(base.subsample(2), 'r--');
title('R^2 of eq. (4) by subsample end', 'FontWeight', 'normal', 'FontSize', 8);
set(gca, 'FontSize', 8);
subplot(1, 3, 2);
plot(ya, log(K(:, 1:4)), 'LineWidth', 1.1);
legend(arrayfun(@(x) sprintf('\\delta = %.2f', x), deltas, 'UniformOutput', false), 'Location', 'southeast', 'FontSize', 6, 'Box', 'off');
title('log real R&D capital stock', 'FontWeight', 'normal', 'FontSize', 8);
set(gca, 'FontSize', 8);
subplot(1, 3, 3);
semilogy(ya, base.benchShare, 'k-', 'LineWidth', 1.2);
xline(1986, 'r--');
title('Benchmark share of stock (\delta = 0.15)', 'FontWeight', 'normal', 'FontSize', 8);
set(gca, 'FontSize', 8);
exportgraphics(fig, fullfile(cfg.figDir, 'fig_rdk.pdf'), 'ContentType', 'vector');
close(fig);

T = array2table([ya Ia K], 'VariableNames', ['year', 'I_real', matlab.lang.makeValidName(strrep(variants(:, 1)', '$', ''))]);
writetable(T, fullfile(cfg.outDir, 'rd_capital_annual.csv'));
writetable(D(:, {'quarter', 'rd', 'rdk'}), fullfile(cfg.outDir, 'rd_capital_quarterly.csv'));

function r = autocorr1(x, l)
c = corrcoef(x(l+1:end), x(1:end-l));
r = c(1, 2);
end

function s = within_share(x, yr)
% Share of the variance of quarterly growth that is within calendar years.
g = yr - yr(1) + 1;
m = accumarray(g, x, [], @mean);
s = 1 - var(m(g))/var(x);
end
