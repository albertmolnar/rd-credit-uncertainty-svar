function out = run_spec(D, cfg, full)
% Estimates one specification. With full = true it also writes every figure and
% table used in the paper for this specification.

[Y, Z, qidx] = select_sample(D, cfg);
vars = cfg.vars;
n = numel(vars);
M = var_estimate(Y, cfg.p, Z);
M.qidx = qidx(cfg.p+1:end);
M.vars = vars;
if strcmp(cfg.ident, 'sign')
    S = make_restrictions(cfg, vars);
    R = sign_identify(M, S, cfg);
else
    R = chol_identify(M, cfg);
end
ird = find(strcmp(vars, cfg.roleRD));
out.M = M;
out.R = R;
out.cfg = cfg;
out.sample = sprintf('%s--%s', quarter_label(M.qidx(1)), quarter_label(M.qidx(end)));
out.summary = summarize_rd(R, ird);
if ~full, return; end

k = numel(R.shockNames);
H = cfg.H;
hh = 0:H;
fmt = @(x) sprintf('%.2f', x);

%% Lag-length criteria
LS = var_lagselect(Y, cfg.pmax, Z);
body = cell(cfg.pmax, 4);
for p = 1:cfg.pmax
    body{p, 1} = sprintf('%d', p);
    for c = 1:3
        s = sprintf('%.3f', LS.crit(p, c));
        if LS.best(c) == p, s = [s '$^{*}$']; end %#ok<AGROW>
        body{p, c+1} = s;
    end
end
write_tex_table(fullfile(cfg.tabDir, 'tab_lagselect.tex'), {'$p$', 'AIC', 'BIC', 'HQ'}, body, 'cccc');
fprintf('Lag selection: AIC %d, BIC %d, HQ %d (baseline uses p = %d)\n', LS.best, cfg.p);
out.lagselect = LS;

%% Residual and stability diagnostics
G = var_diagnostics(M, 4);
body = cell(0, 3);
body(end+1, :) = {'Largest companion root', fmt(G.maxroot), '--'};
for h = 1:4
    body(end+1, :) = {sprintf('LM autocorrelation, order %d', h), fmt(G.lm(h, 1)), fmt(G.lm(h, 3))}; %#ok<AGROW>
end
for i = 1:n
    body(end+1, :) = {sprintf('Jarque--Bera, %s eq.', vars{i}), fmt(G.jb(i, 1)), fmt(G.jb(i, 2))}; %#ok<AGROW>
end
for i = 1:n
    body(end+1, :) = {sprintf('ARCH(4) LM, %s eq.', vars{i}), fmt(G.arch(i, 1)), fmt(G.arch(i, 2))}; %#ok<AGROW>
end
write_tex_table(fullfile(cfg.tabDir, 'tab_diagnostics.tex'), {'Test', 'Statistic', '$p$-value'}, body, 'lcc');
out.diagnostics = G;

%% R&D responses (main figure) and all responses (appendix)
fig = figure('Units', 'inches', 'Position', [1 1 7 2.3], 'Color', 'w');
for s = 1:k
    ax = subplot(1, k, s);
    q = pctl(R.irf(:, :, ird, s), [5 16 50 84 95]);
    plot_band(ax, hh, q(1, :), q(2, :), q(3, :), q(4, :), q(5, :), R.shockNames{s});
    if s == 1, ylabel(ax, 'Real R&D investment (%)', 'FontSize', 8); end
end
exportgraphics(fig, fullfile(cfg.figDir, 'fig_irf_rd.pdf'), 'ContentType', 'vector');
close(fig);

nice = struct('y', 'Real GDP', 'p', 'GDP deflator', 'rd', 'Real R&D', 'rdnom', 'Nominal R&D', ...
    'rdshare', 'R&D / GDP', 'r', '1-year yield', 'ebp', 'EBP', 'baa', 'Baa spread', ...
    'epu', 'EPU', 'vix', 'VIX', 'eq', 'Equipment', 'sw', 'Software');
fig = figure('Units', 'inches', 'Position', [1 1 7 1.4*n], 'Color', 'w');
for i = 1:n
    for s = 1:k
        ax = subplot(n, k, (i-1)*k + s);
        q = pctl(R.irf(:, :, i, s), [5 16 50 84 95]);
        plot_band(ax, hh, q(1, :), q(2, :), q(3, :), q(4, :), q(5, :), ...
            sprintf('%s: %s', R.shockNames{s}, nice.(vars{i})));
    end
end
exportgraphics(fig, fullfile(cfg.figDir, 'fig_irf_all.pdf'), 'ContentType', 'vector');
close(fig);

%% R&D response table: posterior medians, 68% bands, P(response < 0), FEVD
fe = fevd_draws(R, M, H);
body = cell(k, 0);
hs = [0 4 8 20];
for s = 1:k
    row = {R.shockNames{s}};
    for h = hs
        q = pctl(R.irf(:, h+1, ird, s), [16 50 84]);
        row{end+1} = sprintf('%.2f [%.2f, %.2f]', q(2), q(1), q(3)); %#ok<AGROW>
    end
    row{end+1} = sprintf('%.2f', mean(R.irf(:, 9, ird, s) < 0));
    q = pctl(fe(:, 21, ird, s), [16 50 84]);
    row{end+1} = sprintf('%.2f [%.2f, %.2f]', q(2), q(1), q(3));
    body(s, 1:numel(row)) = row;
end
write_tex_table(fullfile(cfg.tabDir, 'tab_rd_responses.tex'), ...
    {'Shock', '$h=0$', '$h=4$', '$h=8$', '$h=20$', '$P(<0)_{h=8}$', 'FEVD$_{h=20}$'}, body, 'lcccccc');

%% Historical decomposition of R&D (median-target draw)
d = median_target(R.irf);
[contrib, baseline, e] = hist_decomp(M, R.B(:, :, d), R.A(:, :, d));
yr = floor(M.qidx/4) + mod(M.qidx, 4)/4;
fig = figure('Units', 'inches', 'Position', [1 1 7 2.4], 'Color', 'w');
ax = axes(fig); hold(ax, 'on');
plot(ax, yr, M.Y(:, ird) - baseline(:, ird), 'k-', 'LineWidth', 1.6);
cols = lines(k + 1);
for s = 1:k
    plot(ax, yr, contrib(:, ird, s), '-', 'Color', cols(s, :), 'LineWidth', 1.1);
end
plot(ax, yr, sum(contrib(:, ird, k+1:end), 3), ':', 'Color', cols(k+1, :), 'LineWidth', 1.1);
yline(ax, 0, 'k-');
legend(ax, [{'Total deviation from baseline'}, R.shockNames, {'Other shocks'}], ...
    'Location', 'southoutside', 'Orientation', 'horizontal', 'FontSize', 7, 'Box', 'off');
ylabel(ax, 'Percent', 'FontSize', 8);
set(ax, 'FontSize', 8, 'Box', 'on');
exportgraphics(fig, fullfile(cfg.figDir, 'fig_hd_rd.pdf'), 'ContentType', 'vector');
close(fig);

%% Persistence: lag-augmented local projections on the identified shocks
Zl = M.X(:, 1:M.ndet);
outcomes = {cfg.roleRD, 'y'};
labels = {'Real R&D investment', 'Real GDP'};
shocksLP = find(ismember(R.shockNames, {'Credit supply', 'Uncertainty'}));
[~, loc] = ismember(M.qidx, D.qidx);
% VAR-implied average response over the same window, per posterior draw, for
% the outcomes that are in the VAR (the LP is the model-free counterpart).
w = cfg.lrWindow;
nd = size(R.B, 3);
lrVAR = zeros(nd, n, k);
% The R&D capital stock is a fixed filter of past investment (Sliker, 2007,
% eq. 1), so its response is computed from the investment response of each
% draw rather than estimated separately:
%   dlnK_h = a dlnK_{h-1} + (1 - dq/2)(I/K) dlnI_h,  a = (1 - dq)/(1 + g_K).
[~, locS] = ismember(M.qidx, D.qidx);
dq = 1 - (1 - cfg.rdDepreciation)^(1/4);
iK = mean(exp(D.rd(locS)/100)/4./exp(D.rdk(locS)/100));
gK = mean(diff(D.rdk(locS)))/100;
aK = (1 - dq)/(1 + gK);
lrK = zeros(nd, k);
irfK = zeros(nd, w(2)+1, k);
for dd = 1:nd
    irfL = var_irf(R.B(:, :, dd), n, M.p, M.ndet, R.A(:, 1:k, dd), w(2));
    lrVAR(dd, :, :) = mean(irfL(w(1)+1:w(2)+1, :, :), 1);
    for s = 1:k
        kp = filter((1 - dq/2)*iK, [1 -aK], irfL(:, ird, s));
        irfK(dd, :, s) = kp;
        lrK(dd, s) = mean(kp(w(1)+1:w(2)+1));
    end
end
out.lrVAR = lrVAR;
out.irfK = irfK;
out.lrK = lrK;

fig = figure('Units', 'inches', 'Position', [1 1 7 4.2], 'Color', 'w');
body = cell(0, 5);
for a = 1:numel(shocksLP)
    s = shocksLP(a);
    sh = e(:, s)/std(e(:, s));
    for b = 1:numel(outcomes)
        yv = D.(outcomes{b})(loc);
        L = local_projection(yv, sh, M.Y, Zl, cfg.lpH, cfg.lpL, cfg.lrWindow);
        ax = subplot(numel(shocksLP), numel(outcomes), (a-1)*numel(outcomes) + b);
        plot_band(ax, L.h, L.beta - L.c90.*L.se, L.beta - L.c68.*L.se, L.beta, L.beta + L.c68.*L.se, ...
            L.beta + L.c90.*L.se, sprintf('%s: %s', R.shockNames{s}, labels{b}));
        iv = find(strcmp(vars, outcomes{b}));
        varCell = '--';
        if ~isempty(iv)
            q = pctl(lrVAR(:, iv, s), [16 50 84]);
            varCell = sprintf('%.2f [%.2f, %.2f]', q(2), q(1), q(3));
        end
        body(end+1, :) = {R.shockNames{s}, strrep(labels{b}, '&', '\&'), ...
            sprintf('%.2f (%.2f)', L.lr.beta, L.lr.se), sprintf('%.3f', L.lr.pval), varCell}; %#ok<AGROW>
        out.lp.(matlab.lang.makeValidName([R.shockNames{s} '_' outcomes{b}])) = L;
    end
    q = pctl(lrK(:, s), [16 50 84]);
    body(end+1, :) = {R.shockNames{s}, 'R\&D capital stock (implied)', '--', '--', ...
        sprintf('%.2f [%.2f, %.2f]', q(2), q(1), q(3))}; %#ok<AGROW>
end
exportgraphics(fig, fullfile(cfg.figDir, 'fig_lp.pdf'), 'ContentType', 'vector');
close(fig);
write_tex_table(fullfile(cfg.tabDir, 'tab_lp_longrun.tex'), ...
    {'Shock', 'Outcome', sprintf('LP avg. $h=%d$--$%d$ (s.e.)', w), '$p$-value', 'VAR avg. [68\%]'}, body, 'llccc');
end

function S = summarize_rd(R, ird)
S = struct();
for s = 1:numel(R.shockNames)
    f = matlab.lang.makeValidName(R.shockNames{s});
    for h = [4 8 20]
        q = pctl(R.irf(:, h+1, ird, s), [16 50 84]);
        S.(f).(sprintf('h%d', h)) = q;
    end
    S.(f).pneg8 = mean(R.irf(:, 9, ird, s) < 0);
end
S.accRate = R.accRate;
end
