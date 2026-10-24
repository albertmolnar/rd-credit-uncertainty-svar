%% Credit supply, uncertainty and U.S. private R&D -- replication script (v2)
% Base MATLAB only (no toolboxes). Run from this folder:  >> main
%   1. download public source data (FRED, Federal Reserve Board) into data_raw/
%   2. build the quarterly dataset and descriptive / unit-root tables
%   3. baseline Bayesian VAR, agnostic sign restrictions (R&D never restricted):
%      lag selection, diagnostics, IRFs, FEVD, historical decomposition
%   4. persistence: lag-augmented local projections on the identified shocks
%   5. robustness specifications
% Outputs go to output/ (tables as .tex, figures as .pdf, log.txt).

clear; close all; clc;
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, 'lib'));
cfg = default_config(here);
rng(cfg.seed);
for d = {cfg.rawDir, cfg.outDir, cfg.figDir, cfg.tabDir}
    if ~exist(d{1}, 'dir'), mkdir(d{1}); end
end
diary(fullfile(cfg.outDir, 'log.txt'));
fprintf('Run started %s\n', datestr(now));

%% 1-2. Data
if isfolder(cfg.localDir), import_local_data(cfg); end
src = fetch_data(cfg);
writetable(src, fullfile(cfg.outDir, 'data_sources.csv'));
D = build_dataset(cfg);
writetable(D, fullfile(cfg.outDir, 'dataset_quarterly.csv'));

q0 = cfg.sampleStart(1)*4 + cfg.sampleStart(2) - 1;
q1 = cfg.sampleEnd(1)*4 + cfg.sampleEnd(2) - 1;
cand = {'y', 'p', 'rd', 'rdshare', 'rdk', 'r', 'ebp', 'baa', 'epu', 'vix', 'eq', 'sw'};
avail = cand(cellfun(@(v) any(~isnan(D.(v)(D.qidx >= q0 & D.qidx <= q1))), cand));

% Descriptive statistics over the baseline sample (growth rates annualized).
in = D.qidx >= q0 & D.qidx <= q1;
desc = {
    'Real R\&D investment growth (\%, ann.)', 4*[NaN; diff(D.rd(in))]
    'Real GDP growth (\%, ann.)',             4*[NaN; diff(D.y(in))]
    'GDP deflator inflation (\%, ann.)',      4*[NaN; diff(D.p(in))]
    'R\&D investment / GDP (\%)',            D.rdshare(in)
    '1-year Treasury yield (\%)',            D.r(in)
    'Excess bond premium (pp)',              D.ebp(in)
    'Baa--10y Treasury spread (pp)',         D.baa(in)
    'EPU index',                             exp(D.epu(in)/100)
    };
body = cell(size(desc, 1), 5);
for i = 1:size(desc, 1)
    x = desc{i, 2};
    x = x(~isnan(x));
    body(i, :) = {desc{i, 1}, sprintf('%.2f', mean(x)), sprintf('%.2f', std(x)), ...
        sprintf('%.2f', min(x)), sprintf('%.2f', max(x))};
end
write_tex_table(fullfile(cfg.tabDir, 'tab_descriptive.tex'), {'Variable', 'Mean', 'Std. dev.', 'Min', 'Max'}, body, 'lcccc');

UR = unitroot_tests(D, avail, cfg.trendVars, q0, q1);
writetable(UR.table, fullfile(cfg.outDir, 'unitroot.csv'));
body = cell(height(UR.table), 6);
for i = 1:height(UR.table)
    t = UR.table(i, :);
    body(i, :) = {char(t.variable), char(t.det), ...
        sprintf('%.2f (%.2f)', t.adf_level, t.adf_level_cv5), sprintf('%.2f (%.2f)', t.adf_diff, t.adf_diff_cv5), ...
        sprintf('%.3f (%.3f)', t.kpss_level, t.kpss_level_cv5), sprintf('%.3f (%.3f)', t.kpss_diff, UR.kpssDiffCV)};
end
write_tex_table(fullfile(cfg.tabDir, 'tab_unitroot.tex'), ...
    {'Variable', 'Det.', 'ADF level', 'ADF diff.', 'KPSS level', 'KPSS diff.'}, body, 'llcccc');

%% 3-4. Baseline
fprintf('\n=== Baseline ===\n');
base = run_spec(D, cfg, true);
save(fullfile(cfg.outDir, 'baseline.mat'), 'base', '-v7.3');

%% 5. Robustness
rob = cfg;
rob.nDraws = 500;
swvars = {'y', 'p', 'rd', 'sw', 'eq', 'r', 'ebp', 'epu'};
swvars = swvars(cellfun(@(v) any(~isnan(D.(v))), swvars));
specs = {
    'Baseline',                                  rob
    'Baa--Treasury spread instead of EBP',       setf(rob, 'vars', {'y', 'p', 'rd', 'r', 'baa', 'epu'}, 'roleCredit', 'baa')
    'VIX instead of EPU (from 1990)',            setf(rob, 'vars', {'y', 'p', 'rd', 'r', 'ebp', 'vix'}, 'roleUnc', 'vix', 'sampleStart', [1990 1])
    'R\&D share of GDP instead of log level',    setf(rob, 'vars', {'y', 'p', 'rdshare', 'r', 'ebp', 'epu'}, 'roleRD', 'rdshare')
    'Sample from 1994Q1',                        setf(rob, 'sampleStart', [1994 1])
    'Sample to latest, COVID dummies',           setf(rob, 'sampleEnd', [], 'covidDummies', true)
    'Narrative restrictions added',              setf(rob, 'useNarrative', true)
    'Restrictions on impact only',               setf(rob, 'Hrestr', 0)
    'Restrictions for $h=0$--$3$',               setf(rob, 'Hrestr', 3)
    'Two lags',                                  setf(rob, 'p', 2)
    'Six lags',                                  setf(rob, 'p', 6)
    'Credit shock: yield unrestricted',          setf(rob, 'creditRateRestr', false)
    'Uncertainty shock: minimal restrictions',   setf(rob, 'uncUseMacro', false)
    'Recursive, financial block last',           setf(rob, 'ident', 'chol')
    'Recursive, financial block first',          setf(rob, 'ident', 'chol', 'vars', {'epu', 'ebp', 'y', 'p', 'rd', 'r'})
    'Software and equipment added',              setf(rob, 'vars', swvars)
    };
hasExt = numel(swvars) > 6;
if ~hasExt, specs = specs(1:end-1, :); end      % software/equipment series not supplied
body = cell(size(specs, 1), 7);
results = cell(size(specs, 1), 1);
for i = 1:size(specs, 1)
    fprintf('\n=== Robustness %d: %s ===\n', i, specs{i, 1});
    rng(cfg.seed + 100 + i);       % own seed per specification: results do not depend on the others
    try
        o = run_spec(D, specs{i, 2}, false);
        results{i} = o;
        S = o.summary;
        c = @(f, h) sprintf('%.2f [%.2f, %.2f]', S.(f).(h)(2), S.(f).(h)(1), S.(f).(h)(3));
        body(i, :) = {specs{i, 1}, o.sample, c('CreditSupply', 'h8'), sprintf('%.2f', S.CreditSupply.pneg8), ...
            c('Uncertainty', 'h8'), sprintf('%.2f', S.Uncertainty.pneg8), c('Uncertainty', 'h20')};
    catch ME
        warning('Specification "%s" failed: %s', specs{i, 1}, ME.message);
        body(i, :) = {specs{i, 1}, 'failed', '', '', '', '', ''};
    end
end
write_tex_table(fullfile(cfg.tabDir, 'tab_robustness.tex'), ...
    {'Specification', 'Sample', 'Credit, $h=8$', '$P(<0)$', 'Uncert., $h=8$', '$P(<0)$', 'Uncert., $h=20$'}, body, 'llccccc');

% Extension: R&D vs software vs equipment responses.
ext = results{end};
if hasExt && ~isempty(ext)
    vv = ext.cfg.vars;
    show = intersect({'rd', 'sw', 'eq'}, vv, 'stable');
    fig = figure('Units', 'inches', 'Position', [1 1 7 2.2*2], 'Color', 'w');
    for a = 1:2
        for b = 1:numel(show)
            ax = subplot(2, numel(show), (a-1)*numel(show) + b);
            q = pctl(ext.R.irf(:, :, strcmp(vv, show{b}), a), [5 16 50 84 95]);
            plot_band(ax, 0:cfg.H, q(1, :), q(2, :), q(3, :), q(4, :), q(5, :), ...
                sprintf('%s: %s', ext.R.shockNames{a}, show{b}));
        end
    end
    exportgraphics(fig, fullfile(cfg.figDir, 'fig_irf_intangibles.pdf'), 'ContentType', 'vector');
    close(fig);
end

save(fullfile(cfg.outDir, 'robustness.mat'), 'results', 'specs', '-v7.3');

%% 6. Size-adjusted inference for the local projections (clears the workspace)
lp_coverage_check;

%% 7. Copy tables/figures into v2/paper and build the Overleaf zip
export_to_paper;
fprintf('\nRun finished %s\n', datestr(now));
diary off;
