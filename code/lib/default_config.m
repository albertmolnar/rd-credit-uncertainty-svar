function cfg = default_config(rootDir)
% Baseline settings. Robustness specifications in main.m override single fields.

cfg.seed   = 20261015;
cfg.rawDir = fullfile(rootDir, 'data_raw');      % downloaded source files (= data vintage)
cfg.outDir = fullfile(rootDir, 'output');
cfg.figDir = fullfile(cfg.outDir, 'figures');
cfg.tabDir = fullfile(cfg.outDir, 'tables');

% ---------------------------------------------------------------------------
% Source series: {local id, URL, description, value column}.
% FRED IDs marked [verify] could not be checked when this file was written:
% open https://fred.stlouisfed.org/series/<ID> and confirm the title before use.
% An empty URL means "not configured": the series is skipped.
% ---------------------------------------------------------------------------
% Optional: manually downloaded FRED files (CSV/XLSX as saved from the FRED
% website, file name = FRED code, e.g. BAA10Y.xlsx). If this folder exists,
% import_local_data copies them into data_raw/ before any download is tried.
cfg.localDir = fullfile(rootDir, 'data_local');

fred = @(id) ['https://fred.stlouisfed.org/graph/fredgraph.csv?id=' id];
cfg.series = {
    'GDP',      fred('GDP'),             'Nominal GDP, bil. USD, SAAR (BEA via FRED: GDP)', 2
    'GDPC1',    '',                      'Real GDP, chained dollars -- not used: the paper uses GDP/GDPDEF', 2
    'GDPDEF',   fred('GDPDEF'),          'GDP implicit price deflator (BEA via FRED: GDPDEF)', 2
    'RD_NOM',   fred('Y006RC1Q027SBEA'), 'Private fixed investment: nonresidential IPP: research and development, nominal, SAAR (BEA via FRED: Y006RC1Q027SBEA)', 2
    'RD_REAL',  '',                      'Real R&D investment -- not used: the paper deflates nominal R&D by GDPDEF', 2
    'RD_PRICE', '',                      'R&D price index -- not used (not available to the author)', 2
    'EQ_REAL',  '',                      'Real private fixed investment: nonresidential equipment -- not supplied', 2
    'SW_REAL',  '',                      'Real private fixed investment: IPP: software -- not supplied', 2
    'EPU',      fred('USEPUINDXD'),      'Economic Policy Uncertainty index, daily news-based (Baker, Bloom and Davis via FRED: USEPUINDXD)', 2
    'EBP',      'https://www.federalreserve.gov/econres/notes/feds-notes/ebp_csv.csv', ...
                                         'Excess bond premium, monthly (Gilchrist and Zakrajsek; Federal Reserve Board)', 'ebp'
    'BAA10Y',   fred('BAA10Y'),          'Moody''s Baa yield minus 10-year Treasury, daily (via FRED: BAA10Y)', 2
    'GS1',      fred('DGS1'),            '1-year Treasury constant maturity yield, daily (via FRED: DGS1)', 2
    'VIXCLS',   fred('VIXCLS'),          'CBOE volatility index VIX, daily close (via FRED: VIXCLS)', 2
    };

% R&D capital (Sliker, 2007): annual geometric depreciation rate. BEA applies
% industry-specific rates (roughly 10-40%); 15% is a common aggregate value.
cfg.rdDepreciation = 0.15;

% ---------------------------------------------------------------------------
% Baseline VAR
% ---------------------------------------------------------------------------
cfg.vars        = {'y', 'p', 'rd', 'r', 'ebp', 'epu'};
cfg.roleRD      = 'rd';      % R&D measure (never restricted)
cfg.roleCredit  = 'ebp';     % credit-conditions indicator
cfg.roleUnc     = 'epu';     % uncertainty indicator
cfg.trendVars   = {'y', 'p', 'rd', 'rdk', 'eq', 'sw'};   % deterministic trend in unit-root tests
cfg.sampleStart = [1985 1];  % first quarter used (the first p quarters are pre-sample)
cfg.sampleEnd   = [2019 4];  % [] = last quarter available for all variables
cfg.covidDummies = false;    % impulse dummies 2020Q1-2021Q4 when the sample includes them
cfg.pmax        = 8;
cfg.p           = 4;

% ---------------------------------------------------------------------------
% Identification
% ---------------------------------------------------------------------------
cfg.ident           = 'sign';   % 'sign' or 'chol'
cfg.Hrestr          = 1;        % sign restrictions hold for h = 0..Hrestr quarters
cfg.creditRateRestr = true;     % credit shock lowers the 1-year yield
cfg.uncUseMacro     = true;     % uncertainty shock lowers output and prices
cfg.useNarrative    = false;
cfg.narrative = {               % {shock, [year quarter], sign of the structural shock}
    'Credit supply', [2008 4], +1    % post-Lehman collapse of credit intermediation
    'Uncertainty',   [2001 3], +1    % September 11 attacks
    'Uncertainty',   [2011 3], +1    % debt-ceiling standoff and US downgrade
    };

cfg.H            = 20;      % IRF horizon (quarters)
cfg.nDraws       = 1000;    % accepted posterior draws
cfg.nRotMax      = 2000;    % rotations tried per reduced-form draw
cfg.maxPostDraws = 50000;   % cap on reduced-form draws

% ---------------------------------------------------------------------------
% Local projections (persistence)
% ---------------------------------------------------------------------------
cfg.lpH      = 28;          % horizons 0..28 quarters
cfg.lpL      = 2;           % control lags (one extra lag is added: lag augmentation);
                            % 2 gave the best coverage in lp_coverage_check.m
cfg.lrWindow = [20 28];     % medium-run window for the average-response test
end
