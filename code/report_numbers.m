%% Numbers quoted in the text of the paper (run after main.m)
clear; clc;
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, 'lib'));
load(fullfile(here, 'output', 'baseline.mat'), 'base');
M = base.M; R = base.R;
fprintf('Posterior: %d accepted draws from %d reduced-form draws; explosive share %.1f%%\n', ...
    R.nKept, R.nTried, 100*R.explosiveShare);
fprintf('Estimation sample %s, T = %d\n', base.sample, M.T);

% Peak responses of R&D
ird = find(strcmp(M.vars, 'rd'));
for s = 1:numel(R.shockNames)
    med = median(R.irf(:, :, ird, s), 1);
    [pk, h] = min(med);
    fprintf('%-16s R&D median trough %.2f%% at h = %d; P(<0) at trough %.2f\n', ...
        R.shockNames{s}, pk, h - 1, mean(R.irf(:, h, ird, s) < 0));
end
% Impact responses used to describe the identified shocks (1 s.d.)
for s = 1:numel(R.shockNames)
    q = median(R.irf(:, 1, :, s), 1);
    fprintf('%-16s impact: y %.2f p %.2f rd %.2f r %.2f ebp %.2f epu %.2f\n', R.shockNames{s}, q);
end

% Historical decomposition: average contribution to R&D over episodes
d = median_target(R.irf);
[contrib, baseline] = hist_decomp(M, R.B(:, :, d), R.A(:, :, d));
dev = M.Y(:, ird) - baseline(:, ird);
ep = {'2001Q1', '2003Q4'; '2008Q3', '2010Q4'; '2011Q1', '2014Q4'};
lab = string(quarter_label(M.qidx));
for e = 1:size(ep, 1)
    in = find(lab == ep{e, 1}):find(lab == ep{e, 2});
    c = squeeze(mean(contrib(in, ird, :), 1));
    fprintf('%s-%s: total dev %.1f; credit %.1f; uncertainty %.1f; supply %.1f; other %.1f\n', ...
        ep{e, 1}, ep{e, 2}, mean(dev(in)), c(1), c(2), c(3), sum(c(4:end)));
end

% Recent R&D intensity (nominal R&D / GDP), from the dataset
D = readtable(fullfile(here, 'output', 'dataset_quarterly.csv'), 'TextType', 'string');
for qq = ["1986Q1", "2000Q4", "2019Q4", "2022Q4", "2025Q4"]
    k = find(D.quarter == qq);
    if ~isempty(k), fprintf('R&D/GDP %s: %.2f%%\n', qq, D.rdshare(k)); end
end
k1 = find(D.quarter == "2019Q4"); k2 = find(D.quarter == "2025Q4");
if ~isempty(k2)
    fprintf('Nominal R&D growth 2019Q4-2025Q4: %.1f%% p.a.\n', (exp((D.rdnom(k2) - D.rdnom(k1))/100/6) - 1)*100);
end
