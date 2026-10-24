function D = build_dataset(cfg)
% Builds the quarterly analysis dataset from the raw files.
% Log levels are scaled by 100, so responses read as percent deviations.
% D.Properties.UserData.notes records any fallback used.

raw = struct();
qmin = inf; qmax = -inf;
for i = 1:size(cfg.series, 1)
    id = cfg.series{i, 1};
    f = fullfile(cfg.rawDir, [id '.csv']);
    if ~isfile(f)
        raw.(id) = [];
        continue
    end
    [q, v] = read_series_quarterly(f, cfg.series{i, 4});
    raw.(id) = [q v];
    qmin = min(qmin, min(q));
    qmax = max(qmax, max(q));
end
notes = strings(0, 1);

g = (qmin:qmax)';
X = @(id) on_grid(raw, id, g);
D = table(g, quarter_label(g), 'VariableNames', {'qidx', 'quarter'});

rgdp = X('GDPC1');
if all(isnan(rgdp))
    rgdp = X('GDP')./X('GDPDEF')*100;
    notes(end+1) = "Real GDP computed as GDP/GDPDEF*100.";
end
rd = X('RD_REAL');
if all(isnan(rd))
    if any(~isnan(X('RD_PRICE')))
        rd = X('RD_NOM')./X('RD_PRICE')*100;
        notes(end+1) = "Real R&D computed as nominal R&D / R&D price index.";
    else
        rd = X('RD_NOM')./X('GDPDEF')*100;
        notes(end+1) = "PROVISIONAL: real R&D = nominal R&D / GDP deflator (no R&D price index supplied).";
    end
end

D.y       = 100*log(rgdp);
D.p       = 100*log(X('GDPDEF'));
D.rd      = 100*log(rd);
D.rdnom   = 100*log(X('RD_NOM'));
D.rdshare = 100*X('RD_NOM')./X('GDP');
D.rdk     = 100*log(rd_capital_quarterly(rd, g, cfg.rdDepreciation));
D.eq      = 100*log(X('EQ_REAL'));
D.sw      = 100*log(X('SW_REAL'));
D.epu     = 100*log(X('EPU'));
D.vix     = 100*log(X('VIXCLS'));
D.ebp     = X('EBP');
D.baa     = X('BAA10Y');
D.r       = X('GS1');
D.Properties.UserData.notes = notes;
for k = 1:numel(notes), fprintf('NOTE: %s\n', notes(k)); end

% Plausibility check: private R&D investment is roughly 2-3% of GDP in the 2010s.
k = find(D.qidx == 2019*4 + 3);
if ~isempty(k) && ~isnan(D.rdshare(k)) && (D.rdshare(k) < 1.5 || D.rdshare(k) > 4)
    warning('R&D share of GDP in 2019Q4 is %.2f%% -- check the R&D series.', D.rdshare(k));
end
end

function x = on_grid(raw, id, g)
x = nan(size(g));
if ~isfield(raw, id) || isempty(raw.(id)), return; end
[tf, loc] = ismember(raw.(id)(:, 1), g);
x(loc(tf)) = raw.(id)(tf, 2);
end

function K = rd_capital_quarterly(Isaar, q, deltaA)
% Quarterly Sliker (2007) stock: flows Isaar/4, quarterly rate 1-(1-deltaA)^(1/4),
% mid-quarter convention, benchmark growth taken from the annual eq. (4) estimate.
K = nan(size(Isaar));
ok = find(~isnan(Isaar));
if numel(ok) < 40, return; end
ok = ok(1):ok(end);
Iq = Isaar(ok)/4;
yr = floor(q(ok)/4);
full = accumarray(yr - yr(1) + 1, 1) == 4;
Ia = accumarray(yr - yr(1) + 1, Iq);
ya = (yr(1):yr(end))';
A = rd_capital_sliker(Ia(full), ya(full), deltaA);
dq = 1 - (1 - deltaA)^(1/4);
gq = (1 + A.gK)^(1/4) - 1;
S = rd_capital_sliker(Iq, (1:numel(Iq))', dq, 'gI', gq);
K(ok) = S.K;
end
