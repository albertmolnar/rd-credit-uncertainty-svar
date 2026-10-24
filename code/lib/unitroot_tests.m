function UR = unitroot_tests(D, names, trendVars, q0, q1)
% ADF (lag length by BIC, max 8) and KPSS tests in levels and first differences.
% Levels of trending variables include a linear trend; differences use a constant.
% ADF 5% critical values from MacKinnon (2010); KPSS 5% critical values from
% Kwiatkowski et al. (1992).

rows = cell(numel(names), 9);
for i = 1:numel(names)
    in = D.qidx >= q0 & D.qidx <= q1;
    y = D.(names{i})(in);
    y = y(~isnan(y));
    tr = any(strcmp(trendVars, names{i}));
    [adfL, cvL] = adf(y, tr, 8);
    [adfD, cvD] = adf(diff(y), false, 8);
    [kpL, kcvL] = kpss(y, tr);
    [kpD, kcvD] = kpss(diff(y), false);
    det = "c";
    if tr, det = "c,t"; end
    rows(i, :) = {string(names{i}), det, adfL, cvL, adfD, cvD, kpL, kcvL, kpD};
    UR.kpssDiffCV = kcvD;
end
UR.table = cell2table(rows, 'VariableNames', ...
    {'variable', 'det', 'adf_level', 'adf_level_cv5', 'adf_diff', 'adf_diff_cv5', ...
     'kpss_level', 'kpss_level_cv5', 'kpss_diff'});
end

function [tstat, cv] = adf(y, trend, pmax)
T = numel(y);
dy = diff(y);
first = pmax + 1;                 % common sample across lag lengths
best = inf; tstat = nan;
for k = 0:pmax
    idx = (first:T-1)';
    X = [ones(numel(idx), 1), y(idx)];
    if trend, X = [X, idx]; end %#ok<AGROW>
    for l = 1:k
        X = [X, dy(idx - l)]; %#ok<AGROW>
    end
    yy = dy(idx);
    b = X\yy;
    e = yy - X*b;
    nobs = numel(yy);
    bic = log(e'*e/nobs) + size(X, 2)*log(nobs)/nobs;
    if bic < best
        best = bic;
        s2 = (e'*e)/(nobs - size(X, 2));
        V = s2*inv(X'*X);
        tstat = b(2)/sqrt(V(2, 2));
    end
end
n = T - first;
if trend
    cv = -3.4126 - 4.039/n - 17.83/n^2;
else
    cv = -2.8621 - 2.738/n - 8.36/n^2;
end
end

function [stat, cv] = kpss(y, trend)
T = numel(y);
X = ones(T, 1);
if trend, X = [X, (1:T)']; end
e = y - X*(X\y);
S = cumsum(e);
l = floor(4*(T/100)^(1/4));
lrv = (e'*e)/T;
for j = 1:l
    lrv = lrv + 2*(1 - j/(l + 1))*(e(j+1:end)'*e(1:end-j))/T;
end
stat = sum(S.^2)/(T^2*lrv);
if trend, cv = 0.146; else, cv = 0.463; end
end
