function LS = var_lagselect(Y, pmax, Z)
% AIC, BIC (Schwarz) and Hannan-Quinn for p = 1..pmax on a common sample.

T = size(Y, 1);
n = size(Y, 2);
Te = T - pmax;
crit = zeros(pmax, 3);
for p = 1:pmax
    rows = (pmax - p + 1):T;
    M = var_estimate(Y(rows, :), p, Z(rows, :));
    ld = log(det(M.SigmaML));
    k = p*n^2;
    crit(p, :) = [ld + 2*k/Te, ld + log(Te)*k/Te, ld + 2*log(log(Te))*k/Te];
end
[~, best] = min(crit);
LS.crit = crit;
LS.best = best;          % [AIC BIC HQ]
LS.names = {'AIC', 'BIC', 'HQ'};
end
