function out = local_projection(y, shock, W, Z, H, L, window)
% Lag-augmented local projections (Jorda, 2005; Montiel Olea and Plagborg-Moller, 2021):
%   y_{t+h} - y_{t-1} = a_h + b_h shock_t + sum_{l=1}^{L+1} G_l W_{t-l} + e_{t+h}.
% Lags of the shock are not added: a VAR shock is a linear combination of current
% and lagged W, so its lags are already spanned by the lags of W.
%
% Inference uses the equal-weighted cosine (EWC) HAC estimator with fixed-b
% Student-t critical values (Lazarus, Lewis, Stock and Watson, 2018). In
% Monte Carlo checks at this sample size (T ~ 136, ~30 controls), EHW and
% Newey-West intervals under-cover badly at long horizons; EWC does not.
%
% Also estimates the average response over h = window(1)..window(2) in one
% regression -- the test of a lasting level effect.
% y: T x 1 (100*log level), shock: T x 1, W: T x q controls, Z: T x d deterministic terms.

T = numel(y);
C = [Z, lagmat(W, L+1)];

out.h = (0:H)';
out.beta = nan(H+1, 1);
out.se = nan(H+1, 1);
out.df = nan(H+1, 1);
out.nobs = nan(H+1, 1);
for h = 0:H
    dep = nan(T, 1);
    dep(2:T-h) = y(2+h:T) - y(1:T-h-1);
    [b, se, nu, nobs] = ewc_slope(dep, shock, C);
    out.beta(h+1) = b;
    out.se(h+1) = se;
    out.df(h+1) = nu;
    out.nobs(h+1) = nobs;
end
% Band multipliers that match the t(nu) reference distribution.
out.c68 = arrayfun(@(v) tquantile(0.84, v), out.df);
out.c90 = arrayfun(@(v) tquantile(0.95, v), out.df);

h1 = window(1); h2 = window(2);
dep = nan(T, 1);
for t = 2:T-h2
    dep(t) = mean(y(t+h1:t+h2)) - y(t-1);
end
[b, se, nu] = ewc_slope(dep, shock, C);
out.lr.beta = b;
out.lr.se = se;
out.lr.df = nu;
out.lr.t = b/se;
out.lr.pval = 2*tcdf_upper(abs(out.lr.t), nu);
out.lr.window = window;
end

function [b, se, nu, Tn] = ewc_slope(y, x, C)
% Slope on x after partialling out C (Frisch-Waugh), with EWC standard error.
ok = all(~isnan([y x C]), 2);
y = y(ok); x = x(ok); C = C(ok, :);
xr = x - C*(C\x);
yr = y - C*(C\y);
b = (xr'*yr)/(xr'*xr);
u = yr - xr*b;
Tn = numel(y);
z = xr.*u;
nu = max(2, floor(0.4*Tn^(2/3)));
t = (1:Tn)';
Om = 0;
for j = 1:nu
    lam = sqrt(2/Tn)*sum(z.*cos(pi*j*(t - 0.5)/Tn));
    Om = Om + lam^2/nu;
end
q = (xr'*xr)/Tn;
k = size(C, 2) + 1;
se = sqrt(Om/q^2/Tn*Tn/(Tn - k));     % small-sample degrees-of-freedom correction
end

function p = tcdf_upper(t, nu)
% P(T > t) for Student t with nu degrees of freedom, t >= 0.
p = 0.5*betainc(nu/(nu + t^2), nu/2, 0.5);
end

function q = tquantile(p, nu)
% Quantile of Student t (p > 0.5) by bisection.
lo = 0; hi = 50;
for it = 1:60
    mid = (lo + hi)/2;
    if 1 - tcdf_upper(mid, nu) < p, lo = mid; else, hi = mid; end
end
q = (lo + hi)/2;
end

function Lm = lagmat(X, L)
[T, q] = size(X);
Lm = nan(T, q*L);
for l = 1:L
    Lm(l+1:end, (l-1)*q + (1:q)) = X(1:end-l, :);
end
end
