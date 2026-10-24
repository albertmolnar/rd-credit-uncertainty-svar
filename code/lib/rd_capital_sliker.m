function S = rd_capital_sliker(I, t, delta, varargin)
% R&D capital stock by geometric perpetual inventory, following Sliker (2007),
% "2007 R&D Satellite Account Methodologies: R&D Capital Stocks and Net Rates of
% Return", BEA/NSF R&D Satellite Account background paper.
%
%   K_t = (1 - d_t) K_{t-1} + (1 - d_t/2) I_t                        eq. (1)
%   K_{t0-1} = I_{t0} (1 - d/2) / (g_K + d)                           eq. (3)
%   ln I_t = b + m t + e_t,   g_I = exp(m) - 1                         eq. (4)
% The regression (4) uses the early subsample that ends at the first peak of the
% R^2 sequence after its initial plunge. With time-varying depreciation,
% d_t = c + n t (eq. 9), and g_K solves the quadratic (7), positive root (8).
%
% I      real investment flow per period (column, contiguous)
% t      period index (e.g. calendar years)
% delta  depreciation rate per period: scalar, or vector (time-varying variant)
% Options: 'benchmark' 'r2' (default) | 'full' | 'zero';  'gI' to impose the
%          investment growth rate (e.g. from the annual estimate);  'minLen'
%          minimum subsample length for the R^2 search (default 8).

p = inputParser;
addParameter(p, 'benchmark', 'r2');
addParameter(p, 'gI', []);
addParameter(p, 'minLen', 8);
parse(p, varargin{:});
o = p.Results;

I = I(:); t = t(:);
N = numel(I);
if isscalar(delta), d = repmat(delta, N, 1); else, d = delta(:); end

% --- growth of investment, eq. (4) -----------------------------------------
y = log(I);
r2 = nan(N, 1); mm = nan(N, 1); se = nan(N, 1);
for k = 3:N
    X = [ones(k, 1), t(1:k)];
    b = X\y(1:k);
    e = y(1:k) - X*b;
    r2(k) = 1 - (e'*e)/sum((y(1:k) - mean(y(1:k))).^2);
    V = (e'*e)/(k - 2)*inv(X'*X);
    mm(k) = b(2);
    se(k) = sqrt(V(2, 2));
end
switch o.benchmark
    case 'full'
        kstar = N;
    otherwise
        kstar = first_post_plunge_peak(r2, o.minLen);
end
S.r2 = r2;
S.kstar = kstar;
S.subsample = [t(1) t(kstar)];
S.m = mm(kstar);
S.m_se = se(kstar);
S.R2 = r2(kstar);
if isempty(o.gI)
    gI = exp(S.m) - 1;
else
    gI = o.gI;
end
S.gI = gI;

% --- growth of capital: constant or time-varying depreciation ---------------
d1 = d(1);
X = [ones(N, 1), t];
nb = X\d;
n = nb(2);                                         % eq. (9) slope
if abs(n) < 1e-12
    gK = gI;
else
    a = (gI - d1) - n/(2 - d1);
    disc = a^2/4 + d1*gI - 2*n/(2 - d1);
    if disc < 0, error('Depreciation trend too steep for eq. (8).'); end
    gK = a/2 + sqrt(disc);                         % eq. (8)
end
S.n = n;
S.gK = gK;

% --- benchmark and accumulation, eqs. (3) and (1) ---------------------------
if strcmp(o.benchmark, 'zero')
    K0 = 0;
else
    K0 = I(1)*(1 - d1/2)/(gK + d1);
end
K = zeros(N, 1);
B = zeros(N, 1);                                   % surviving benchmark capital
Kprev = K0; Bprev = K0;
for i = 1:N
    K(i) = (1 - d(i))*Kprev + (1 - d(i)/2)*I(i);
    B(i) = (1 - d(i))*Bprev;
    Kprev = K(i); Bprev = B(i);
end
S.K0 = K0;
S.K = K;
S.benchShare = B./K;
S.D = d.*([K0; K(1:end-1)] + I/2);                 % depreciation flow, D_t = d (K_{t-1} + I_t/2)
S.delta = d;
end

function k = first_post_plunge_peak(r2, minLen)
% Sliker: the R^2 sequence first plunges from the two-observation default of 1;
% the subsample stops at the first peak after that plunge. A minimum length
% guards against spurious early peaks.
N = numel(r2);
k = N;
[~, trough] = min(r2(3:min(N, 3 + 2*minLen)));
trough = trough + 2;
for j = max(trough + 1, minLen):N-1
    if r2(j) > r2(j-1) && r2(j) >= r2(j+1)
        k = j;
        return
    end
end
end
