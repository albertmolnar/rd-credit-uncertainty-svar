function [irf, Psi] = var_irf(B, n, p, ndet, Aimp, H)
% Reduced-form MA coefficients Psi(:,:,h+1) and structural responses
% irf(h+1, variable, shock) = Psi_h * Aimp for h = 0..H.

Al = zeros(n, n, p);
for l = 1:p
    Al(:, :, l) = B(ndet + (l-1)*n + (1:n), :)';
end
Psi = zeros(n, n, H+1);
Psi(:, :, 1) = eye(n);
for h = 1:H
    acc = zeros(n);
    for l = 1:min(h, p)
        acc = acc + Al(:, :, l)*Psi(:, :, h+1-l);
    end
    Psi(:, :, h+1) = acc;
end
if isempty(Aimp)
    irf = [];
    return
end
k = size(Aimp, 2);
irf = zeros(H+1, n, k);
for h = 0:H
    irf(h+1, :, :) = reshape(Psi(:, :, h+1)*Aimp, [1 n k]);
end
end
