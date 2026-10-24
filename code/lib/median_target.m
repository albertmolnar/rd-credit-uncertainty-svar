function i = median_target(irfs)
% Index of the single draw whose IRFs are closest to the pointwise posterior
% median (Fry and Pagan, 2011), used for historical decompositions and LP shocks.
med = median(irfs, 1);
sdv = std(irfs, 0, 1);
sdv(sdv == 0) = 1;
z = (irfs - med)./sdv;
[~, i] = min(sum(reshape(z.^2, size(irfs, 1), []), 2));
end
