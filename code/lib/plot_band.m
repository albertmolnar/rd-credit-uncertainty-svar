function plot_band(ax, h, lo90, lo68, mid, hi68, hi90, ttl)
% Median (or point estimate) with 68% and 90% bands.
h = h(:)';
hold(ax, 'on');
fill(ax, [h fliplr(h)], [lo90(:)' fliplr(hi90(:)')], [0.82 0.88 0.96], 'EdgeColor', 'none');
fill(ax, [h fliplr(h)], [lo68(:)' fliplr(hi68(:)')], [0.60 0.73 0.90], 'EdgeColor', 'none');
plot(ax, h, mid(:)', 'LineWidth', 1.5, 'Color', [0.05 0.20 0.50]);
yline(ax, 0, 'k-');
xlim(ax, [h(1) h(end)]);
title(ax, ttl, 'FontWeight', 'normal', 'FontSize', 9);
xlabel(ax, 'Quarters after shock', 'FontSize', 8);
set(ax, 'FontSize', 8, 'Box', 'on');
end
