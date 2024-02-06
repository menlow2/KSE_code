function calm_convergence_plotter(nametag, T, pts)

if ~(exist('T','var')) % Final time
    T = 10;
end
if ~(exist('pts','var')) % Number of epsilon points
    pts = 6;
end


epsilons = linspace(-14, -4, pts);
epsilons = 10.^(epsilons);

base_name = ['Convergence_T_' num2str(T)];
base_name = strrep(base_name,'.','p');

save_dir = makeFolder(mfilename, base_name, 'convergence_plots');

fh1 = figure;
set(fh1,'units','normalized','position',[0,0,.8,.8]);
Markers = ['o','+','*','s','d','v','^','<','>','p','x'];

load([nametag  '.mat'], 'LIL2_error', 'L2H2_error', 'LILI_error');

for i = 1:6
tag = ['ConvergenceTest' num2str(i)];

MIN = min([LIL2_error(i,:) L2H2_error(i,:) LILI_error(i,:)]);
MAX = max([LIL2_error(i,:) L2H2_error(i,:) LILI_error(i,:)]);
xlim( [epsilons(1) epsilons(end) ] );
ylim([eps MAX]);
if (i == 1) || (i == 4)
loglog(epsilons, 10^(MAX)*epsilons(), 'k-.', 'linewidth', 1); hold on;
else
loglog(epsilons(), 10^(MAX+8)*epsilons().^2, 'k-.', 'linewidth', 1); hold on;
end
loglog(epsilons, LIL2_error(i,:), 'MarkerIndices',floor(0.95*length(epsilons)),'marker',Markers(1),'MarkerSize',10, 'linewidth', 2, 'color', 'r'); hold on;
loglog(epsilons, L2H2_error(i,:), 'MarkerIndices',floor(0.85*length(epsilons)),'marker',Markers(2),'MarkerSize',10, 'linewidth', 2, 'color', 'g'); hold on;
loglog(epsilons, LILI_error(i,:), 'MarkerIndices',floor(0.75*length(epsilons)),'marker',Markers(3),'MarkerSize',10, 'linewidth', 2, 'color', 'b'); hold on;

if (i == 1) || (i == 4)
legend( '$\epsilon$',  ...
    '$L^\infty(L^2) $', ...
    '$L^2(H^2) $', ...
    '$L^\infty(L^\infty) $', ...
    'FontSize',16,'Location','northwest','Interpreter','Latex');
else
legend( '$\epsilon^2$',  ...
    '$L^\infty(L^2) $', ...
    '$L^2(H^2) $', ...
    '$L^\infty(L^\infty) $', ...
    'FontSize',16,'Location','northwest','Interpreter','Latex');
end

xlabel('epsilon');
ylabel('error');
title(['Convergence rates on time interval $[0,' num2str(T) ']$'], 'interpreter', 'latex')
axis('tight');
xlim( [epsilons(1) epsilons(end) ] );
ylim([eps MAX]);
savefig(fh1, [save_dir tag])
print(fh1,'-djpeg', [ save_dir tag '.png'] );
clf(fh1);

end

end