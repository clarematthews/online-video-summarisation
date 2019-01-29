% Plot all 2-dimensional (plus 2 noise dimensions) data sets
%
% -------------------------------------------------------------------------
clear, clc, close all

figi = 0;
for dataopt = [1 2 5 6 7]
    [mu, sig, clustersize, noise] = define_data(dataopt);
    [data, labels] = generate_data_noise(mu, sig, clustersize,...
        noise(1), noise(2), true);
    ncluster = size(mu, 1);
    mcol = rand(ncluster, 3);
    
    figure, hold on, grid on, axis equal

    N = size(data,1);
    for i = 2:N
        col = i/N;
        plot([data(i-1, 1), data(i, 1)],[data(i-1, 2), data(i, 2)],...
            'k.-','markersize', 25, 'color', [1 1 1]*(1 - col))
    end
    
    for j = 1:numel(clustersize)
        fj = kf_choose_the_most_similar(data(labels == j,:));
        lfj = find(labels == j);
        plot(data(lfj(fj), 1), data(lfj(fj), 2), 'r+', 'markersize', 15,...
            'linewidth', 1.8)
        plot(data(lfj(fj), 1), data(lfj(fj), 2), 'ro', 'markersize', 8,...
            'linewidth', 1.8)
    end
    set(gca,'YTickLabel',[], 'XTickLabel',[]);
end
