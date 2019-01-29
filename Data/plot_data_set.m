function plot_data_set(num)
% Plot synthetic data set
%
% -------------------------------------------------------------------------
close all

% Generate data

dataopt = num;
[mu, sig, clustersize, noise] = define_data(dataopt);
[data, labels] = generate_data_noise(mu, sig, clustersize,...
    noise(1), noise(2), true);

figure, hold on, grid on, axis equal
for i = 2:size(data,1)
    plot([data(i-1, 1), data(i, 1)],[data(i-1, 2), data(i, 2)],...
        'k.-','markersize', 20,'color',[1 1 1]*(1-i/size(data,1)))
end
axis([-8 9 -4 8])
set(gca,'FontName','Candara', 'FontSize',14,'Xtick',-8:2:9)

nlabel = max(labels);
for i = 1:nlabel
    fi = kf_choose_the_most_similar(data(labels == i,:));
    lfi = find(labels == i);
    
    plot(data(lfi(fi),1),data(lfi(fi),2),'r+','markersize',18,...
        'linewidth',1.8)
    plot(data(lfi(fi),1),data(lfi(fi),2),'ro','markersize',10,...
        'linewidth',1.8)
end

xlabel('{\it x}','FontName','Times')
ylabel('{\it y}','FontName','Times')
end