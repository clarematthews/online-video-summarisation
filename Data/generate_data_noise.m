function [data, labels] = generate_data_noise(mu, sig, csize,...
    noisedim, noise, seed)
% mu = cluster centres (num cluster x data dim)
% sig = covariance matrix of each cluster
% csize = number of points in each cluster
% noisedim = number of additional noise dimensions
% noise = level of noise (std)
% seed = use seed for random numbers

if seed
    rng(1959)
end

[ncluster, ndim] = size(mu);
npoints = sum(csize);
data = zeros(npoints, ndim + noisedim);
labels = zeros(1, npoints);
for i = 1:ncluster
    ci = mvnrnd(mu(i, :), sig(:, :, i), csize(i));
    crange = sum(csize(1:i - 1)) + 1:sum(csize(1:i));
    data(crange, 1:ndim) = ci;
    if noisedim > 0
        ni = mvnrnd(zeros(1, noisedim), noise*eye(noisedim), csize(i));
        data(crange, ndim + 1:ndim + noisedim) = ni;
    end
    
    labels(crange) = i;
end
end