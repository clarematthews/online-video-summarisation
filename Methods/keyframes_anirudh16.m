function S = keyframes_anirudh16(X, K, C, B, randupdate)
% K = number of keyframes
% C = normalising factor Ci for frame number i
% B = balance between error and diversity penalties
% randupdate = % frames to use randomly as update to avoid local minima

% -------------------------------------------------------------------------

[nframes, ndim] = size(X);

% Initialisation

S = 1:K; % keyframes
mu = X(1:K, :);

% Keyframes must be unique to avoid singular matrix
mu = unique(mu, 'rows');
n = K + 1;
while size(mu, 1) < K 
    mu = [mu; X(n, :)]; %#ok<AGROW>
    mu = unique(mu, 'rows');
    n = n + 1;
end

% Convex hull measure requires the number of clusters to be greater than
% the dimensions of the feature space. Use PCA to reduce dimensionality if
% necessary

if K <= ndim
    mu = zscore(mu);
    coeff = pca(mu);
    X = X*coeff;
    X = X(:, 1:K - 1);
    mu = mu*coeff;
    mu = mu(:, 1:K - 1);
    mu = double(mu);
end
[~, zeta] = convhulln(mu);


for i = n:nframes
    x = X(i, :); % current point
    
    d = zeros(K, 1);
    div = zeros(K, 1);
    identicalframe = ismember(x, mu, 'rows');
    if identicalframe
        continue
    end
    for k = 1:K
        muik = mu;
        muik(k, :) = x;

        [~, divscore] = convhulln(muik);

        div(k) = divscore;
        d(k) = B*norm(x - mu(k, :)) - i*C*(1 - B)*(divscore - zeta);
    end
    [~, dmin] = min(d);

    if (div(dmin) > zeta) || (randi(100) <= randupdate)
        zeta = div(dmin);
        mu(dmin, :) = x;
        S(dmin) = i;
    end
end
end
