function S = keyframes_gmm_ou14(X, K, a, T, var0, weight0)
% K = number of components
% a = learning rate alpha
% T = threshold for frame selection

% -------------------------------------------------------------------------

thresh = 3;
sigma = 1;

[nframes, ndim] = size(X);
S = []; % keyframes

% The initialisation is not specified in [Ou14]

mu = X(1:K, :); % take the first K objects as the cluster centres

sig2 = ones(size(mu, 1), 1)*sigma;
w = ones(1,K)/K; % equal weights

for i = K + 1:nframes
    x = X(i, :); % current point
    delta = pdist2(x, mu)./sqrt(sig2');
    [minDelta, jstar] = min(delta);
    if minDelta < sigma*thresh
        % Found closest component
        rho = a*mvnpdf(x, mu(jstar, :), eye(ndim)*sig2(jstar));
        mu(jstar, :) = (1 - rho) * mu(jstar, :) + rho * x;
        sig2(jstar) = (1 - rho) * sig2(jstar) + rho *...
            (x - mu(jstar, :)) * (x - mu(jstar, :))';
        w(jstar) = (1 - a) * w(jstar) + a;
        w(setxor(1:K, jstar)) = (1 - a) * w(setxor(1:K, jstar));
        
        % Check whether we should keep the frame
        [~, li] = sort(w./sqrt(sig2'), 'descend');
        k = find(li == jstar); % identify position of winning component
        if sum(w(li(1:k)) >= T) % keep the frame          
            S = [S i]; %#ok<AGROW>
        end
    else
        % 'large' variance and 'small' weight -- whatever this means
        [~, jnotstar] = max(delta);
        mu(jnotstar, :) = x;
        sig2(jnotstar) = var0;
        w(jnotstar) = weight0;
        w = w/sum(w);
    end
end
end




