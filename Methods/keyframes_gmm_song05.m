function S = keyframes_gmm_song05(X, M, alpha)
% M = frame batch size
% alpha = significance level for merging

% -------------------------------------------------------------------------

[nframes, ndim] = size(X);

% The initialisation is not specified in [Song05], assume just estimate GMM
% for first batch

M = min(M, nframes);
Xi = X(1:M, :); % take the first M objects as the initial batch

[~, model] = gmmest(Xi'); % components for initial batch
gmu = model.mu;
gcov = model.Sigma;
gweights = model.w;
Kg = size(gmu, 2);

dists = pdist2(gmu', Xi);
[dmin, idx] = min(dists, [], 2);
S = idx';
dstar = dmin';

for i = M + 1:M:nframes
    Xi = X(i:min(i + M - 1, nframes), :);
    
    [alabels, model] = gmmest(Xi'); % components for new batch
    amu = model.mu;
    acov = model.Sigma;
    Ka = size(amu, 2);
    
    % Check if any new components are equivalent to existing components
    equivcomps = [];
    for k = 1:Ka
        Xk = Xi(alabels == k, :);
        equivcompsk = [];
        for j = 1:Kg
            match = comparecomp(Xk, gmu(:, j)', gcov(:, :, j), alpha);
            if match
                equivcompsk = [equivcompsk; j]; %#ok<AGROW>
            end
        end
        % If equivalent to more than one existing component, get best match
        if numel(equivcompsk) > 1
            maxllh = -Inf;
            bestj = 0;
            for jidx = 1:numel(equivcompsk)
                j = equivcompsk(jidx);
                llh = sum(logmvnpdf(Xk, gmu(:, j)', gcov(:, :, j)));
                if llh > maxllh
                    maxllh = llh;
                    bestj = j;
                end
            end
            equivcompsk = bestj;
        end
        if ~isempty(equivcompsk)
            equivcomps = [equivcomps; k, equivcompsk(1)]; %#ok<AGROW>
        end
    end
    
    % Update components
    nmerge = size(equivcomps, 1);
    newcomps = 1:Ka; % new components with no equivalents
    oldcomps = 1:Kg; % existing components to remain unchanged
    if ~isempty(equivcomps)
        newcomps = setxor(newcomps, equivcomps(:, 1));
        oldcomps = setxor(oldcomps, equivcomps(:, 2));
    end
    
    ncluster = Kg + Ka - nmerge;
    mustar = zeros(ndim, ncluster);
    covstar = zeros(ndim, ndim, ncluster);
    weightsstar = zeros(1, ncluster);
    alabelsstar = zeros(1, numel(newcomps));
    for j = 1:Kg
        if ismember(j, oldcomps) % component unchanged (except for weight)
            mustar(:, j) = gmu(:, j);
            covstar(:, :, j) = gcov(:, :, j);
            nj = (i - 1)*gweights(j);
            weightsstar(j) = nj/(i - 1 + M);          
        else % merge new component with equivalent existing component
            k = equivcomps(1, equivcomps(:, 2) == j);
            nj = (i - 1)*gweights(j);
            nk = sum(alabels == k);
            mu = (nj*gmu(:, j) + nk*amu(:, k))/(nj + nk);
            sigma = (nj*gcov(:, :, j) + nk*acov(:, :, k))/(nj + nk) +...
                (nj*gmu(:, j)*gmu(:, j)' +...
                nk*amu(:, k)*amu(:, k)')/(nj + nk) + mu*mu';
            weight = (nj + nk)/(i - 1 + M);
            mustar(:, j) = mu;
            covstar(:, :, j) = sigma;
            weightsstar(j) = weight;
            alabelsstar(alabels == k) = j;
        end        
    end
    count = 1;
    for knew = 1:numel(newcomps) % new components
        k = newcomps(knew);
        mustar(:, Kg + count) = amu(:, k);
        covstar(:, :, Kg + count) = acov(:, :, k);
        nk = sum(alabels == k);
        weightsstar(Kg + count) = nk/(i - 1 + M);
        alabelsstar(alabels == k) = Kg + count;
        count = count + 1;
    end
    
    gmu = mustar;
    gcov = covstar;
    gweights = weightsstar;
    Kg = size(gmu, 2);
            
    % Select keyframe for component if better representative than current
    % component keyframe
    dstar = [dstar ones(1, ncluster - numel(dstar))*Inf]; %#ok<AGROW>
    dists = pdist2(gmu', Xi);
    [dmin, idx] = min(dists, [], 2);
    for k = 1:ncluster
        if dmin(k) < dstar(k)
            dstar(k) = dmin(k);
            S(k) = i + idx(k) - 1;
        end
    end
end
end


function [labels, model] = gmmest(X)
[ndim, M] = size(X);
Ks = 1:5;
minbic = Inf;
for i = 1:numel(Ks)
    [labeli, modeli, llhi] = mixGaussEm(X, Ks(i));
    bici = bic(M, ndim, Ks(i), llhi(end));
    if bici < minbic
        minbic = bici;
        labels = labeli;
        model = modeli;
    end
end
end


function val = bic(n, ndim, k, llh)
nu = k*(ndim + 1)*(ndim + 2)/2 - 1;
val = log(n)*nu - 2*llh;
end


function match = comparecomp(X, mu0, sigma0, alpha)
match = false;
equivcov = w(X, sigma0, alpha);
if equivcov
    equivmu = t2(X, mu0, alpha);
    if equivmu
        match = true;
    end
end
end


function match = w(X, Sigma0, alpha)
[nsample, ndim] = size(X);
L0 = chol(Sigma0, 'lower');
Y = L0\X';
Sy = cov(Y');
W = 1/ndim*trace((Sy - eye(ndim))^2) -...
    ndim/nsample*(1/ndim*trace(Sy))^2 + ndim/nsample;
testval = W*nsample*ndim/2;
chisq = chi2inv(1 - alpha, ndim*(ndim + 1)/2);
match = false;
if testval <= chisq
    match = true;
end
end


function match = t2(X, mu0, alpha)
[nsample, ndim] = size(X);
mux = mean(X);
Sx = cov(X);
T2 = nsample*(mux - mu0)*(Sx\(mux - mu0)');
testval = (nsample - ndim)/(ndim*(nsample - 1))*T2;
f = finv(1 - alpha, ndim, nsample - ndim);
match = false;
if testval <= f
    match = true;
end
end


function [logp] = logmvnpdf(x, mu, Sigma)
% outputs log likelihood array for observations x  where x_n ~ N(mu, Sigma)
% x is NxD, mu is 1xD, Sigma is DxD

[~, ndim] = size(x);
const = -0.5 * ndim * log(2*pi);

xc = bsxfun(@minus, x, mu);

U = chol(Sigma);
y = 2*sum(log(diag(U)));

term1 = -0.5 * sum((xc / Sigma) .* xc, 2);
term2 = const - 0.5 * y;
logp = term1' + term2;

end
