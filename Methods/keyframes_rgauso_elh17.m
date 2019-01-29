function S = keyframes_rgauso_elh17(X, B, lambda)
% B = batch size
% la = regularisation parameter

% -------------------------------------------------------------------------

nframes = size(X, 1);
S = []; % keyframes

for i = 1:B:nframes
    end_batch = min(nframes, i + B - 1);
    batch_index = i:end_batch;
    batch = X(batch_index, :);
    nbatch = size(batch, 1); % actual batch size
    
    % distances from example i to all other examples
    di = squareform(pdist(batch)); 
    [~, l_star] = min(sum(di)); % element with minimum sum of distances
    lambda_max = max(sum(abs(di - repmat(di(:, l_star), 1, nbatch))))/2;
    l = lambda_max*lambda;
    
    % Here starts the part of the algortithm given in the paper
    % [Elhamifar17]
    
    D2 = squareform(pdist(batch)); % distances within the batch
    P = [];
    Q = 1:nbatch;
    if ~isempty(S)
        D1 = pdist2(X(S, :), batch); % distances to selected
    else
        D1 = [];
    end
    for j = 1:nbatch
        
        al = max(0, crit(D1, D2, P, l) - crit(D1, D2, [P j], l));
        bl = max(0, crit(D1, D2, Q, l) - crit(D1, D2, setxor(Q, j), l));

        if (~al && ~bl) || rand < al/(al+bl)
            P = [P j]; %#ok<AGROW>
        else
            Q = setxor(Q, j);
        end
    end
    S = [S batch_index(P)]; %#ok<AGROW>
end
end

function f = crit(D1, D2, candidate, lambda)
% D1 = distances between K elements of sample 2 (D_n) and M elements of
%    sample 1 (Epsilon_0) (M-by-K)
% D2 = distances between K elements of sample 2 and K elements of sample 2 
%    (K-by-K)
% candidate = indices of selected objects from sample 2
% lambda = regularisation parameter

if isempty(D1) && isempty(candidate)
    f = sum(sum(D2));
elseif isempty(D1) && ~isempty(candidate)
    m2 = min(D2(candidate, :), [], 1);
    f = sum(m2) + lambda*numel(candidate);
elseif ~isempty(D1) && isempty(candidate)
    m1 = min(D1, [], 1);
    f = sum(m1);
else
    m1 = min(D1, [], 1);
    m2 = min(D2(candidate, :), [], 1);
    f = sum(min([m1; m2])) + lambda*numel(candidate);
end
end


