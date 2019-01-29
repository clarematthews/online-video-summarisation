function S = keyframes_scc_truong07(X, e, fnc)
% X = data
% e = epsilon (threshold)
% fnc = distance function: 1: Euclidean, 2: Minkowski, 3: Cosine

% -------------------------------------------------------------------------
nframes = size(X, 1);
S = 1; % initial frame (always?)

for i = 2:nframes
    Xi = X(i, :); % current frame
    candidate = find_next_frame(X(S(end), :), Xi, fnc, e);
    if candidate % add the frame
        S = [S, i];  %#ok<AGROW>
    end
end
end

function candidate = find_next_frame(x, y, fnc, e)
switch fnc
    case 1
        di = pdist2(x, y);
    case 2
        di = pdist2(x, y, 'minkowski', 1);
    case 3
        di = pdist2(x, y, 'cosine');
end
candidate = di > e;
end



