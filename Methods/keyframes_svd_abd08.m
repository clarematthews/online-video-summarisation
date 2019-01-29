function S = keyframes_svd_abd08(X, M, tau)
% M = window size
% tau = threshold

% -------------------------------------------------------------------------

minrank = 1;

nframes = size(X, 1);
S = []; % keyframes
Q = []; % shot frames

R = zeros(1, nframes);
for i = M:nframes
    window = X(i - M + 1:i, :);
    [~, si] = svd(window);
    v = diag(si);
    rank = sum(v/max(v) >= tau);
    R(i) = rank;
        
    if i > M && R(i - 1) == minrank && R(i) > minrank % i - 1 = end of shot
        if isempty(Q)
            % add frame in middle of opening shot as keyframe
            keyframe = round(i/2);
        else
            % add shot frame with greatest rank to keyframes
            [~, ind] = max(R(Q));
            keyframe = Q(ind);
        end        
        S = [S keyframe];  %#ok<AGROW>
        Q = i;
    elseif rank > minrank
        % shot buid-up
        Q = [Q i];  %#ok<AGROW>
    end
end

if ~isempty(Q) % get keyframe from final shot
    [~, ind] = max(R(Q));
    S = [S Q(ind)];
end



