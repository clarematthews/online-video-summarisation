function S = keyframes_dynamic_controlchart(X, buffer, minseq, thresh,...
    budget, varargin)
% X = frame features
% buffer = buffer size
% minseq = minimum segment length
% thresh = baseline keyframe similarity threshold for merging shots
% budget = maximum number of keyframes
 
% -------------------------------------------------------------------------
 
[nframes, ndim] = size(X);
buffer = min(buffer, nframes);
if size(varargin) > 0
    framepath = varargin{1};
    simdim = 16; % Dimensions of feature space used for similiarity check
else
    framepath = '';
    % For normalising synthetic data
    maxvals = max(X);
    minvals = min(X);
    simdim = ndim; % Dimensions of feature space used for similiarity check
end
 
% Initialisation

curridx = 1:buffer;
S = []; % key frame set


% Calculate mean and std of buffer frames

pairdis = diag(pdist2(X(1:buffer - 1, :), X(2:buffer, :)));
mu = mean(pairdis); 
sig = std(pairdis);
meancount = buffer - 1;
dsumsq = sum(pairdis.^2);

% -------------------------------------------------------------------------
% Process frames

for i = buffer + 1:nframes
    
    % Pairwise distance between frames
    d = pdist2(X(i, :), X(i - 1, :));
    
    % Check the value of pairwise distances
    if d <= mu + 3*sig % no new shot detected
        % Build up current shot
        curridx = [curridx i]; %#ok<AGROW>
        
        % Update mean and std
        meancount = meancount + 1;
        mu = (mu*(meancount - 1) + d)/meancount;
        dsumsq = dsumsq + d^2;
        sig = sqrt((dsumsq/meancount - mu^2)*meancount/(meancount - 1));
        continue
    end
    
    % New shot analysis
    if numel(curridx) < minseq % shot is too short
        % Reset current shot
        curridx = i;
        continue
    end
        
    % Shot keyframe selection
    shot = X(curridx, :);
    kfnum = knnsearch(shot, mean(shot, 1));
    kfnum = curridx(kfnum);
    if isempty(S)
        S = kfnum;
    else
        % Check keyframe similarity
        lastkf = S(end);
        if framepath % real video
            kfsim = check_similarity(lastkf, kfnum, framepath);
        else % synthetic data
            lastnorm = (X(lastkf, :) - minvals)./(maxvals - minvals);
            currnorm = (X(kfnum, :) - minvals)./(maxvals - minvals);
            kfsim = sum(abs(lastnorm - currnorm));
        end
        dynthresh = dynamic_threshold(thresh, numel(S), budget, i,...
            nframes, simdim);
        if kfsim > dynthresh % include keyframe
            S = [S kfnum]; %#ok<AGROW>
        end
    end
        
    % Reset current shot
    curridx = i;
    
end

% Include the last shot
if numel(curridx) >= minseq
    shot = X(curridx, :);
    kfnum = knnsearch(shot, mean(shot, 1));
    kfnum = curridx(kfnum);
    if isempty(S)
        S = kfnum;
    else
        lastkf = S(end);
        if framepath % real video
            kfsim = check_similarity(lastkf, kfnum, framepath);
        else % synthetic data
            lastnorm = (X(lastkf, :) - minvals)./(maxvals - minvals);
            currnorm = (X(kfnum, :) - minvals)./(maxvals - minvals);
            kfsim = sum(abs(lastnorm - currnorm));
        end
        
        % Two consecutive kfs are similar
        dynthresh = dynamic_threshold(thresh, numel(S), budget,...
            nframes, nframes, simdim);
        if kfsim > dynthresh
            S = [S kfnum];
        end
    end   
    
end

end

function thresh = dynamic_threshold(base, numkf, budget, time,...
    totaltime, ndim)
expkf = budget*time/totaltime;
if expkf == budget
    thresh = (numkf >= budget)*ndim;
else
    thresh = (base*(budget - numkf) + ndim*(numkf - expkf))/...
        (budget - expkf);
end
end

