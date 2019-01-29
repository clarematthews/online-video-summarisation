function S = keyframes_controlchart(X, buffer, minseq, thresh, varargin)
% X = frame features
% buffer = buffer size
% minseq = minimum segment length
% thresh = keyframe similarity threshold for merging shots
 
% -------------------------------------------------------------------------
 
nframes = size(X, 1);
buffer = min(buffer, nframes);
if size(varargin) > 0
    framepath = varargin{1};
else
    framepath = '';
end
 
% Initialisation

curridx = 1:buffer;
S = []; % key frame set
maxvals = max(X);
minvals = min(X);

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
        lastidx = curridx;
    else
        % Check keyframe similarity
        lastkf = S(end);
        if framepath
            kfsim = check_similarity(lastkf, kfnum, framepath);
        else
            lastnorm = (X(lastkf, :) - minvals)./(maxvals - minvals);
            currnorm = (X(kfnum, :) - minvals)./(maxvals - minvals);
            kfsim = sum(abs(lastnorm - currnorm));
        end
        if kfsim < thresh % keyframes too similar
            % Concatenate two shot boundaries
            S(end) = [];
            combidx = [lastidx curridx];
            kfnum = combined_shot_kf(X, combidx);
            lastidx = combidx;
        else
            lastidx = curridx;
        end
    end
    S = [S kfnum]; %#ok<AGROW>
    
    % Reset current shot
    curridx = i;
    
end

% Include the last shot
if numel(curridx) >= minseq
    shot = X(curridx, :);
    kfnum = knnsearch(shot, mean(shot, 1));
    kfnum = curridx(kfnum);
    if ~isempty(S)
        lastkf = S(end);
        if framepath
            kfsim = check_similarity(lastkf, kfnum, framepath);
        else
            lastnorm = (X(lastkf, :) - minvals)./(maxvals - minvals);
            currnorm = (X(kfnum, :) - minvals)./(maxvals - minvals); 
            kfsim = sum(abs(lastnorm - currnorm));
        end
        
        % Two consecutive kfs are similar
        if kfsim < thresh
            % Concatenate two shot boundaries
            S(end) = [];
            combidx = [lastidx curridx];
            kfnum = combined_shot_kf(X, combidx);
        end
    end
    
    S = [S kfnum];
end

end

function kf = combined_shot_kf(data, idx)
shot = data(idx, :);
% Recalculate keyframe
newkf = knnsearch(shot, mean(shot, 1));
kf = idx(newkf);
end

