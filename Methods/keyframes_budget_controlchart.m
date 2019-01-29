function S = keyframes_budget_controlchart(X, buffer, minseq, thresh,...
    budget, timeweight, varargin)
% X = frame features
% buffer = buffer size
% minseq = minimum segment length
% thresh = baseline keyframe similarity threshold for merging shots
% budget = maximum number of keyframes
% timeweight = weight to give to time in difference function
% -------------------------------------------------------------------------
nframes = size(X, 1);
buffer = min(buffer, nframes);
if size(varargin) > 0
    framepath = varargin{1};
else
    framepath = '';
    % For normalising synthethic data for diff calculations
    maxvals = max(X);
    minvals = min(X);
    datanorm = (X - minvals)./(maxvals - minvals);
end

% -------------------------------------------------------------------------
% Initialise

dt = 1;
S = 1:budget; % Add first frames as keyframes
t = S*dt;

% Calculate pairwise similarities between initial keyframes
simtx = ones(budget);
for i = 1:budget - 1
    for j = i + 1:budget
        if framepath % real videos
            simtx(i, j) = simfns(X, [S(i); S(j)],...
                [t(i); t(j)], nframes, timeweight, framepath);
        else % synthethic data
            simtx(i, j) = simfns(datanorm, [S(i); S(j)],...
                [t(i); t(j)], nframes, timeweight);
        end
        
        simtx(j, i) = simtx(i, j);
    end
end

% -------------------------------------------------------------------------
% Keyframes selected by CC method

if framepath % real videos
    Scc = keyframes_controlchart(X(budget + 1: end, :), buffer, minseq,...
        thresh, framepath);
else % synthetic data
    Scc = keyframes_controlchart(X(budget + 1: end, :), buffer, minseq,...
        thresh);
end

Scc = Scc + budget;

% -------------------------------------------------------------------------
% Process newly selected keyframes - for each keyframe selected by the CC
% method, ignore if too similar to existing keyframes, or else swap with
% the most redundant existing keyframe

for kf = Scc
    simframe = most_sim_frame(simtx); % most redundant current keyframe
    
    % Calculate similarity of potential new keyframe to existing keyframes
    simvec = zeros(1, budget);
    for i = 1:budget
        if framepath % real video
            simvec(i) = simfns(X, [kf; S(i)],...
                [kf, t(i)], nframes, timeweight, framepath);
        else % synthetic data
            datanorm = (X - minvals)./(maxvals - minvals);
            simvec(i) = simfns(datanorm, [kf; S(i)],...
                [kf, t(i)], nframes, timeweight);
        end
    end
    meansim = mean(simvec);
    if meansim < simframe.sim % proposed new keyframe better
        % Replace most redundant existing keyframe with new keyframe and
        % update pairwise similarities for keyframe set
        mostsim = simframe.idx;
        S(mostsim) = kf;
        t(mostsim) = kf;
        simtx(mostsim, 1:mostsim - 1) = simvec(1:mostsim - 1);
        simtx(mostsim, mostsim + 1:end) = simvec(mostsim + 1:end);
        simtx(:, mostsim) = simtx(mostsim, :)';
    end
end

end

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

function frame = most_sim_frame(simtx)
nframes = size(simtx, 1);
% Find most similar pair of keyframes
maxsim = max(max(triu(simtx, 1)));
[i, j] = find(triu(simtx, 1) == maxsim);
i = i(1);
j = j(1);

% Find which one of most similar pair is most similar to others, on average
meani = mean(simtx(i, setdiff(1:nframes, i)));
meanj = mean(simtx(j, setdiff(1:nframes, j)));
if meani > meanj
    mostsim = i;
    meansim = meani;
else
    mostsim = j;
    meansim = meanj;
end
frame.idx = mostsim;
frame.sim = meansim;
end


function value = simfns(data, frames, times, totaltime, timeweight,...
    varargin)
t1 = times(1);
t2 = times(2);
if size(varargin) > 0 % real video
    % Feature similarity
    dx = check_similarity(frames(1), frames(2), varargin{1});
else % synthetic data
    ndim = size(frames, 2);
    x1 = data(frames(1), :);
    x2 = data(frames(2), :);
    % Feature similarity
    dx = norm(x1 - x2)/ndim;
end
% Time similiarity
dt = abs(t1 - t2)/totaltime;

% Combined similarity measure
value = 1 - (timeweight*dt + (1 - timeweight)*dx);
end
