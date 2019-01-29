function S = keyframes_rasheed03(X, Tcolour, Th)
% Tcolour = threshold for shot boundaries
% Th = threshold for keyrame selection

% -------------------------------------------------------------------------

nframes = size(X, 1);

% Initialisation

S = []; % keyframes

framesi = X(1, :);
shotst = 1;

for i = 2:nframes
    x = X(i, :); % current point

    lastframe = framesi(end, :);
    D = sum(bsxfun(@min, lastframe, x));
    
    if D < Tcolour % new shot
        keyframes = shotkeyframes(framesi, Th);
        S = [S keyframes + shotst - 1]; %#ok<AGROW>
        
        framesi = x;
        shotst = i;
        continue
    end
    
    framesi = [framesi; x]; %#ok<AGROW>
    
end

keyframes = shotkeyframes(framesi, Th);
S = [S keyframes + shotst - 1];

end

function S = shotkeyframes(shotframes, thresh)
nframes = size(shotframes, 1);
S = ceil(nframes/2);
keyframes = shotframes(S, :);

for i = setxor(1:nframes, S)
    Dmax = 0;
    x = shotframes(i, :);
    for k = 1:numel(S)
        Dk = sum(bsxfun(@min, keyframes(k, :), x));
        Dmax = max(Dmax, Dk);
    end
    if Dmax < thresh
        keyframes = [keyframes; x]; %#ok<AGROW>
        S = [S i]; %#ok<AGROW>
    end
end
end