function S = keyframes_msr_mei15(X, Tpor, maxframes)
% Tpor = threshold for frame selection
% maxframes = maximum number of frames in reconstruction set
 
% -------------------------------------------------------------------------
 
[nframes, ndim] = size(X);
maxframes = min(maxframes, ndim);
 
% Initialisation
 
S = 1;
F = X(1, :);
R = F;
 
for i = 2:nframes
    nkeyframes = size(F, 1);
    x = X(i, :); % current point
     
    rvec = R'*((R*R')\R)*x';
    rnorm = norm(rvec);
    por = rnorm/norm(x);
 
    if por < Tpor
        F = [F; x]; %#ok<AGROW>
        S = [S i]; %#ok<AGROW>
        if nkeyframes < maxframes % keep all keyframes for reconstruction
            R = F;
        else
            errors = zeros(nkeyframes, 1);
            for f = 1:nkeyframes
                errors(f) = recon_error(F(f, :), F);
            end
            [~, idx] = sort(errors, 'descend');
            R = F(idx(1:maxframes - 1), :);
            R = [R; x]; %#ok<AGROW>
        end
    end
 
end
end
 
 
function error = recon_error(frame, keyframes)
projector = eye(size(frame, 2)) - frame'*((frame*frame')\frame);
error = norm(projector*keyframes');
end