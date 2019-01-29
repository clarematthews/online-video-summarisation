function S = keyframes_almeida13(X, e, lambda)
% e = threshold for peak detection
% lambda = duration threshold for keeping a sequence

% This function is based on online video summarisation- keyframe which is
% proposed by Jurandy Almeida 2013 and 2012 
% VISION: Video Summarisation for Online applications 2012
% Online video summarisation on compressed domain 2013

% ++++++++++++++++++++++++++++++++++++++++++++++++++++
% we can use the same with xcorr
% ZNCC(1) = xcorr(X(1,:),X(1,:),0,'coeff');
% template = X(1,:);
% for i= 2:d(1)
%     ZNCC(i) = xcorr(X(i,:),template,0,'coeff');
%     template = X(i,:);
% end
% ++++++++++++++++++++++++++++++++++++++++++++++++++++
% Jurandy's formula
% ZNCC is between -1 to 1. Therefore, 1-ZNCC is between 0 and 2
% -------------------------------------------------------------------------

thresh = 2*(1 - e);

nframes = size(X, 1);
template = X(1, :) - mean(X);
peaks = [];
ZNCC = zeros(nframes, 1);

S =[]; % keyframes

for i = 1:nframes
    X_tem = X(i, :) - mean(X);
    ZNCC(i) = sum(X_tem.*template)/(sqrt(sum(template.^2)*sum(X_tem.^2)));
    template = X_tem;
    if i > 3 && i <= nframes
        % The original equation is (ZNCC(i)/ZNCC(i-1))+(ZNCC(i)/ZNCC(i+1))
        % We put i = i-1 to overcome the i+1 argument
        thr = (ZNCC(i - 1)/ZNCC(i - 2))+(ZNCC(i - 1)/ZNCC(i));
        % 0.05 < E < 0.15  based on the papers
        % there is no upper limit for this equation however, if we plot the
        % thr(i) the upper limit would be about 2!
        if thr > thresh
            peaks = [peaks; i - 1]; %#ok<AGROW>
        end
    end
end
% OTHER CONCERN??????
% if we know the total number of frames the algorithm is not really online
% 0.5% * total_number_objects < lamda < 2% * total_number_objects

% For labelling, assume peak represents the start of a group

stframe = 1;
for k = 1:numel(peaks)
    endframe = peaks(k);
    if endframe - stframe < lambda*nframes/100
        stframe = endframe;
        continue
    end
    keyframe = round((stframe - 1 + endframe)/2);
    S = [S keyframe]; %#ok<AGROW>
    stframe = endframe;
end
if (nframes - stframe) >= lambda*nframes/100
    keyframe = round((stframe - 1 + nframes)/2);
    S = [S keyframe];
end
if isempty(S)
    S = round(nframes/2);
end
end