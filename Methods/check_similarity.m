function sc = check_similarity(kf1, kf2, imgdr)
% this function needs to have access to the vid library
% It calculates the similarity score between two key frames, given by their
% frame numbers.
% The similarity score calculates based on De Avila paper, using Histogram
% with 16 bins in HSV feature space.

%--------------------------------------------------------------------------

img1 = vid_extract_chosen_frames(imgdr,kf1,1);
x1 = vid_get_features(img1{1,1},'H',1,16); x1 = x1/sum(x1);

img2 = vid_extract_chosen_frames(imgdr,kf2,1);
x2 = vid_get_features(img2{1,1},'H',1,16); x2 = x2/sum(x2);

sc = sum(abs(x1 - x2));