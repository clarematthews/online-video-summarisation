function fnum = kf_choose_the_most_similar(features)
mdata = mean(features, 1);
fnum = knnsearch(features, mdata); % find the closest instance
end