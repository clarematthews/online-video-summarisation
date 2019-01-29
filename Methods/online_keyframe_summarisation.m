function S = online_keyframe_summarisation(data, method, params)

% Run specified method on data to extract keyframes

% =========================================================================

addpath Methods

% Run summarisation

switch method
    
    % [Abd08] - SVD
    case 'abd'
        if numel(params) ~= 2
            err(method)
        end
        M = params(1); % window size
        tau = params(2); % threshold
        
        S = keyframes_svd_abd08(data, M, tau);
        
        % [Almeida12] - Change detection
    case 'almeida'
        if numel(params) ~= 2
            err(method)
        end
        e = params(1); % threshold for peak detection
        lambda = params(2); % threshold for dropping sequence (%)
        
        S = keyframes_almeida13(data, e, lambda);
        
        % [Anirudh16] - Diversity promotion
    case 'anirudh'
        if numel(params) ~= 4
            err(method)
        end
        
        K = params(1);
        C = params(2);
        B = params(3);
        randupdate = params(4);
        
        S = keyframes_anirudh16(data, K, C, B, randupdate);
        
        % [Elhamifar17] - Submodular optimisation
    case 'elhamifar'
        if numel(params) ~= 2
            err(method)
        end
        B = params(1); % batch size
        la = params(2); % lambda
        
        S = keyframes_rgauso_elh17(data, B, la);
        
        % [Mei15] - Minimum sparse reconstruction
    case 'mei'
        if numel(params) ~= 1 && numel(params) ~= 2
            err(method)
        end
        Tpor = params(1);
        if numel(params) == 2
            maxframes = min(params(2), size(data, 2));
        else
            maxframes = size(data, 2);
        end
        
        S = keyframes_msr_mei15(data, Tpor, maxframes);
        
        % [Ou14] - GMM
    case 'ou'
        if numel(params) ~= 5
            err(method)
        end
        K = params(1);
        a = params(2);
        T = params(3);
        var0 = params(4); % "large" variance for new component
        weight0 = params(5); % "small" weight for new component
        
        S = keyframes_gmm_ou14(data, K, a, T, var0, weight0);
        
        % [Rasheed03] - Change detection
    case 'rasheed'
        if numel(params) ~= 2
            err(method)
        end
        Tcolour = params(1); % shot boundary threshold;
        Th = params(2); % keyframe threshold;
        
        % Require feature vectors to be positive and normalised
        datashift = data - min(data(:));
        datanorm = datashift./repmat(sum(datashift, 2), 1, size(data, 2));
        
        S = keyframes_rasheed03(datanorm, Tcolour, Th);
        
        % [Song05] - GMM merging
    case 'song'
        if numel(params) ~= 2
            err(method)
        end
        M = params(1); % batch size
        alpha = params(2); % significance level for merging components
        
        S = keyframes_gmm_song05(data, M, alpha);
        
        % [Truong07] - Significant content change
    case 'truong'
        if numel(params) ~= 2
            err(method)
        end
        e = params(1); % threshold for change
        fn = params(2); % distance fn 1: Euclidean, 2: Minkowski, 3: Cosine
        
        S = keyframes_scc_truong07(data, e, fn);
        
        % Shewhart - Control chart
    case 'shewhart'
        if numel(params) ~= 3 && numel(params) ~= 4
            err(method)
        end
        
        if numel(params) == 3
            minseg = params(1); % minimum segment length
            bsize = params(2); % buffer size
            thresh = params(3); % threshold for keyframe similarity
            S = keyframes_controlchart(data, minseg, bsize, thresh);
        else
            minseg = params{1}; % minimum segment length
            bsize = params{2}; % buffer size
            thresh = params{3}; % threshold for keyframe similarity
            framepath = params{4};
            S = keyframes_controlchart(data, minseg, bsize, thresh,...
                framepath);
        end
        
        % Control chart with budget
    case 'budget-cc'
        if numel(params) ~= 5 && numel(params) ~= 6
            err(method)
        end
        
        if numel(params) == 3
            minseg = params(1); % minimum segment length
            bsize = params(2); % buffer size
            thresh = params(3); % threshold for keyframe similarity
            budget = params(4); % total number of keyframes
            timeweight = params(5); % importance of time in difference fn
            S = keyframes_budget_controlchart(data, minseg, bsize,...
                thresh, budget, timeweight);
        else
            minseg = params{1}; % minimum segment length
            bsize = params{2}; % buffer size
            thresh = params{3}; % threshold for keyframe similarity
            budget = params{4}; % total number of keyframes
            timeweight = params{5}; % importance of time in difference fn
            framepath = params{6};
            S = keyframes_budget_controlchart(data, minseg, bsize,...
                thresh, budget, timeweight, framepath);
        end
        
        % Control chart with dynamic threshold parameter
    case 'dynamic-cc'
        if numel(params) ~= 4 && numel(params) ~= 5
            err(method)
        end
        
        if numel(params) == 3
            minseg = params(1); % minimum segment length
            bsize = params(2); % buffer size
            thresh = params(3); % threshold for keyframe similarity
            budget = params(4); % maximum number of keyframes
            S = keyframes_dynamic_controlchart(data, minseg, bsize,...
                thresh, budget);
        else
            minseg = params{1}; % minimum segment length
            bsize = params{2}; % buffer size
            thresh = params{3}; % threshold for keyframe similarity
            budget = params{4}; % total number of keyframes
            framepath = params{5};
            S = keyframes_dynamic_controlchart(data, minseg, bsize,...
                thresh, budget, framepath);
        end
        
        % Control chart with dynamic threshold parameter
    case 'dynamic-replace-cc'
        if numel(params) ~= 4 && numel(params) ~= 5
            err(method)
        end
        
        if numel(params) == 3
            minseg = params(1); % minimum segment length
            bsize = params(2); % buffer size
            thresh = params(3); % threshold for keyframe similarity
            budget = params(4); % maximum number of keyframes
            S = keyframes_dynamic_replace_controlchart(data, minseg,...
                bsize, thresh, budget);
        else
            minseg = params{1}; % minimum segment length
            bsize = params{2}; % buffer size
            thresh = params{3}; % threshold for keyframe similarity
            budget = params{4}; % total number of keyframes
            framepath = params{5};
            S = keyframes_dynamic_replace_controlchart(data, minseg,...
                bsize, thresh, budget, framepath);
        end
        
    otherwise
        error(['No match found for method ' method]);
end
end

function err(method)
error(['Incorrect number of parameters for ' method])
end
