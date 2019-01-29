function [mu, cov, clustersize, noise] = define_data(idx)
% mean, covariance, cluster size and noise definition for different data
% examples 
%
% -------------------------------------------------------------------------

switch idx
    
    case 1 % Simple example      
        mu = [0 5; -5 0; 5 0];
        cov = repmat(eye(2), [1, 1, 3]);
        clustersize = 30*ones(1, 3);
        noisedim = 2;
        noisevar = 0.5;
        
    case 2 % Higher variance & closer clusters
        mu = [0 5; -4 0; 3 2];
        cov = repmat(1.5*eye(2), [1, 1, 3]);
        cov(:, :, 2) = 2*eye(2);
        clustersize = 30*ones(1, 3);
        noisedim = 2;
        noisevar = 0.5;
        
    case 3 % More dimensions - 6
        mu = [0 5 0 1 1 6; -5 0 -1 -1 0 3; 5 0 2 2 -2 -2];
        cov = repmat(eye(6), [1, 1, 3]);
        clustersize = 30*ones(1, 3);
        noisedim = 0;
        noisevar = 0;
        
    case 4 % More dimensions - 8
        mu = [0 5 5 0 0 1 1 6; -5 0 3 1 -1 -1 0 3; 5 0 2 -1 2 2 -2 -2];
        cov = repmat(eye(8), [1, 1, 3]);
        clustersize = 50*ones(1, 3);
        noisedim = 0;
        noisevar = 0;

    case 5 % More & variable sized clusters
        mu = [1 5; -5 0; 5 1; 1 -5; -3 3; 4 -3; 0 0];
        cov = repmat(0.5*eye(2), [1, 1, 7]);
        clustersize = [30 10 30 5 10 20 15];
        noisedim = 2;
        noisevar = 0.5;
        
    case 6 % More & variable sized clusters
        mu = [1 5; -5 0; 5 1; 1 -5; 4 -3; 0 0];
        cov = repmat(0.5*eye(2), [1, 1, 6]);
        clustersize = [7 36 51 14 6 20];
        noisedim = 2;
        noisevar = 0.5;
        
    case 7 % Non-symmetric variance & non-zero covariance
        mu = [0 5; -5 0; -3 0; 4 1];
        cov = [0.2 0; 0 2];
        cov(:, :, 2) = [1 1.4; 1.4 2];
        cov(:, :, 3) = [1 1.2; 1.2 2];
        cov(:, :, 4) = [2 -1.4; -1.4 1];
        clustersize = [51, 78, 17, 104];
        noisedim = 2;
        noisevar = 0.5;
     
    case 8 % Repeated cluster
        mu = [0 5; -5 0; 5 0; -5 0];
        cov = repmat(eye(2), [1, 1, 4]);
        clustersize = 30*ones(1, 4);
        noisedim = 2;
        noisevar = 0.5;

    case 9 % Noisy (ego-centric-esque)
        mu = [0 5; 2 2; -6 0; 0 0];
        cov = repmat(eye(2), [1, 1, 4]);
        clustersize = [40, 10, 50, 20];
        noisedim = 0;
        noisevar = 0;
end
noise = [noisedim, noisevar];
end

