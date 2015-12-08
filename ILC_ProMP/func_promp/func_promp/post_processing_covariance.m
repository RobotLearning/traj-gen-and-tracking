function [covOut, cholOut] = post_processing_covariance( covIn, nBasis, nDemo)
% Numerical treatment on the covariance.
% If not done sometimes covariance values are huge.


   %% bound the covariance
   
   % Using canonical arguments   
   minCov  = zeros(1, nBasis);
   minCorr = 1;

   covTemp1 = boundCovariance( covIn, minCov,  minCorr);

   
   
   %% regularize the covariance
   
   % Using canonical arguments
   priorCovariance = eye(size(covTemp1,1));
   priorCovarianceWeight = 1e-16;
   
   
   %> regularize covariance
   [covOut, cholOut] = regularizeCovariance( covTemp1, ...
                                  priorCovariance, nDemo,...
                                  priorCovarianceWeight);

end
