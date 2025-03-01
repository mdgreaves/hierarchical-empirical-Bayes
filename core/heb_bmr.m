function [F,sE,sC] = heb_bmr(alpha, beta, C, P)
% HEB_BMR Function to perform data-guided (or constrained) Bayesian model 
% reduction.
%
% DESCRIPTION:
% This function performs Bayesian model reduction (BMR) using a linear 
% mapping function to update prior variances based on input parameters. It 
% then returns the log-evidence for the reduced model, in addition to the 
% posterior expectations and covariance.
%
% INPUTS:
%   alpha   - Scalar value representing the intercept in a linear mapping.
%   beta    - Scalar value representing the slope in a linear mapping.
%   C       - Normalized [0,1] structural connectivity data used to inform 
%             third-level empirical priors. In principle, this could be 
%             other relevant data.
%   P       - PEB or DCM model structure, including fields:
%               .Ep - posterior expectations,
%               .Cp - posterior covariance,
%               .M.pE - prior expectations,
%               .M.pC - prior covariance.
%
% OUTPUTS:
%   F       - Log model evidence for the reduced model.
%   sE      - Posterior expectations.
%   sC      - Posterior covariance.
%
% REQUIREMENTS:
%   - SPM12 must be installed and added to the MATLAB path.
%
% USAGE:
%   [F, sE, sC] = heb_bmr(alpha, beta, C, P)
%
% Example:
%   alpha = 0.1; % Intercept
%   beta = 0.4;  % Slope
%   C = ...;     % Define normalized structural connectivity data
%   P = ...;     % Define PEB or DCM model structure
%   [F, sE, sC] = heb_bmr(alpha, beta, C, P);

% Apply (linear) mapping function to data to obtain prior variances
variance = beta .* C + alpha;

% Add updated variances to reduced covariance matrix
rC = P.M.pC; 
indx = find(~eye(length(C)));
rC(sub2ind(size(P.M.pC), indx', indx')) = variance(indx);

% Perform Bayesian model reduction and return sufficient statistics for 
% posteriors and the relative log model evidence
[F,sE,sC] = spm_log_evidence_reduce(P.Ep, P.Cp, P.M.pE, P.M.pC,...
    P.M.pE, rC);

end