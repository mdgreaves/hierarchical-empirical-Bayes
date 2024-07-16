function heb_study(P, C, network, HEB)
% HEB_STUDY Function to implement hierarchical empirical Bayes (HEB) 
% modeling of effective connectivity.
%
% This function is associated with the study: Greaves et al. (2024).
% DOI: https://doi.org/10.1101/2024.04.03.587831
%
% DESCRIPTION: This function explores HEB models of effective connectivity 
% and tests their reliability and out-of-sample validity, following the 
% procedures described in the associated study.
%
% INPUTS:
%   P        -  Cell array of DCMs (which constitute the first level of the 
%               HEB model). 
%   C        -  Normalized [0,1] structural connectivity data used to 
%               inform third-level empirical priors. In principle, this 
%               could be other relevant data.
%   network  -  String, name of the network used for output naming. 
%   HEB      -  (Optional) File path to a pre-existing group-level model
%               structure containing Bayesian model average (BMA) 
%               parameters of a data-to-variance mapping that one wishes to 
%               "validate".
%
% REQUIREMENTS:
%   - SPM12 must be installed and added to the MATLAB path. 
%
% USAGE:
%   heb_study(P, C, network, HEB)
%
% Example:
%   P = {...}; % Define subject-level DCM structures 
%   C = ...;   % Define normalized structural connectivity data 
%   network = 'default_network';
%   heb_study(P, C, network);

% Ensure dimensions of C match those of the network in DCM, and ensure that
% C is in the correct range (i.e., [0, 1]).
assert(all(size(C) == [P{1}.n, P{1}.n]),...
    'Dimensions of C must match the dimensions of the DCM network');
assert(all(C(:) >= 0 & C(:) <= 1), ['To follow procedures ',...
        'in the associated study, Values in C must be within the ',...
        'range [0, 1].']);

% Validate subject-level priors for state transitions
for i = length(P)
    p = P{i}.M.pC(sub2ind(size(P{1}.M.pC), (1:P{1}.n^2)', (1:P{1}.n^2)'));
    assert(all(1./p(~eye(P{1}.n)) == 2), ['To follow procedures ',...
        'in the associated study, the precision of subject-level ',...
        'priors over state transitions, when i ~= j, should equal 2.']);
    assert(all(1./p(~(~eye(P{1}.n))) == 64), ['To follow procedures ',...
        'in the associated study, the precision of subject-level ',...
        'priors over state transitions, when i == j, should equal 64.']);
end

% If no pre-existing group-level model is provided as input, the function 
% proceeds to explore hierarchical empirical Bayes models with third-level 
% empirical priors post hoc, per the face validation phase described in the 
% associated study.

if ~exist('HEB', 'var')
    % Specify parameter structure.

    % Parameter combinations for data-to-variance mapping to explore
    % Define number of points for sampling
    num_points = 30;

    % Generate alpha (intercept) values ranging from -1/2 to 1/2
    params.alphas = linspace(-p(2), p(2), num_points);

    % Initialize beta (slope) matrix
    params.betas = zeros(num_points, num_points);

    % For each alpha, generate beta ensuring that alpha + beta <= 1/2
    for i = 1:num_points
        params.betas(i, :) = linspace(0, p(2) - params.alphas(i),...
            num_points);
    end
    
    % Specify the field(s) to be treated as random effects
    params.field  = 'A';

    % Name of the network (used to save output)
    params.name   = network;

    % Implement the two-step hierarchical empirical Bayes model.
    HEB = heb_explore(P, C, params);

    save(fullfile(pwd,...
        sprintf('HEB_explore_%s.mat', HEB.params.name)), 'HEB', '-v7.3');

% If a pre-existing group-level model is provided as input, the function 
% proceeds to assess the consistency of the model's prior constraints, 
% per the tests of reliability and out-of-sample validity described in the 
% associated study.

elseif exist('HEB', 'var') && isfile(HEB)
    HEB_exp = getField(HEB, 'HEB');

    % Invert hierarchical empirical Bayes model without third-level
    % empirical priors
    [HEB_group_null, HEB_subs_null] = heb_capture(P, struct(),...
        HEB_exp.params.field);

    % Obtain the free energy estimates for each (reduced) DCM resulting 
    % from hierarchical model without third-level empirical priors
    F_subs_null = arrayfun(@(p1) HEB_subs_null{p1, 1}.F,...
        1:length(HEB_subs_null))';

    % Tests of reliability and out-of-sample validity:
    % Obtain prior variance anew using the Bayesian model average (BMA)
    % parameters of the data-to-variance mapping.
    variance = HEB_exp.winning.beta_bma .* C...
        + HEB_exp.winning.alpha_bma;

    % Add data-informed variances to reduced covariance matrix.
    rC = HEB_exp.HEB_null.M.pC;
    indx = find(~eye(length(C)));
    rC(sub2ind(size(HEB_exp.HEB_null.M.pC), indx', indx')) = variance(indx);

    % Add reduced covariances to model structure
    M = HEB_exp.HEB_null.M;
    M.bC = rC;
    M.bE = HEB_exp.HEB_null.M.pE;

    % Invert hierarchical empirical Bayes model with third-level empirical,
    % leveraging the BMA data-to-variance mapping obtained out of session 
    % (or out of sample).
    [HEB_group_val, HEB_subs_val] = heb_capture(P, M,...
        HEB_exp.params.field);

    % Obtain the free energy estimates for each (reduced) DCM resulting 
    % from hierarchical model with third-level empirical priors
    F_subs_val = arrayfun(@(p1) HEB_subs_val{p1, 1}.F,...
        1:length(HEB_subs_val))';

    save(fullfile(pwd,...
        sprintf('HEB_validate_%s.mat', HEB_exp.params.name)),...
        'HEB_group_null', 'HEB_group_val', 'F_subs_null', 'F_subs_val',...
        'HEB_exp', '-v7.3');

    % Report the log Bayes factor for the group-level model
    fprintf(['\nInformed vs. uninformed (group-level) model: ',...
        'log BF %d\n\n'], HEB_group_val.F-HEB_group_null.F);

    % Report the mean log Bayes factor for the subject-level models
    fprintf(['Informed vs. uninformed (subject-level) model: ',...
        'mean log BF: %d\n\n'], mean(F_subs_val-F_subs_null))
end
end

% Helper function for loading files.
function Output = getField(filepath, var)
loadedStruct = load(filepath, var);
fields = fieldnames(loadedStruct);
Output = loadedStruct.(fields{1});
end