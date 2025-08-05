function heb_sim(snr)
% HEB_SIM - Demonstrates structurally informed methods for characterizing  
%           directed connectivity.  
%  
% Syntax: heb_sim(snr)  
%  
% Description:  
%   This function simulates subject-level effective connectivity matrices  
%   as perturbations of a ground-truth (group-level) matrix. It compares  
%   two methods for estimating structurally informed directed
%   connectivity:  
%     1. Multivariate autoregressive (MVAR) modeling using structural  
%        connectivity (SC) as a mask.  
%     2. Hierarchical empirical Bayes (HEB) applied to dynamic causal  
%        models (DCMs).  
%  
%   The simulation pipeline includes:  
%     - Generation of ground-truth group-level structural and effective  
%       connectivity matrices (`heb_connectivity`).  
%     - Subject-level variability with Gaussian noise (`heb_error`).  
%     - MVAR estimation (`heb_init_mvar`).  
%     - DCM inversion and HEB updates (`heb_init_dcm`, `heb_study`).  
%  
% Input:  
%   snr - Desired signal-to-noise ratio, controlling the level of noise  
%         added to simulated signals.  
%  
% Example:  
%   heb_sim(1); % Runs the simulation with an SNR of 1  
%  
% See also: heb_connectivity, heb_error, heb_init_mvar, heb_init_dcm,  
%           heb_study  

% Simulate ground truth
%--------------------------------------------------------------------------

% Add software to path (ensure SPM is on the path)
assert(exist('spm', 'file') == 2, 'SPM is not on the MATLAB path.');
addpath('../core/');

% Set random seed for reproducibility
rng('default');

% Options
dlt = true;  % Remove default output from core functions

% Global parameters
n   = 6;     % Number of regions
s   = 50;    % Number of instantiations (simulated subjects)
sp  = 1/3;   % Sparsity of structural connectivity
T   = 1e3;   % Number of observations (scans)
TR  = 1;     % Repetition time
tau = ...    % Transit times
    randn(n,1)*exp(-4);

% Group-level parameters
a   = 1/16;  % Baseline variance of off-diagonal effective connectivity
b   = 1/10;  % Sensitivity of variance to structural connectivity
v   = 1/64;  % Variance of diagonal effective connectivity

% Subject-level parameters
ve  = v;     % Variance of off-diagonal effective connectivity
vi  = v;     % Variance of diagonal effective connectivity
pv  = 1/2;   % Prior variance of off-diagonal effective connectivity

% Filename
network = ...
    sprintf('sim_n%d_s%d_sp%d_snr%d_T%d_TR%d_a%d_b%d_v%d',...
    n, s, 1/sp, snr, T, TR, 1/a, 1/b, 1/v);

% Generate ground-truth effective connectivity
[Ag, SC]    = heb_connectivity(n, a, b, (1-sp), v, T, TR, tau);
[Y, As]     = heb_error(Ag, s, ve, vi, T, TR, snr, tau); %#ok

% Fit/invert models
%--------------------------------------------------------------------------
% MVAR
% Obtain (structurally informed) autoregressive parameters
W   = nan(n,n,s);
for i = 1:s
    [W(:,:,i), ~] = heb_init_mvar(Y{i}', SC);
end

% DCM
% Infer (non-informed) subject-level effective connectivity
P = cell(s,1);

% Start parallel pool
parpool('local');

parfor i = 1:s
    fprintf('\n\nInverting instantiation: %.0f\n', i);
    DCM = heb_init_dcm(Y{i}, pv, TR);
    P{i} = DCM{1};
end
delete([P{s}.name, '.mat'])

% Shut down the parallel pool
delete(gcp('nocreate'));

% HEB (explore)
% Invert (structurally informed) effective connectivity
heb_study(P, SC, network);
EXP = load(dir(fullfile(pwd, '**',...
    sprintf('*explore_%s*.mat', network))).name, 'HEB').HEB;
RFX = EXP.HEB_null;

% HEB (step 1: apply BMA prior-variance transformation)
rV = EXP.winning.beta_bma .* SC + EXP.winning.alpha_bma;

% Add updated variances to reduced covariance matrix
pE = RFX.M.pE;                  pC = RFX.M.pC;
Ep = RFX.Ep;                    Cp = RFX.Cp;
rE = pE;                        rC = pC;
indx = find(~eye(n));
rC(sub2ind(size(rC), indx', indx')) = rV(indx);

% Reduce group-level model
HEB = struct;
[~, HEB.Ep, HEB.Cp] = spm_log_evidence_reduce(Ep, Cp, pE, pC, rE, rC);

% HEB (step 2: re-evaluate DCMs)
DCM = cell(s,1);
for i = 1:s
    rE = HEB.Ep;                rC = HEB.Cp;
    Ep = sparse(P{i}.Ep.A);     Cp = full(P{1}.Cp(1:n^2,1:n^2));
    [~, DCM{i}.Ep, DCM{i}.Cp] = spm_log_evidence_reduce(Ep, Cp, pE, pC,...
        rE, rC);
end

% Delete and save files
%--------------------------------------------------------------------------
if dlt == true
    delete(fullfile(pwd, dir(fullfile(pwd, '**',...
        sprintf('*explore_%s*.mat', network))).name));
end

if ~exist(fullfile(pwd, 'output'), 'dir')
    mkdir('output')
end
save(fullfile('./output', network), '-v7.3');

end
