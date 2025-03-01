function P = heb_init_dcm(y, v, dt)
% heb_init_dcm - Initializes and inverts a DCM for fMRI data using SPM.
%
% Syntax:  P = heb_init_dcm(y, v, dt)
%
% Inputs:
%   y  - Time-series (T x n), where T is time points and n is regions.
%   v  - Variance of Gaussian priors for interregional connectivity.
%   dt - Repetition time (seconds) of the fMRI acquisition.
%
% Outputs:
%   P  - Posterior distribution and model structure.
%
% Description:
%   This function initializes a dynamic causal model (DCM) structure for
%   fMRI data, specifying prior covariance, endogenous connections, and
%   analysis options. The model is then inverted using `spm_dcm_fit` to 
%   estimate effective connectivity. Based on SPM's `spm_dcm_fmri_csd`.
%
% Example:
%   P = heb_init_dcm(timeSeries, 1/2, 2);
%
% See also: spm_dcm_fmri_csd, spm_dcm_fmri_priors, spm_dcm_fit

% DCM structure
% -------------------------------------------------------------------------
DCM.v       = size(y,1);                % number of scans
DCM.n       = size(y,2);                % number of regions
DCM.a       = ones(DCM.n,DCM.n);        % switch on endogenous connections
DCM.b       = zeros(DCM.n,DCM.n,0);     % switch on bilinear modulations
DCM.c       = zeros(DCM.n,0);           % switch on exogenous connections
DCM.d       = zeros(DCM.n,DCM.n,0);     % switch on exogenous connections
DCM.Y.y     = y;                        % responses over time
DCM.Y.dt    = dt;                       % repetition time

DCM.options.maxnodes    = DCM.n;    % maximum number of nodes
DCM.options.precision   = log(64);  % log precision on intraregional priors
DCM.options.two_state   = 0;        % one or two states per region
DCM.options.stochastic  = 0;        % exogenous or endogenous fluctuations
DCM.options.centre      = 0;        % mean-centre inputs
DCM.options.analysis    = 'CSD';    % type of analysis
DCM.options.order       = 8;        % precision spectral observation noise
DCM.options.nograph     = false;    % graphical display
DCM.options.maxit       = 128;      % maximum number of iterations
DCM.options.induced     = 1;        % switch for CSD data features
DCM.options.nonlinear   = 0;        % interactions among hidden states
% -------------------------------------------------------------------------

% Update covariance of interregional priors
[DCM.M.pE, DCM.M.pC, ~, ~] = spm_dcm_fmri_priors(DCM.a,DCM.b,DCM.c,...
    DCM.d,DCM.options);
indx = find(DCM.a==1 & ~eye(DCM.n));
DCM.M.pC(sub2ind(size(DCM.M.pC), indx', indx')) = v;

% DCM inversion
P = spm_dcm_fit(DCM);

end
