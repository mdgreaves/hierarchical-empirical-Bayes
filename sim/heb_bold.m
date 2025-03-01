function y = heb_bold(T, TR, A_mat, snr, tau)
% heb_bold - Simulates BOLD time-series data with specified SNR for 
%            observation noise.
%
% Syntax:  y = heb_bold(T, TR, A_mat, snr, tau)
%
% Inputs:
%   T      - Number of time points (observations) for the simulation.
%   TR     - Repetition time (seconds), defining the temporal resolution.
%   A_mat  - Effective connectivity matrix (n x n).
%   snr    - Desired signal-to-noise ratio.
%   tau    - Transit times.
%
% Outputs:
%   y      - Simulated BOLD time-series data (T x n), where T is the number
%            of time points and n is the number of regions.
%
% Description:
%   This function simulates neuronal and haemodynamic dynamics based on an
%   effective connectivity matrix `A_mat`. The neuronal activity is
%   integrated using the SPM Jacobian-based integrator `spm_int_J`, and the
%   haemodynamic response is observed via `spm_gx_fmri`. The scheme is
%   based on the demonstration in SPM's `DEM_demo_connectivity_fMRI`.
%
%   Key Features:
%   - Uses priors from `spm_dcm_fmri_priors` to specify system properties.
%   - Integrates neural activity with `spm_int_J`.
%   - Observes haemodynamic responses with `spm_gx_fmri`.
%   - Adds biologically realistic observation noise.
%
% Example:
%   A_mat = randn(6, 6) * 0.2;  % Example connectivity matrix
%   T = 1000;                   % 1000 time points
%   TR = 2;                     % Repetition time of 2 seconds
%   y = heb_bold(T, TR, A_mat); % Simulated BOLD time series
%
% See also:
%   spm_int_J, spm_gx_fmri, spm_dcm_fmri_priors, DEM_demo_connectivity_fMRI

% Options and priors
n                   = size(A_mat,1);
options.maxnodes    = n;
options.nonlinear   = 0;
options.two_state   = 0;
options.stochastic  = 0;
options.centre      = 1;
options.induced     = 1;
A                   = ones(n,n);
B                   = zeros(n,n,0);
C                   = zeros(n,n);
D                   = zeros(n,n,0);
[pP, ~, M.x, ~]     = spm_dcm_fmri_priors(A,B,C,D,options);

% True parameters
pP.A = A_mat;
pP.C = eye(n,n);
pP.transit = tau;

% Integrate states
U.u                 = spm_rand_mar(T,n,1/2)/4;
U.dt                = TR;
M.f                 = 'spm_fx_fmri';
x                   = spm_int_J(pP,M,U);

% Blood-oxygen-level-dependent (BOLD) signals
y                   = zeros(T,n);
for i = 1:T
    y(i,:) = spm_gx_fmri(spm_unvec(x(i,:),M.x),[],pP)';
end

% Calculate region-specific noise variances based on SNR
signal_variance     = var(y, 0, 1);             % Variance of BOLD signal
noise_variance      = signal_variance ./ snr;   % Desired noise variance

% Generate region-specific MAR observation noise
e                   = zeros(T, n);
for region = 1:n
    e_raw           = spm_rand_mar(T, 1, 1/2);
    e(:, region)    = e_raw * sqrt(noise_variance(region) / var(e_raw)); 
end

% Add region-specific noise to BOLD signal
y                   = y + e;

end