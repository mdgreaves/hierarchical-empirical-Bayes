function out = heb_perm_pipeline_null_pseudo(P, C, network, varargin)
% =========================================================================
% heb_perm_pipeline_null_pseudo.m
% =========================================================================
% Template code for a full-pipeline permutation null under scrambled
% structural connectivity.
%
% Context:
% The full analysis reported in the Supporting Information was run on an
% HPC system because it required 1000 full-pipeline permutations for each
% of 17 networks. The raw DCM/HEB structures are too large to distribute
% with this repository, so this file provides a concrete template showing
% how the analysis was performed with user-supplied first-level DCMs and a
% matched structural-connectivity matrix.
%
% "Pseudo" here means that data loading is intentionally left to the user.
% The model-scoring and re-inversion steps are concrete MATLAB code.
%
% Pipeline:
%   1) Invert the corresponding uninformed hierarchical model.
%   2) Score empirical structure-based priors over a grid of alpha,beta
%      mappings using Bayesian model reduction.
%   3) Bayesian-model-average the grid to obtain an evidence-weighted
%      prior-variance transformation.
%   4) Re-invert the hierarchical model de novo under that informed prior.
%   5) Repeat steps 2-4 after scrambling off-diagonal structural-
%      connectivity values while preserving symmetry.
%   6) Compare each permuted model's free energy to the uninformed model.
%
% Usage:
%   % P is a cell array of inverted first-level DCMs.
%   % C is the normalized [0,1] structural-connectivity matrix.
%   out = heb_perm_pipeline_null_pseudo(P, C, 'ContA', ...
%       'B', 1000, ...
%       'OutputDir', fullfile(pwd, 'pipeline_null_ContA'));
%
% HPC usage:
%   % For one permutation per Slurm array task:
%   b = str2double(getenv('SLURM_ARRAY_TASK_ID'));
%   out = heb_perm_pipeline_null_pseudo(P, C, 'ContA', ...
%       'B', 1000, ...
%       'PermutationIndices', b, ...
%       'OutputDir', fullfile(pwd, 'pipeline_null_ContA'));
%
% Outputs:
%   - pipeline_null_baseline_<network>.mat
%   - null_pipe_BFs_<network>_itt_<b>.mat for each permutation b
%
% Dependencies:
% SPM12 must already be on the MATLAB path. This script also uses core
% functions in ../core and adds that directory to the MATLAB path.
%
% =========================================================================

repo_root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(repo_root, 'core'));

check_dependencies();
validate_inputs(P, C);

opts = parse_options(varargin{:});
if isnan(opts.SampleSize)
    opts.SampleSize = numel(P);
end
network = char(network);
if isempty(opts.PermutationIndices)
    opts.PermutationIndices = permutation_indices_from_environment(opts.B);
end
if ~exist(opts.OutputDir, 'dir')
    mkdir(opts.OutputDir);
end

params = default_grid_params(P, network, opts);
baseline = load_or_compute_baseline(P, C, network, params, opts);

result_template = empty_result();
results = repmat(result_template, numel(opts.PermutationIndices), 1);

for ii = 1:numel(opts.PermutationIndices)
    b = opts.PermutationIndices(ii);
    out_file = fullfile(opts.OutputDir, ...
        sprintf('null_pipe_BFs_%s_itt_%04d.mat', network, b));

    if opts.ReuseExisting && exist(out_file, 'file')
        loaded = load(out_file, 'result');
        results(ii) = loaded.result;
        fprintf('Using existing result: %s\n', out_file);
        continue
    end

    fprintf('\nPermutation %d/%d for %s\n', b, opts.B, network);
    C_perm = scramble_symmetric_connectivity(C, opts.Seed + b);
    winning = score_mapping_grid(baseline.PEB_uninformed, C_perm, params);
    M_perm = informed_model_structure(baseline.PEB_uninformed.M, ...
        C_perm, winning.alpha_bma, winning.beta_bma);

    [PEB_perm, ~] = heb_capture(P, M_perm, params.field);

    result = result_template;
    result.net = network;
    result.itt = b;
    result.B = opts.B;
    result.sample_size = opts.SampleSize;
    result.F_uninformed_existing = baseline.F_uninformed_existing;
    result.F_empirical_existing = baseline.F_empirical_existing;
    result.perm_PEB_F = PEB_perm.F;
    result.perm_logBF = PEB_perm.F - baseline.F_uninformed_existing;
    result.perm_scaled_logBF = result.perm_logBF / opts.SampleSize;
    result.empirical_logBF = baseline.F_empirical_existing - ...
        baseline.F_uninformed_existing;
    result.empirical_scaled_logBF = result.empirical_logBF / ...
        opts.SampleSize;
    result.winning = winning;
    result.empirical_winning = baseline.empirical_winning;
    result.SC_empirical = C;
    result.SC_perm = C_perm;
    result.perm_prior_variance_modulation_index = ...
        prior_variance_modulation_index(C_perm, ...
            winning.alpha_bma, winning.beta_bma);
    result.empirical_prior_variance_modulation_index = ...
        baseline.empirical_prior_variance_modulation_index;

    save(out_file, 'result', '-v7.3');
    fprintf('Saved %s\n', out_file);
    results(ii) = result;
end

out = struct();
out.network = network;
out.options = opts;
out.params = params;
out.baseline = baseline;
out.results = results;

end

function opts = parse_options(varargin)
ip = inputParser;
ip.addParameter('B', 1000, @(x) isnumeric(x) && isscalar(x) && x >= 1);
ip.addParameter('PermutationIndices', [], @(x) isempty(x) || isnumeric(x));
ip.addParameter('OutputDir', fullfile(pwd, 'pipeline_null'), ...
    @(x) ischar(x) || isstring(x));
ip.addParameter('SampleSize', [], @(x) isempty(x) || ...
    (isnumeric(x) && isscalar(x) && x >= 1));
ip.addParameter('Seed', 5489, @(x) isnumeric(x) && isscalar(x));
ip.addParameter('Field', 'A', @(x) ischar(x) || isstring(x));
ip.addParameter('NumGridPoints', 30, @(x) isnumeric(x) && ...
    isscalar(x) && x >= 2);
ip.addParameter('Alphas', [], @(x) isempty(x) || isnumeric(x));
ip.addParameter('Betas', [], @(x) isempty(x) || isnumeric(x));
ip.addParameter('ReuseExisting', true, @(x) ...
    (islogical(x) || isnumeric(x)) && isscalar(x));
ip.parse(varargin{:});
opts = ip.Results;
opts.OutputDir = char(opts.OutputDir);
opts.Field = char(opts.Field);
opts.B = round(opts.B);
opts.NumGridPoints = round(opts.NumGridPoints);
opts.ReuseExisting = logical(opts.ReuseExisting);
if isempty(opts.SampleSize)
    opts.SampleSize = nan;
end
opts.PermutationIndices = opts.PermutationIndices(:)';
end

function idx = permutation_indices_from_environment(B)
task_id = str2double(getenv('SLURM_ARRAY_TASK_ID'));
if isfinite(task_id) && task_id >= 1
    idx = task_id;
else
    idx = 1:B;
end
end

function baseline = load_or_compute_baseline(P, C, network, params, opts)
baseline_file = fullfile(opts.OutputDir, ...
    sprintf('pipeline_null_baseline_%s.mat', network));

if opts.ReuseExisting && exist(baseline_file, 'file')
    loaded = load(baseline_file, 'baseline');
    baseline = loaded.baseline;
    fprintf('Using existing baseline: %s\n', baseline_file);
    return
end

fprintf('\nComputing uninformed and empirical baselines for %s\n', network);
[PEB_uninformed, ~] = heb_capture(P, struct(), params.field);

empirical_winning = score_mapping_grid(PEB_uninformed, C, params);
M_empirical = informed_model_structure(PEB_uninformed.M, C, ...
    empirical_winning.alpha_bma, empirical_winning.beta_bma);
[PEB_empirical, ~] = heb_capture(P, M_empirical, params.field);

baseline = struct();
baseline.net = network;
baseline.PEB_uninformed = PEB_uninformed;
baseline.F_uninformed_existing = PEB_uninformed.F;
baseline.F_empirical_existing = PEB_empirical.F;
baseline.empirical_winning = empirical_winning;
baseline.empirical_prior_variance_modulation_index = ...
    prior_variance_modulation_index(C, empirical_winning.alpha_bma, ...
        empirical_winning.beta_bma);
baseline.sample_size = opts.SampleSize;

save(baseline_file, 'baseline', '-v7.3');
fprintf('Saved baseline: %s\n', baseline_file);
end

function params = default_grid_params(P, network, opts)
params = struct();
params.field = opts.Field;
params.name = network;

if ~isempty(opts.Alphas) && ~isempty(opts.Betas)
    params.alphas = opts.Alphas(:);
    params.betas = opts.Betas(:);
    return
end

n = P{1}.n;
diag_idx = sub2ind(size(P{1}.M.pC), (1:n^2)', (1:n^2)');
state_prior_variances = full(P{1}.M.pC(diag_idx));
state_prior_variances = reshape(state_prior_variances, n, n);
max_offdiag_variance = max(state_prior_variances(~eye(n)));

alpha_values = linspace(-max_offdiag_variance, ...
    max_offdiag_variance, opts.NumGridPoints);

alphas = [];
betas = [];
for i = 1:numel(alpha_values)
    beta_values = linspace(0, max_offdiag_variance - alpha_values(i), ...
        opts.NumGridPoints);
    alphas = [alphas; repmat(alpha_values(i), numel(beta_values), 1)]; %#ok<AGROW>
    betas = [betas; beta_values(:)]; %#ok<AGROW>
end

epsilon = 1e-5;
valid = alphas >= epsilon & ...
    alphas <= max_offdiag_variance & ...
    alphas + betas >= epsilon & ...
    alphas + betas <= max_offdiag_variance;

params.alphas = alphas(valid);
params.betas = betas(valid);
end

function winning = score_mapping_grid(PEB_uninformed, C, params)
n_models = numel(params.alphas);
Fs = nan(n_models, 1);

fprintf('Scoring %d alpha,beta mappings via BMR', n_models);
for i = 1:n_models
    if mod(i, max(1, round(n_models / 20))) == 0
        fprintf('.');
    end
    Fs(i) = heb_bmr(params.alphas(i), params.betas(i), C, ...
        PEB_uninformed);
end
fprintf('\n');

[maxF, max_idx] = max(Fs);
if any(~isfinite(Fs))
    error('Grid scoring returned non-finite free-energy values.');
end

weights = exp(Fs - maxF);
weights = weights ./ sum(weights);

winning = struct();
winning.maxF = maxF;
winning.alpha = params.alphas(max_idx);
winning.beta = params.betas(max_idx);
winning.alpha_bma = sum(weights .* params.alphas);
winning.beta_bma = sum(weights .* params.betas);
winning.alpha_weighted_mean = winning.alpha_bma;
winning.beta_weighted_mean = winning.beta_bma;
winning.Fs = Fs;
winning.weights = weights;
end

function M = informed_model_structure(M0, C, alpha, beta)
n = size(C, 1);
variance = beta .* C + alpha;
offdiag = find(~eye(n));

rC = M0.pC;
rC(sub2ind(size(rC), offdiag(:), offdiag(:))) = variance(offdiag);

M = M0;
M.bC = rC;
M.bE = M0.pE;
end

function C_perm = scramble_symmetric_connectivity(C, seed)
rng(seed, 'twister');
n = size(C, 1);
lower_idx = tril(true(n), -1);
values = C(lower_idx);
values = values(randperm(numel(values)));

C_perm = zeros(n);
C_perm(lower_idx) = values;
C_perm = C_perm + C_perm';
C_perm(1:n+1:end) = diag(C);
end

function pvmi = prior_variance_modulation_index(C, alpha, beta)
variance = beta .* C + alpha;
offdiag = variance(~eye(size(C, 1)));
mu = mean(offdiag(:));
if isempty(offdiag) || ~isfinite(mu) || abs(mu) < eps
    pvmi = nan;
else
    pvmi = std(offdiag(:)) / mu;
end
end

function result = empty_result()
result = struct( ...
    'net', '', ...
    'itt', nan, ...
    'B', nan, ...
    'sample_size', nan, ...
    'F_uninformed_existing', nan, ...
    'F_empirical_existing', nan, ...
    'perm_PEB_F', nan, ...
    'perm_logBF', nan, ...
    'perm_scaled_logBF', nan, ...
    'empirical_logBF', nan, ...
    'empirical_scaled_logBF', nan, ...
    'winning', struct(), ...
    'empirical_winning', struct(), ...
    'SC_empirical', [], ...
    'SC_perm', [], ...
    'perm_prior_variance_modulation_index', nan, ...
    'empirical_prior_variance_modulation_index', nan);
end

function validate_inputs(P, C)
assert(iscell(P) && ~isempty(P), 'P must be a non-empty cell array.');
assert(isfield(P{1}, 'n'), 'Each DCM must contain field .n.');
n = P{1}.n;
assert(isequal(size(C), [n n]), ...
    'C must be an n-by-n matrix matching the first-level DCMs.');
assert(all(isfinite(C(:))), 'C must contain finite values.');
assert(all(C(:) >= 0 & C(:) <= 1), ...
    'C should be normalized to the interval [0,1].');
end

function check_dependencies()
required = {'heb_capture', 'heb_bmr', 'spm_dcm_peb', ...
    'spm_log_evidence_reduce'};
missing = required(cellfun(@(f) isempty(which(f)), required));
if ~isempty(missing)
    error('Missing required function(s): %s', strjoin(missing, ', '));
end
end
