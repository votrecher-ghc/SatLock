function [result_row, out_dir] = run_simulation_workflow_best(obs_filepath, nav_filepath, target_letter, user_cfg)
% RUN_SIMULATION_WORKFLOW_BEST
% One-click simulation workflow:
%   parse -> inject ideal gesture -> preprocess -> recognize -> evaluate.
%
% Example:
%   [row, out_dir] = run_simulation_workflow_best( ...
%       'data/1_8/A_1_8_1.obs', 'data/1_8/2026_1_8.nav', 'A');

if nargin < 1 || isempty(obs_filepath)
    obs_filepath = fullfile('data', '1_8', 'A_1_8_1.obs');
end
if nargin < 2 || isempty(nav_filepath)
    nav_filepath = fullfile('data', '1_8', '2026_1_8.nav');
end
if nargin < 3 || isempty(target_letter)
    target_letter = 'A';
end
if nargin < 4
    user_cfg = struct();
end

t = benchmark_sim_oracle_vs_baselines(obs_filepath, nav_filepath, {target_letter}, user_cfg);
result_row = t(1, :);
pp = result_row.plot_path;
if iscell(pp)
    pp = pp{1};
elseif isstring(pp)
    pp = char(pp(1));
elseif ischar(pp)
    % keep as is
else
    pp = char(string(pp));
end
out_dir = fileparts(pp);

fprintf('\n== Simulation Workflow Result ==\n');
fprintf('Letter        : %s\n', string(result_row.letter));
fprintf('Baseline RMSE : %.4f m\n', result_row.baseline_rmse_m);
fprintf('Known RMSE    : %.4f m (%.2f%%)\n', result_row.known_rmse_m, result_row.known_err_pct);
fprintf('Known pass20  : %d\n', result_row.known_pass20);
fprintf('Auto RMSE     : %.4f m (%.2f%%)\n', result_row.auto_rmse_m, result_row.auto_err_pct);
fprintf('Auto pass20   : %d\n', result_row.auto_pass20);
fprintf('Auto letter   : %s\n', string(result_row.auto_best_letter));
fprintf('Output dir    : %s\n', out_dir);
end
