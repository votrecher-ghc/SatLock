function summary_tbl = benchmark_sim_oracle_vs_baselines(obs_filepath, nav_filepath, letters, user_cfg)
% BENCHMARK_SIM_ORACLE_VS_BASELINES
% Compare baseline tracker and simulation-oracle tracker on ideal synthetic data.
%
% Default:
%   summary_tbl = benchmark_sim_oracle_vs_baselines();
%
% Notes:
%   - "known-letter" mode uses the generated target label.
%   - "auto-letter" mode lets run_gesture_analysis_sim_oracle infer the label.

if nargin < 1 || isempty(obs_filepath)
    obs_filepath = fullfile('data', '1_8', 'A_1_8_1.obs');
end
if nargin < 2 || isempty(nav_filepath)
    nav_filepath = fullfile('data', '1_8', '2026_1_8.nav');
end
if nargin < 3 || isempty(letters)
    letters = {'A', 'B', 'M', 'Star', 'L', 'X', 'Z', 'N'};
end
if nargin < 4 || isempty(user_cfg)
    user_cfg = struct();
end
if ischar(letters) || isstring(letters)
    letters = cellstr(letters);
end

cfg = default_cfg();
cfg = merge_cfg(cfg, user_cfg);

this_dir = fileparts(mfilename('fullpath'));
repo_dir = fileparts(this_dir);
addpath(genpath(repo_dir));

obs_filepath = resolve_data_path(obs_filepath, repo_dir);
nav_filepath = resolve_data_path(nav_filepath, repo_dir);

fprintf('\n== Sim-Oracle Benchmark ==\n');
fprintf('obs: %s\n', obs_filepath);
fprintf('nav: %s\n', nav_filepath);
fprintf('letters: %s\n\n', strjoin(upper(letters), ', '));

out_root = fullfile(this_dir, 'results');
if ~exist(out_root, 'dir')
    mkdir(out_root);
end
stamp = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
out_dir = fullfile(out_root, ['sim_oracle_', stamp]);
mkdir(out_dir);

if ~cfg.figure_visible
    old_fig_vis = get(0, 'DefaultFigureVisible');
    cleanup_obj = onCleanup(@() set(0, 'DefaultFigureVisible', old_fig_vis)); %#ok<NASGU>
    set(0, 'DefaultFigureVisible', 'off');
end

obs_base = parse_rinex_obs(obs_filepath);
nav_data = parse_rinex_nav_multi_gnss(nav_filepath);

n_case = numel(letters);
results = repmat(struct( ...
    'letter', '', ...
    'baseline_rmse_m', NaN, ...
    'baseline_start_err_m', NaN, ...
    'baseline_end_err_m', NaN, ...
    'baseline_coverage', NaN, ...
    'known_rmse_m', NaN, ...
    'known_start_err_m', NaN, ...
    'known_end_err_m', NaN, ...
    'known_coverage', NaN, ...
    'known_err_pct', NaN, ...
    'known_pass20', false, ...
    'auto_rmse_m', NaN, ...
    'auto_start_err_m', NaN, ...
    'auto_end_err_m', NaN, ...
    'auto_coverage', NaN, ...
    'auto_err_pct', NaN, ...
    'auto_pass20', false, ...
    'auto_best_letter', '', ...
    'plot_path', ''), n_case, 1);

for i = 1:n_case
    letter_raw = strtrim(char(letters{i}));
    letter = upper(letter_raw);
    sim_letter = letter;
    if strcmp(letter, 'STAR')
        sim_letter = 'Star';
    end

    fprintf('\n[Case %d/%d] Letter %s\n', i, n_case, letter);
    close all;
    rng(1000 + i, 'twister');

    obs_sim = generate_ideal_multi_shape(obs_base, nav_data, sim_letter);
    [obs_clean, step1_res] = gesture_analysis_baseline_gvi(obs_sim);
    [obs_waveform, step1_res_shaped] = waveform_drop_recovery_reshaping(obs_clean, step1_res);
    if cfg.run_baseline
        [obs_aligned, step1_aligned] = gesture_analysis_align_duration(obs_waveform, step1_res_shaped);
    else
        obs_aligned = [];
        step1_aligned = struct();
    end

    if cfg.run_baseline
        [base_x, base_y, base_t, ~] = run_gesture_analysis_continuous_track_line(obs_aligned, nav_data, step1_aligned);
    else
        base_x = []; base_y = []; base_t = [];
    end

    cfg_known = cfg;
    cfg_known.sim.target_letter = sim_letter;
    cfg_known.sim.auto_detect_letter = false;
    [k_x, k_y, k_t, k_conf, ~] = run_gesture_analysis_sim_oracle(obs_waveform, nav_data, step1_res_shaped, cfg_known);

    if cfg.run_auto_mode
        cfg_auto = cfg;
        cfg_auto.sim.target_letter = '';
        cfg_auto.sim.auto_detect_letter = true;
        [a_x, a_y, a_t, a_conf, dbg_auto] = run_gesture_analysis_sim_oracle(obs_waveform, nav_data, step1_res_shaped, cfg_auto);
        auto_best_letter = string(dbg_auto.best_letter);
    else
        a_x = []; a_y = []; a_t = []; a_conf = [];
        auto_best_letter = "";
    end

    N = numel(step1_res.t_grid);
    [gt_x, gt_y, gt_pen_down] = build_ground_truth(letter, N);

    [base_fx, base_fy] = to_full_series(base_x, base_y, base_t, N);
    [k_fx, k_fy] = to_full_series(k_x, k_y, k_t, N);
    [k_cf, ~] = to_full_series(k_conf, k_conf, k_t, N);
    [a_fx, a_fy] = to_full_series(a_x, a_y, a_t, N);
    [a_cf, ~] = to_full_series(a_conf, a_conf, a_t, N);

    if cfg.run_baseline
        base_met = evaluate_against_gt(base_fx, base_fy, gt_x, gt_y, gt_pen_down, cfg.max_shift_eval);
    else
        base_met = make_nan_metrics();
    end
    k_met = evaluate_against_gt(k_fx, k_fy, gt_x, gt_y, gt_pen_down, cfg.max_shift_eval);
    if cfg.run_auto_mode
        a_met = evaluate_against_gt(a_fx, a_fy, gt_x, gt_y, gt_pen_down, cfg.max_shift_eval);
    else
        a_met = make_nan_metrics();
    end

    k_pct = estimate_error_percent(letter, k_met.rmse_m);
    if cfg.run_auto_mode
        a_pct = estimate_error_percent(letter, a_met.rmse_m);
    else
        a_pct = NaN;
    end

    fig_path = fullfile(out_dir, sprintf('compare_%s.png', letter));
    make_compare_plot(letter, gt_x, gt_y, gt_pen_down, ...
        base_fx, base_fy, base_met, ...
        k_fx, k_fy, k_cf, k_met, ...
        a_fx, a_fy, a_cf, a_met, ...
        cfg.run_auto_mode, cfg.figure_visible, fig_path);

    if cfg.run_baseline
        fprintf('    Baseline RMSE=%.4f | Known RMSE=%.4f (%.2f%%)\n', ...
            base_met.rmse_m, k_met.rmse_m, k_pct);
    else
        fprintf('    Known RMSE=%.4f (%.2f%%)\n', k_met.rmse_m, k_pct);
    end
    if cfg.run_auto_mode
        fprintf('    Auto RMSE=%.4f (%.2f%%) | best=%s\n', ...
            a_met.rmse_m, a_pct, auto_best_letter);
    end

    results(i).letter = letter;
    results(i).baseline_rmse_m = base_met.rmse_m;
    results(i).baseline_start_err_m = base_met.start_err_m;
    results(i).baseline_end_err_m = base_met.end_err_m;
    results(i).baseline_coverage = base_met.coverage;
    results(i).known_rmse_m = k_met.rmse_m;
    results(i).known_start_err_m = k_met.start_err_m;
    results(i).known_end_err_m = k_met.end_err_m;
    results(i).known_coverage = k_met.coverage;
    results(i).known_err_pct = k_pct;
    results(i).known_pass20 = (k_pct <= 20);
    results(i).auto_rmse_m = a_met.rmse_m;
    results(i).auto_start_err_m = a_met.start_err_m;
    results(i).auto_end_err_m = a_met.end_err_m;
    results(i).auto_coverage = a_met.coverage;
    results(i).auto_err_pct = a_pct;
    results(i).auto_pass20 = (cfg.run_auto_mode && a_pct <= 20);
    results(i).auto_best_letter = auto_best_letter;
    results(i).plot_path = fig_path;
end

summary_tbl = struct2table(results);
csv_path = fullfile(out_dir, 'summary.csv');
writetable(summary_tbl, csv_path);
save(fullfile(out_dir, 'summary.mat'), 'summary_tbl', 'results', 'out_dir');

fprintf('\nBenchmark complete. Output: %s\n', out_dir);
disp(summary_tbl(:, {'letter', 'known_rmse_m', 'known_err_pct', 'known_pass20', 'auto_rmse_m', 'auto_err_pct', 'auto_pass20', 'auto_best_letter'}));
end

% -------------------------------------------------------------------------
function cfg = default_cfg()
cfg = struct();
cfg.run_auto_mode = false;
cfg.run_baseline = false;
cfg.figure_visible = false;
cfg.max_shift_eval = 75;
cfg.sim = struct();
end

% -------------------------------------------------------------------------
function met = make_nan_metrics()
met = struct( ...
    'rmse_m', NaN, ...
    'best_shift_samples', NaN, ...
    'start_err_m', NaN, ...
    'end_err_m', NaN, ...
    'coverage', NaN);
end

% -------------------------------------------------------------------------
function p = resolve_data_path(p_in, repo_dir)
p = char(p_in);
if exist(p, 'file')
    return;
end
cands = {
    fullfile(repo_dir, p), ...
    fullfile(repo_dir, 'data', '1_8', p)
};
for i = 1:numel(cands)
    if exist(cands{i}, 'file')
        p = cands{i};
        return;
    end
end
error('File not found: %s', p_in);
end

% -------------------------------------------------------------------------
function [fx, fy] = to_full_series(x, y, t_idx, N)
fx = nan(N, 1);
fy = nan(N, 1);
if isempty(x) || isempty(y) || isempty(t_idx)
    return;
end
x = x(:); y = y(:); t_idx = round(t_idx(:));
n = min([numel(x), numel(y), numel(t_idx)]);
x = x(1:n); y = y(1:n); t_idx = t_idx(1:n);
keep = isfinite(x) & isfinite(y) & isfinite(t_idx) & t_idx >= 1 & t_idx <= N;
if ~any(keep)
    return;
end
x = x(keep); y = y(keep); t_idx = t_idx(keep);
[u_idx, ia] = unique(t_idx, 'stable');
fx(u_idx) = x(ia);
fy(u_idx) = y(ia);
end

% -------------------------------------------------------------------------
function met = evaluate_against_gt(est_x, est_y, gt_x, gt_y, gt_mask, max_shift)
met = struct('rmse_m', inf, 'best_shift_samples', 0, 'start_err_m', inf, 'end_err_m', inf, 'coverage', 0);
N = numel(gt_x);
idx_gt = find(gt_mask & isfinite(gt_x) & isfinite(gt_y));
if isempty(idx_gt)
    met.rmse_m = NaN; met.start_err_m = NaN; met.end_err_m = NaN; met.coverage = NaN;
    return;
end
valid_est = find(isfinite(est_x) & isfinite(est_y));
if isempty(valid_est)
    return;
end
met.coverage = nnz(isfinite(est_x(idx_gt)) & isfinite(est_y(idx_gt))) / numel(idx_gt);
if numel(valid_est) >= 2
    ex = interp1(valid_est, est_x(valid_est), 1:N, 'pchip', 'extrap')';
    ey = interp1(valid_est, est_y(valid_est), 1:N, 'pchip', 'extrap')';
else
    ex = repmat(est_x(valid_est(1)), N, 1);
    ey = repmat(est_y(valid_est(1)), N, 1);
end

best_rmse = inf;
best_shift = 0;
for sh = -max_shift:max_shift
    idx_est = idx_gt + sh;
    keep = idx_est >= 1 & idx_est <= N;
    if nnz(keep) < max(8, round(0.2 * numel(idx_gt)))
        continue;
    end
    g = idx_gt(keep); e = idx_est(keep);
    rmse = sqrt(mean((ex(e) - gt_x(g)).^2 + (ey(e) - gt_y(g)).^2));
    if rmse < best_rmse
        best_rmse = rmse;
        best_shift = sh;
    end
end

met.rmse_m = best_rmse;
met.best_shift_samples = best_shift;

sg = idx_gt(1); eg = idx_gt(end);
se = sg + best_shift; ee = eg + best_shift;
if se >= 1 && se <= N
    met.start_err_m = hypot(ex(se) - gt_x(sg), ey(se) - gt_y(sg));
else
    met.start_err_m = NaN;
end
if ee >= 1 && ee <= N
    met.end_err_m = hypot(ex(ee) - gt_x(eg), ey(ee) - gt_y(eg));
else
    met.end_err_m = NaN;
end
end

% -------------------------------------------------------------------------
function pct = estimate_error_percent(letter, rmse_m)
[dx, dy] = letter_span(letter);
diag_len = hypot(dx, dy);
pct = 100 * rmse_m / max(diag_len, eps);
end

% -------------------------------------------------------------------------
function [dx, dy] = letter_span(letter)
switch upper(string(letter))
    case "A",    dx = 0.80; dy = 1.00;
    case "L",    dx = 0.60; dy = 0.80;
    case "Z",    dx = 1.00; dy = 0.80;
    case "N",    dx = 0.60; dy = 0.80;
    case "STAR", dx = 0.96; dy = 1.00;
    case "B",    dx = 1.90; dy = 2.00;
    case "M",    dx = 0.80; dy = 0.80;
    case "X",    dx = 0.60; dy = 0.80;
    otherwise,   dx = 1.00; dy = 1.00;
end
end

% -------------------------------------------------------------------------
function [gt_x, gt_y, gt_pen_down] = build_ground_truth(letter, num_samples)
sampling_rate = 25;
stages = letter_stages(letter);
gt_x = nan(num_samples, 1);
gt_y = nan(num_samples, 1);
gt_pen_down = false(num_samples, 1);

total_dur = 0;
for k = 1:size(stages, 1)
    total_dur = total_dur + stages{k, 3};
end
start_idx = round(num_samples * 0.3);
end_idx = min(num_samples, start_idx + round(total_dur * sampling_rate));

for t_idx = start_idx:end_idx
    dt = (t_idx - start_idx) / sampling_rate;
    elapsed = 0;
    for k = 1:size(stages, 1)
        p1 = stages{k, 1}; p2 = stages{k, 2}; dur = stages{k, 3}; pen = stages{k, 4};
        if dt <= (elapsed + dur)
            alpha = max(0, min(1, (dt - elapsed) / max(dur, eps)));
            pos = p1 + (p2 - p1) * alpha;
            gt_x(t_idx) = pos(1);
            gt_y(t_idx) = pos(2);
            gt_pen_down(t_idx) = logical(pen);
            break;
        end
        elapsed = elapsed + dur;
    end
end
end

% -------------------------------------------------------------------------
function stages = letter_stages(letter)
switch upper(string(letter))
    case "A"
        P1 = [-0.40, -0.50]; P2 = [0.00, 0.50]; P3 = [0.40, -0.50];
        P4 = [-0.20, -0.10]; P5 = [0.20, -0.10];
        stages = {P1,P2,1.5,true; P2,P2,3.0,true; P2,P3,1.5,true; P3,P4,3.0,false; P4,P5,1.5,true};
    case "B"
        P1 = [-0.40,-1]; P2 = [-0.40,1]; P3 = [1.5,0.40];
        P4 = [-0.40,0]; P5 = [1.5,0]; P6 = [-0.40,-1];
        stages = {P1,P2,1.5,true; P2,P2,3.0,true; P2,P3,1.5,true; P3,P3,3.0,true; ...
                  P3,P4,3.0,true; P4,P4,3.0,true; P4,P5,1.5,true; P5,P5,3.0,true; P5,P6,1.5,true};
    case "M"
        P1=[-0.40,-0.40]; P2=[-0.40,0.40]; P3=[0.00,0.00]; P4=[0.40,0.40]; P5=[0.40,-0.40];
        stages = {P1,P2,1.5,true; P2,P2,3.0,true; P2,P3,1.5,true; P3,P3,3.0,true; ...
                  P3,P4,1.5,true; P4,P4,3.0,true; P4,P5,1.5,true};
    case "STAR"
        P1=[-0.30,-0.45]; P2=[0.00,0.55]; P3=[0.30,-0.45]; P4=[-0.48,0.15]; P5=[0.48,0.15];
        stages = {P1,P2,1.5,true; P2,P2,3.0,true; P2,P3,1.5,true; P3,P3,3.0,true; ...
                  P3,P4,1.5,true; P4,P4,3.0,true; P4,P5,1.5,true; P5,P5,3.0,true; P5,P1,1.5,true};
    case "L"
        P1=[-0.3,0.4]; P2=[-0.3,-0.4]; P3=[0.3,-0.4];
        stages = {P1,P2,1.5,true; P2,P2,3.0,true; P2,P3,1.5,true};
    case "X"
        P1=[-0.3,0.4]; P2=[0.3,-0.4]; P3=[0.3,0.4]; P4=[-0.3,-0.4];
        stages = {P1,P2,1.5,true; P2,P3,3.0,false; P3,P4,1.5,true};
    case "Z"
        P1=[-0.3,0.4]; P2=[0.3,0.4]; P3=[-0.3,-0.4]; P4=[0.7,-0.4];
        stages = {P1,P2,1.5,true; P2,P2,3.0,true; P2,P3,1.5,true; P3,P3,3.0,true; P3,P4,1.5,true};
    case "N"
        P1=[-0.3,-0.4]; P2=[-0.3,0.4]; P3=[0.3,-0.4]; P4=[0.3,0.4];
        stages = {P1,P2,1.5,true; P2,P2,3.0,true; P2,P3,1.5,true; P3,P3,3.0,true; P3,P4,1.5,true};
    otherwise
        error('Unsupported letter: %s', char(letter));
end
end

% -------------------------------------------------------------------------
function make_compare_plot(letter, gt_x, gt_y, gt_mask, ...
    base_x, base_y, base_met, ...
    k_x, k_y, k_conf, k_met, ...
    a_x, a_y, a_conf, a_met, ...
    run_auto, figure_visible, out_path)

fig_vis = 'off';
if figure_visible
    fig_vis = 'on';
end
f = figure('Color', 'w', 'Position', [100, 100, 860, 700], 'Visible', fig_vis);
tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
hold on; grid on; axis equal;
idx_pen = gt_mask & isfinite(gt_x) & isfinite(gt_y);
idx_up = ~gt_mask & isfinite(gt_x) & isfinite(gt_y);

plot(gt_x(idx_up), gt_y(idx_up), '--', 'Color', [0.75 0.75 0.75], 'LineWidth', 1.0, 'DisplayName', 'GT Pen Up');
plot(gt_x(idx_pen), gt_y(idx_pen), '-', 'Color', [0.10 0.45 0.95], 'LineWidth', 3.0, 'DisplayName', 'GT Pen Down');
if any(isfinite(base_x)) && any(isfinite(base_y))
    plot(base_x, base_y, ':', 'Color', [0.95 0.55 0.10], 'LineWidth', 2.0, 'DisplayName', sprintf('Baseline (RMSE %.3f)', base_met.rmse_m));
end
plot(k_x, k_y, '-', 'Color', [0.85 0.20 0.20], 'LineWidth', 2.0, 'DisplayName', sprintf('Oracle Known (RMSE %.3f)', k_met.rmse_m));
if run_auto
    plot(a_x, a_y, '-', 'Color', [0.10 0.60 0.20], 'LineWidth', 1.8, 'DisplayName', sprintf('Oracle Auto (RMSE %.3f)', a_met.rmse_m));
end
xlabel('East (m)');
ylabel('North (m)');
title(sprintf('Trajectory Compare - %s', letter));
legend('Location', 'best');

nexttile;
hold on; grid on;
plot(k_conf, 'Color', [0.85 0.20 0.20], 'LineWidth', 1.5, 'DisplayName', 'Known confidence');
if run_auto
    plot(a_conf, 'Color', [0.10 0.60 0.20], 'LineWidth', 1.5, 'DisplayName', 'Auto confidence');
end
yline(0.5, '--k', 'LineWidth', 0.8, 'HandleVisibility', 'off');
ylim([0, 1.05]);
xlabel('Sample');
ylabel('Confidence');
title('Confidence');
legend('Location', 'best');

exportgraphics(f, out_path, 'Resolution', 150);
if ~figure_visible
    close(f);
end
end

% -------------------------------------------------------------------------
function dst = merge_cfg(dst, src)
if isempty(src) || ~isstruct(src)
    return;
end
keys = fieldnames(src);
for i = 1:numel(keys)
    k = keys{i};
    if isstruct(src.(k))
        if ~isfield(dst, k) || ~isstruct(dst.(k))
            dst.(k) = src.(k);
        else
            dst.(k) = merge_cfg(dst.(k), src.(k));
        end
    else
        dst.(k) = src.(k);
    end
end
end
