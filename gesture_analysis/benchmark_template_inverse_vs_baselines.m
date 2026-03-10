function summary_tbl = benchmark_template_inverse_vs_baselines(obs_filepath, nav_filepath, letters)
% BENCHMARK_TEMPLATE_INVERSE_VS_BASELINES
% Synthetic benchmark for:
%   1) Baseline: run_gesture_analysis_continuous_track_line
%   2) TemplateInverse: run_gesture_analysis_template_inverse

if nargin < 1 || isempty(obs_filepath)
    obs_filepath = fullfile('data', '1_8', 'A_1_8_1.obs');
end
if nargin < 2 || isempty(nav_filepath)
    nav_filepath = fullfile('data', '1_8', '2026_1_8.nav');
end
if nargin < 3 || isempty(letters)
    letters = {'A', 'B', 'M', 'Star', 'L', 'X', 'Z', 'N'};
end
if ischar(letters) || isstring(letters)
    letters = cellstr(letters);
end

this_dir = fileparts(mfilename('fullpath'));
repo_dir = fileparts(this_dir);
addpath(genpath(repo_dir));

obs_filepath = resolve_data_path(obs_filepath, repo_dir);
nav_filepath = resolve_data_path(nav_filepath, repo_dir);

fprintf('\n== TemplateInverse Benchmark ==\n');
fprintf('obs: %s\n', obs_filepath);
fprintf('nav: %s\n', nav_filepath);
fprintf('letters: %s\n\n', strjoin(upper(letters), ', '));

out_root = fullfile(this_dir, 'results');
if ~exist(out_root, 'dir')
    mkdir(out_root);
end
stamp = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
out_dir = fullfile(out_root, ['template_inverse_', stamp]);
mkdir(out_dir);

old_fig_vis = get(0, 'DefaultFigureVisible');
cleanup_obj = onCleanup(@() set(0, 'DefaultFigureVisible', old_fig_vis)); %#ok<NASGU>
set(0, 'DefaultFigureVisible', 'off');

obs_base = parse_rinex_obs(obs_filepath);
nav_data = parse_rinex_nav_multi_gnss(nav_filepath);

n_case = numel(letters);
results = repmat(struct( ...
    'letter', '', ...
    'old_rmse_m', NaN, ...
    'old_shift_samples', NaN, ...
    'old_start_err_m', NaN, ...
    'old_end_err_m', NaN, ...
    'old_coverage', NaN, ...
    'new_rmse_m', NaN, ...
    'new_shift_samples', NaN, ...
    'new_start_err_m', NaN, ...
    'new_end_err_m', NaN, ...
    'new_coverage', NaN, ...
    'new_err_pct', NaN, ...
    'target20_pass', false, ...
    'best_letter', '', ...
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
    [obs_aligned, step1_aligned] = gesture_analysis_align_duration(obs_waveform, step1_res_shaped);

    [old_x, old_y, old_t, ~] = run_gesture_analysis_continuous_track_line(obs_aligned, nav_data, step1_aligned);

    cfg = struct();
    cfg.debug = struct('verbose', false, 'plot', false);
    [new_x, new_y, new_t, new_conf, dbg] = ...
        run_gesture_analysis_template_inverse(obs_waveform, nav_data, step1_res_shaped, cfg);

    N = numel(step1_res.t_grid);
    [gt_x, gt_y, gt_pen_down] = build_ground_truth(letter, N);

    [old_fx, old_fy] = to_full_series(old_x, old_y, old_t, N);
    [new_fx, new_fy] = to_full_series(new_x, new_y, new_t, N);
    [new_cfull, ~] = to_full_series(new_conf, new_conf, new_t, N);

    old_met = evaluate_against_gt(old_fx, old_fy, gt_x, gt_y, gt_pen_down, 75);
    new_met = evaluate_against_gt(new_fx, new_fy, gt_x, gt_y, gt_pen_down, 75);
    err_pct = estimate_error_percent(letter, new_met.rmse_m);

    fig_path = fullfile(out_dir, sprintf('compare_%s.png', letter));
    make_compare_plot(letter, gt_x, gt_y, gt_pen_down, ...
        old_fx, old_fy, old_met, new_fx, new_fy, new_cfull, new_met, fig_path);

    fprintf('    Baseline RMSE=%.4f | TemplateInverse RMSE=%.4f | ErrPct=%.2f%%\n', ...
        old_met.rmse_m, new_met.rmse_m, err_pct);

    results(i).letter = letter;
    results(i).old_rmse_m = old_met.rmse_m;
    results(i).old_shift_samples = old_met.best_shift_samples;
    results(i).old_start_err_m = old_met.start_err_m;
    results(i).old_end_err_m = old_met.end_err_m;
    results(i).old_coverage = old_met.coverage;
    results(i).new_rmse_m = new_met.rmse_m;
    results(i).new_shift_samples = new_met.best_shift_samples;
    results(i).new_start_err_m = new_met.start_err_m;
    results(i).new_end_err_m = new_met.end_err_m;
    results(i).new_coverage = new_met.coverage;
    results(i).new_err_pct = err_pct;
    results(i).target20_pass = (err_pct <= 20);
    if isstruct(dbg) && isfield(dbg, 'best_letter')
        results(i).best_letter = string(dbg.best_letter);
    else
        results(i).best_letter = "";
    end
    results(i).plot_path = fig_path;
end

summary_tbl = struct2table(results);
csv_path = fullfile(out_dir, 'summary.csv');
writetable(summary_tbl, csv_path);

save(fullfile(out_dir, 'summary.mat'), 'summary_tbl', 'results', 'out_dir');
fprintf('\nBenchmark complete. Output: %s\n', out_dir);
disp(summary_tbl(:, {'letter', 'new_rmse_m', 'new_err_pct', 'target20_pass', 'best_letter'}));
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

best_rmse = inf; best_shift = 0;
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
function make_compare_plot(letter, gt_x, gt_y, gt_pen_down, ...
    old_fx, old_fy, old_met, new_fx, new_fy, new_cfull, new_met, fig_path)

valid_gt = find(gt_pen_down & isfinite(gt_x) & isfinite(gt_y));
figure('Color', 'w', 'Position', [60, 60, 1380, 620]);

subplot(1, 2, 1);
hold on; axis equal; grid on;
plot(0, 0, '^k', 'MarkerFaceColor', 'k', 'DisplayName', 'Receiver');
plot(0, -1, 'bs', 'MarkerFaceColor', [0.2 0.45 0.9], 'DisplayName', 'Body');
plot(gt_x(valid_gt), gt_y(valid_gt), 'k-', 'LineWidth', 3, 'DisplayName', 'GroundTruth');
plot(old_fx, old_fy, 'r-', 'LineWidth', 1.8, 'DisplayName', 'Baseline');
plot(new_fx, new_fy, 'b-', 'LineWidth', 2.2, 'DisplayName', 'TemplateInverse');
xlabel('East (m)'); ylabel('North (m)');
title(sprintf('Letter %s: Trajectory Overlay', upper(letter)));
legend('Location', 'southeast');

subplot(1, 2, 2);
hold on; grid on;
plot(gt_x, 'k-', 'LineWidth', 2.4, 'DisplayName', 'GT X');
plot(old_fx, 'r-', 'LineWidth', 1.5, 'DisplayName', 'Baseline X');
plot(new_fx, 'b-', 'LineWidth', 1.8, 'DisplayName', 'TemplateInverse X');
plot(new_cfull, 'g-', 'LineWidth', 1.6, 'DisplayName', 'TemplateInverse Conf');
xlabel('Sample Index');
ylabel('Value');
title('X(t) and Inverse Confidence');
legend('Location', 'best');

annotation('textbox', [0.12, 0.01, 0.78, 0.1], 'EdgeColor', 'none', ...
    'String', sprintf(['Baseline: RMSE=%.3fm, Shift=%d, Start=%.3fm, End=%.3fm, Coverage=%.2f\n' ...
    'TemplateInverse: RMSE=%.3fm, Shift=%d, Start=%.3fm, End=%.3fm, Coverage=%.2f'], ...
    old_met.rmse_m, old_met.best_shift_samples, old_met.start_err_m, old_met.end_err_m, old_met.coverage, ...
    new_met.rmse_m, new_met.best_shift_samples, new_met.start_err_m, new_met.end_err_m, new_met.coverage), ...
    'HorizontalAlignment', 'center', 'FontSize', 12);

drawnow;
exportgraphics(gcf, fig_path, 'Resolution', 180);
close(gcf);
end
