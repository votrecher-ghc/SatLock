function [summary_tbl, out_dir] = run_research_plot_suite(user_cfg)
% RUN_RESEARCH_PLOT_SUITE
% Standalone research plotting suite for trajectory comparison and paper figures.

if nargin < 1 || isempty(user_cfg)
    user_cfg = struct();
end

repo_dir = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(repo_dir));

cfg = default_cfg(repo_dir);
cfg = merge_cfg(cfg, user_cfg);
cfg.template_names = canonicalize_templates(cfg.template_names);
cfg.sat_curve.template_names = canonicalize_templates(cfg.sat_curve.template_names);

[out_dir, traj_dir, paper_dir] = prepare_output_dirs(repo_dir);

fprintf('\n== Research Plot Suite ==\n');
fprintf('obs: %s\n', cfg.obs_filepath);
fprintf('nav: %s\n', cfg.nav_filepath);
fprintf('templates: %s\n', strjoin(cfg.template_names, ', '));
fprintf('output: %s\n\n', out_dir);

obs_base = parse_rinex_obs(cfg.obs_filepath);
nav_data = parse_rinex_nav_multi_gnss(cfg.nav_filepath);

[cases, rows] = run_genuine_cases(obs_base, nav_data, cfg, traj_dir);
summary_tbl = struct2table(rows);
writetable(summary_tbl, fullfile(out_dir, 'summary_metrics.csv'));

plot_all_templates_overview(cases, fullfile(out_dir, 'all_gestures_compare.png'), cfg);
plot_plane_error_curve(cases, paper_dir, cfg);
plot_feature_space(obs_base, nav_data, paper_dir, cfg);
plot_satellite_count_curve(obs_base, nav_data, paper_dir, cfg);

save(fullfile(out_dir, 'research_plot_suite.mat'), 'summary_tbl', 'cases', 'cfg', 'out_dir');
end

function cfg = default_cfg(repo_dir)
cfg = struct();
cfg.obs_filepath = fullfile(repo_dir, 'data', '1_8', 'A_1_8_1.obs');
cfg.nav_filepath = fullfile(repo_dir, 'data', '1_8', '2026_1_8.nav');
cfg.template_names = default_templates();
cfg.random_seed = 20260321;
cfg.show_figures = false;
cfg.save_resolution = 260;

cfg.span_cfg = struct();
cfg.span_cfg.max_span_x = 0.50;
cfg.span_cfg.max_span_y = 0.50;

cfg.sim_cfg = struct();
cfg.sim_cfg.enable = true;
cfg.sim_cfg.plot = false;
cfg.sim_cfg.max_span_x = cfg.span_cfg.max_span_x;
cfg.sim_cfg.max_span_y = cfg.span_cfg.max_span_y;

cfg.inverse_cfg = struct();
cfg.inverse_cfg.debug = struct('verbose', false, 'plot', false);

cfg.data_cfg = struct();
cfg.data_cfg.debug = struct('verbose', false, 'plot', false);

cfg.embedding = struct();
cfg.embedding.enable = true;
cfg.embedding.attack_modes = {'genuine', 'REPLAY', 'SDR'};
cfg.embedding.max_templates = inf;
cfg.embedding.perplexity = 8;

cfg.sat_curve = struct();
cfg.sat_curve.enable = true;
cfg.sat_curve.template_names = {'A', 'C', 'M', 'Star', 'Rectangle'};
cfg.sat_curve.keep_counts = [4 6 8 10 12];
cfg.sat_curve.num_trials = 1;
cfg.sat_curve.max_templates = 5;
end

function [out_dir, traj_dir, paper_dir] = prepare_output_dirs(repo_dir)
stamp = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss_SSS'));
base_dir = fullfile(repo_dir, 'gesture_analysis', 'results');
out_dir = fullfile(base_dir, ['gesture_research_', stamp]);
traj_dir = fullfile(out_dir, 'by_template');
paper_dir = fullfile(out_dir, 'paper_figures');

ensure_dir(out_dir);
ensure_dir(traj_dir);
ensure_dir(paper_dir);
end

function [cases, rows] = run_genuine_cases(obs_base, nav_data, cfg, traj_dir)
n_case = numel(cfg.template_names);
cases = repmat(struct( ...
    'template', '', ...
    't_grid', [], ...
    'gt_x', [], ...
    'gt_y', [], ...
    'gt_pen', [], ...
    'num_visible_sats', NaN, ...
    'inverse', empty_alg_case(), ...
    'data_driven', empty_alg_case()), n_case, 1);
rows = repmat(empty_row(), n_case * 2, 1);
row_idx = 0;

for i = 1:n_case
    template_name = cfg.template_names{i};
    rng(cfg.random_seed + i, 'twister');

    obs_sim = simulate_template(obs_base, nav_data, template_name, cfg.sim_cfg);
    [~, step1_res, obs_waveform, step1_res_shaped, ~, ~] = run_preprocess_pipeline(obs_sim);
    t_grid = resolve_t_grid(step1_res, step1_res_shaped);
    [gt_x, gt_y, gt_pen] = build_ground_truth_local(template_name, numel(t_grid), cfg.span_cfg);

    inv_case = run_algorithm_case('inverse_beam_legacy', obs_waveform, nav_data, step1_res_shaped, t_grid, gt_x, gt_y, gt_pen, cfg);
    data_case = run_algorithm_case('data_driven_new', obs_waveform, nav_data, step1_res_shaped, t_grid, gt_x, gt_y, gt_pen, cfg);

    cases(i).template = template_name;
    cases(i).t_grid = t_grid;
    cases(i).gt_x = gt_x;
    cases(i).gt_y = gt_y;
    cases(i).gt_pen = gt_pen;
    cases(i).num_visible_sats = numel(step1_res_shaped.valid_sats);
    cases(i).inverse = inv_case;
    cases(i).data_driven = data_case;

    row_idx = row_idx + 1;
    rows(row_idx) = make_row(template_name, "inverse_beam_legacy", inv_case, cases(i).num_visible_sats);
    row_idx = row_idx + 1;
    rows(row_idx) = make_row(template_name, "data_driven_new", data_case, cases(i).num_visible_sats);

    save_template_compare(cases(i), fullfile(traj_dir, sprintf('%s_compare.png', template_name)), cfg);
end

rows = rows(1:row_idx);
end

function obs_sim = simulate_template(obs_base, nav_data, template_name, sim_cfg)
sim_cfg_local = sim_cfg;
sim_cfg_local.enable = true;
sim_cfg_local.target_letter = template_name;
sim_cfg_local.plot = false;
obs_sim = generate_ideal_multi_shape(obs_base, nav_data, template_name, sim_cfg_local);
end

function t_grid = resolve_t_grid(step1_res, step1_res_shaped)
if isfield(step1_res_shaped, 't_grid') && ~isempty(step1_res_shaped.t_grid)
    t_grid = step1_res_shaped.t_grid;
else
    t_grid = step1_res.t_grid;
end
end

function alg_case = run_algorithm_case(alg_name, obs_waveform, nav_data, step1_res_shaped, t_grid, gt_x, gt_y, gt_pen, cfg)
alg_case = empty_alg_case();
alg_case.name = alg_name;

try
    switch alg_name
        case 'inverse_beam_legacy'
            [x, y, t, conf] = run_gesture_analysis_inverse_beam( ...
                obs_waveform, nav_data, step1_res_shaped, cfg.inverse_cfg);
        case 'data_driven_new'
            [x, y, t, conf] = run_gesture_analysis_data_driven( ...
                obs_waveform, nav_data, step1_res_shaped, cfg.data_cfg);
        otherwise
            error('Unsupported algorithm: %s', alg_name);
    end
    alg_case.status = "ok";
catch ME
    warning('run_research_plot_suite:AlgorithmFailure', '%s failed: %s', alg_name, ME.message);
    x = [];
    y = [];
    t = [];
    conf = [];
    alg_case.status = "failed";
    alg_case.error_message = string(ME.message);
end

alg_case.x = x;
alg_case.y = y;
alg_case.t = t;
alg_case.conf = conf;
[alg_case.plot_x, alg_case.plot_y] = order_plot_series(x, y, t, t_grid);
[alg_case.full_x, alg_case.full_y] = to_full_series(x, y, t, numel(t_grid), t_grid);
alg_case.metrics = evaluate_reconstruction(alg_case.full_x, alg_case.full_y, gt_x, gt_y, gt_pen, 75);
alg_case.metrics.mean_conf = mean(conf, 'omitnan');
alg_case.feature_vector = feature_vector_from_metrics(alg_case.metrics);
end

function save_template_compare(case_item, out_path, cfg)
fig_vis = on_off(cfg.show_figures);
f = figure('Visible', fig_vis, 'Color', 'w', 'Position', [100, 80, 1000, 760]);
ax = axes(f);
hold(ax, 'on');
grid(ax, 'on');
axis(ax, 'equal');

idx_up = ~case_item.gt_pen & isfinite(case_item.gt_x) & isfinite(case_item.gt_y);
idx_dn = case_item.gt_pen & isfinite(case_item.gt_x) & isfinite(case_item.gt_y);
plot(ax, case_item.gt_x(idx_up), case_item.gt_y(idx_up), '--', ...
    'Color', [0.70 0.70 0.70], 'LineWidth', 1.0, 'DisplayName', 'GT pen-up');
plot(ax, case_item.gt_x(idx_dn), case_item.gt_y(idx_dn), '-', ...
    'Color', [0.10 0.45 0.95], 'LineWidth', 2.8, 'DisplayName', 'GT pen-down');
plot(ax, case_item.inverse.plot_x, case_item.inverse.plot_y, '-', ...
    'Color', [0.83 0.18 0.17], 'LineWidth', 2.0, 'DisplayName', 'Inverse Beam');
plot(ax, case_item.data_driven.plot_x, case_item.data_driven.plot_y, '-', ...
    'Color', [0.12 0.12 0.12], 'LineWidth', 2.0, 'DisplayName', 'Data-Driven');

title(ax, sprintf('Template %s', case_item.template), 'Interpreter', 'none');
xlabel(ax, 'East (m)');
ylabel(ax, 'North (m)');
style_axes(ax, 11);
fit_axis(ax, case_item.gt_x, case_item.gt_y, case_item.inverse.plot_x, case_item.inverse.plot_y);
fit_axis(ax, case_item.gt_x, case_item.gt_y, case_item.data_driven.plot_x, case_item.data_driven.plot_y);

legend(ax, 'Location', 'southoutside', 'Orientation', 'horizontal');
summary_text = sprintf([ ...
    'Inverse: RMSE %.3f  DTW %.3f  Cov %.2f\n', ...
    'Data:    RMSE %.3f  DTW %.3f  Cov %.2f\n', ...
    'Visible Sats: %d'], ...
    case_item.inverse.metrics.rmse_m, case_item.inverse.metrics.dtw_m, case_item.inverse.metrics.coverage, ...
    case_item.data_driven.metrics.rmse_m, case_item.data_driven.metrics.dtw_m, case_item.data_driven.metrics.coverage, ...
    round(case_item.num_visible_sats));
annotate_box(ax, summary_text);
save_figure(f, out_path, cfg.save_resolution, cfg.show_figures);
end

function plot_all_templates_overview(cases, out_path, cfg)
if isempty(cases)
    return;
end

n_case = numel(cases);
n_col = 4;
n_row = ceil(n_case / n_col);
f = figure('Visible', on_off(cfg.show_figures), 'Color', 'w', 'Position', [60, 60, 1650, 940]);
tl = tiledlayout(f, n_row, n_col, 'TileSpacing', 'compact', 'Padding', 'compact');

for i = 1:n_case
    ax = nexttile(tl);
    hold(ax, 'on');
    grid(ax, 'on');
    axis(ax, 'equal');

    c = cases(i);
    idx_up = ~c.gt_pen & isfinite(c.gt_x) & isfinite(c.gt_y);
    idx_dn = c.gt_pen & isfinite(c.gt_x) & isfinite(c.gt_y);
    plot(ax, c.gt_x(idx_up), c.gt_y(idx_up), '--', 'Color', [0.72 0.72 0.72], 'LineWidth', 0.8);
    plot(ax, c.gt_x(idx_dn), c.gt_y(idx_dn), '-', 'Color', [0.10 0.45 0.95], 'LineWidth', 2.0);
    plot(ax, c.inverse.plot_x, c.inverse.plot_y, '-', 'Color', [0.83 0.18 0.17], 'LineWidth', 1.5);
    plot(ax, c.data_driven.plot_x, c.data_driven.plot_y, '-', 'Color', [0.12 0.12 0.12], 'LineWidth', 1.5);
    title(ax, c.template, 'Interpreter', 'none');
    style_axes(ax, 9);
    fit_axis(ax, c.gt_x, c.gt_y, c.inverse.plot_x, c.inverse.plot_y);
    fit_axis(ax, c.gt_x, c.gt_y, c.data_driven.plot_x, c.data_driven.plot_y);
end

title(tl, 'All Gesture Templates: GT vs Inverse Beam vs Data-Driven');
save_figure(f, out_path, cfg.save_resolution, cfg.show_figures);
end

function plot_plane_error_curve(cases, paper_dir, cfg)
curve_n = 100;
inv_curves = nan(numel(cases), curve_n);
data_curves = nan(numel(cases), curve_n);

for i = 1:numel(cases)
    inv_curves(i, :) = resample_error_curve(cases(i).inverse.metrics.point_errors_m, curve_n);
    data_curves(i, :) = resample_error_curve(cases(i).data_driven.metrics.point_errors_m, curve_n);
end

f = figure('Visible', on_off(cfg.show_figures), 'Color', 'w', 'Position', [120, 80, 980, 720]);
ax = axes(f);
hold(ax, 'on');
grid(ax, 'on');
style_axes(ax, 11);

x = linspace(0, 100, curve_n);
plot_mean_band(ax, x, inv_curves, [0.83 0.18 0.17], 'Inverse Beam');
plot_mean_band(ax, x, data_curves, [0.12 0.12 0.12], 'Data-Driven');
legend(ax, 'Location', 'northeast');
xlabel(ax, 'Normalized Gesture Progress (%)');
ylabel(ax, 'Hand-Plane Error (m)');
title(ax, 'Gesture Plane Error Curve');

save_figure(f, fullfile(paper_dir, 'gesture_plane_error_curve.png'), cfg.save_resolution, cfg.show_figures);
end

function plot_feature_space(obs_base, nav_data, paper_dir, cfg)
if ~cfg.embedding.enable
    return;
end

template_names = cfg.template_names;
if isfinite(cfg.embedding.max_templates)
    template_names = template_names(1:min(numel(template_names), cfg.embedding.max_templates));
end

features = zeros(0, 12);
labels = strings(0, 1);
algorithms = strings(0, 1);

for i = 1:numel(template_names)
    template_name = template_names{i};
    rng(cfg.random_seed + 100 + i, 'twister');

    obs_genuine = simulate_template(obs_base, nav_data, template_name, cfg.sim_cfg);
    for m = 1:numel(cfg.embedding.attack_modes)
        mode_name = cfg.embedding.attack_modes{m};
        if strcmpi(mode_name, 'genuine')
            obs_case = obs_genuine;
        else
            obs_case = simulate_gnss_spoofing(obs_genuine, nav_data, upper(mode_name));
        end

        [~, step1_res, obs_waveform, step1_res_shaped, ~, ~] = run_preprocess_pipeline(obs_case);
        t_grid = resolve_t_grid(step1_res, step1_res_shaped);
        [gt_x, gt_y, gt_pen] = build_ground_truth_local(template_name, numel(t_grid), cfg.span_cfg);

        inv_case = run_algorithm_case('inverse_beam_legacy', obs_waveform, nav_data, step1_res_shaped, t_grid, gt_x, gt_y, gt_pen, cfg);
        data_case = run_algorithm_case('data_driven_new', obs_waveform, nav_data, step1_res_shaped, t_grid, gt_x, gt_y, gt_pen, cfg);

        if inv_case.status == "ok" && isfinite(inv_case.metrics.rmse_m) && all(isfinite(inv_case.feature_vector))
            features(end + 1, :) = inv_case.feature_vector; %#ok<AGROW>
            labels(end + 1, 1) = string(mode_name); %#ok<AGROW>
            algorithms(end + 1, 1) = "Inverse Beam"; %#ok<AGROW>
        end
        if data_case.status == "ok" && isfinite(data_case.metrics.rmse_m) && all(isfinite(data_case.feature_vector))
            features(end + 1, :) = data_case.feature_vector; %#ok<AGROW>
            labels(end + 1, 1) = string(mode_name); %#ok<AGROW>
            algorithms(end + 1, 1) = "Data-Driven"; %#ok<AGROW>
        end
    end
end

if size(features, 1) < 6
    warning('run_research_plot_suite:EmbeddingSkipped', 'Not enough samples for PCA/t-SNE.');
    return;
end

z = zscore(features);
z(:, all(~isfinite(z), 1)) = 0;
z(~isfinite(z)) = 0;

[~, score, ~, ~, explained] = pca(z);
plot_embedding(score(:, 1:2), labels, algorithms, ...
    sprintf('Feature Space PCA (PC1 %.1f%%, PC2 %.1f%%)', explained(1), explained(2)), ...
    fullfile(paper_dir, 'feature_space_pca.png'), cfg);

perplexity = min(cfg.embedding.perplexity, max(3, floor((size(z, 1) - 1) / 3)));
map = tsne(z, 'NumDimensions', 2, 'Perplexity', perplexity, 'Standardize', false);
plot_embedding(map, labels, algorithms, 'Feature Space t-SNE', ...
    fullfile(paper_dir, 'feature_space_tsne.png'), cfg);
end

function plot_satellite_count_curve(obs_base, nav_data, paper_dir, cfg)
if ~cfg.sat_curve.enable
    return;
end

template_names = cfg.sat_curve.template_names;
if isfinite(cfg.sat_curve.max_templates)
    template_names = template_names(1:min(numel(template_names), cfg.sat_curve.max_templates));
end

keep_counts = cfg.sat_curve.keep_counts(:).';
inv_vals = nan(numel(template_names), numel(keep_counts), cfg.sat_curve.num_trials);
data_vals = nan(numel(template_names), numel(keep_counts), cfg.sat_curve.num_trials);

for i = 1:numel(template_names)
    template_name = template_names{i};
    rng(cfg.random_seed + 200 + i, 'twister');
    obs_sim = simulate_template(obs_base, nav_data, template_name, cfg.sim_cfg);
    sat_list = collect_valid_sats(obs_sim);
    if numel(sat_list) < min(keep_counts)
        continue;
    end

    for c = 1:numel(keep_counts)
        keep_n = keep_counts(c);
        if keep_n > numel(sat_list)
            continue;
        end

        for trial = 1:cfg.sat_curve.num_trials
            rng(cfg.random_seed + 3000 + 100 * i + 10 * c + trial, 'twister');
            ord = randperm(numel(sat_list));
            keep_sats = sat_list(ord(1:keep_n));
            drop_sats = setdiff(sat_list, keep_sats, 'stable');
            obs_sub = del_sat(obs_sim, drop_sats);

            [~, step1_res, obs_waveform, step1_res_shaped, ~, ~] = run_preprocess_pipeline(obs_sub);
            t_grid = resolve_t_grid(step1_res, step1_res_shaped);
            [gt_x, gt_y, gt_pen] = build_ground_truth_local(template_name, numel(t_grid), cfg.span_cfg);

            inv_case = run_algorithm_case('inverse_beam_legacy', obs_waveform, nav_data, step1_res_shaped, t_grid, gt_x, gt_y, gt_pen, cfg);
            data_case = run_algorithm_case('data_driven_new', obs_waveform, nav_data, step1_res_shaped, t_grid, gt_x, gt_y, gt_pen, cfg);
            if inv_case.status == "ok" && isfinite(inv_case.metrics.rmse_m)
                inv_vals(i, c, trial) = inv_case.metrics.rmse_m;
            end
            if data_case.status == "ok" && isfinite(data_case.metrics.rmse_m)
                data_vals(i, c, trial) = data_case.metrics.rmse_m;
            end
        end
    end
end

f = figure('Visible', on_off(cfg.show_figures), 'Color', 'w', 'Position', [120, 80, 980, 720]);
ax = axes(f);
hold(ax, 'on');
grid(ax, 'on');
style_axes(ax, 11);

plot_mean_band(ax, keep_counts, reshape(inv_vals, [], numel(keep_counts)), [0.83 0.18 0.17], 'Inverse Beam');
plot_mean_band(ax, keep_counts, reshape(data_vals, [], numel(keep_counts)), [0.12 0.12 0.12], 'Data-Driven');
legend(ax, 'Location', 'northeast');
xlabel(ax, 'Number of Retained Satellites');
ylabel(ax, 'RMSE (m)');
title(ax, 'Impact of Visible Satellite Count');

save_figure(f, fullfile(paper_dir, 'satellite_count_impact.png'), cfg.save_resolution, cfg.show_figures);
end

function plot_embedding(map, labels, algorithms, title_str, out_path, cfg)
f = figure('Visible', on_off(cfg.show_figures), 'Color', 'w', 'Position', [120, 80, 980, 720]);
ax = axes(f);
hold(ax, 'on');
grid(ax, 'on');
style_axes(ax, 11);

classes = unique(labels, 'stable');
colors = lines(numel(classes));
for i = 1:numel(classes)
    mask_class = labels == classes(i);
    for alg_name = ["Inverse Beam", "Data-Driven"]
        mask = mask_class & algorithms == alg_name;
        if ~any(mask)
            continue;
        end
        scatter(ax, map(mask, 1), map(mask, 2), 58, ...
            'Marker', marker_for_algorithm(alg_name), ...
            'MarkerFaceColor', colors(i, :), ...
            'MarkerEdgeColor', 'k', ...
            'DisplayName', sprintf('%s - %s', classes(i), alg_name));
    end
end

title(ax, title_str, 'Interpreter', 'none');
xlabel(ax, 'Dim 1');
ylabel(ax, 'Dim 2');
legend(ax, 'Location', 'eastoutside');
save_figure(f, out_path, cfg.save_resolution, cfg.show_figures);
end

function plot_mean_band(ax, x, samples, color_rgb, label_name)
if isempty(samples)
    return;
end
mu = mean(samples, 1, 'omitnan');
sd = std(samples, 0, 1, 'omitnan');
valid = isfinite(mu);
if ~any(valid)
    return;
end

xv = x(valid);
mu = mu(valid);
sd = sd(valid);
fill(ax, [xv, fliplr(xv)], [mu - sd, fliplr(mu + sd)], color_rgb, ...
    'FaceAlpha', 0.14, 'EdgeColor', 'none', 'HandleVisibility', 'off');
plot(ax, xv, mu, '-o', 'LineWidth', 2.0, 'Color', color_rgb, ...
    'MarkerFaceColor', color_rgb, 'MarkerSize', 6, 'DisplayName', label_name);
end

function curve = resample_error_curve(err, n_out)
curve = nan(1, n_out);
err = err(:);
err = err(isfinite(err));
if isempty(err)
    return;
end
if numel(err) == 1
    curve(:) = err;
    return;
end
x = linspace(0, 1, numel(err));
xi = linspace(0, 1, n_out);
curve = interp1(x, err, xi, 'linear');
end

function [gt_x, gt_y, gt_pen] = build_ground_truth_local(template_name, n_samples, span_cfg)
if exist('gesture_template_library', 'file') == 2
    [gt_x, gt_y, gt_pen] = gesture_template_library('groundtruth', template_name, n_samples, span_cfg);
    return;
end
error('run_research_plot_suite:MissingTemplateLibrary', ...
    'gesture_template_library.m is required for the research plotting suite.');
end

function [plot_x, plot_y] = order_plot_series(x, y, t_idx, t_grid)
plot_x = [];
plot_y = [];
if isempty(x) || isempty(y)
    return;
end
x = x(:);
y = y(:);
n = min(numel(x), numel(y));
x = x(1:n);
y = y(1:n);
if isempty(t_idx)
    keep = isfinite(x) & isfinite(y);
    plot_x = x(keep);
    plot_y = y(keep);
    return;
end
t_idx = normalize_time_index(t_idx, numel(t_grid), t_grid);
t_idx = t_idx(1:min(numel(t_idx), n));
n = min([numel(x), numel(y), numel(t_idx)]);
x = x(1:n);
y = y(1:n);
t_idx = t_idx(1:n);
keep = isfinite(x) & isfinite(y) & isfinite(t_idx);
x = x(keep);
y = y(keep);
t_idx = t_idx(keep);
[~, ord] = sort(t_idx, 'ascend');
plot_x = x(ord);
plot_y = y(ord);
end

function [fx, fy] = to_full_series(x, y, t_idx, N, t_grid)
fx = nan(N, 1);
fy = nan(N, 1);
if isempty(x) || isempty(y) || isempty(t_idx)
    return;
end
x = x(:);
y = y(:);
t_idx = normalize_time_index(t_idx, N, t_grid);
n = min([numel(x), numel(y), numel(t_idx)]);
x = x(1:n);
y = y(1:n);
t_idx = t_idx(1:n);
keep = isfinite(x) & isfinite(y) & isfinite(t_idx) & t_idx >= 1 & t_idx <= N;
x = x(keep);
y = y(keep);
t_idx = t_idx(keep);
if isempty(t_idx)
    return;
end
[u_idx, ia] = unique(t_idx, 'stable');
fx(u_idx) = x(ia);
fy(u_idx) = y(ia);
end

function idx = normalize_time_index(t, N, t_grid)
if isnumeric(t)
    idx = round(t(:));
    if all(isfinite(idx)) && all(idx >= 1) && all(idx <= N)
        return;
    end
end
if isdatetime(t)
    tg = posixtime(t_grid(:));
    tt = posixtime(t(:));
    idx = round(interp1(tg, 1:N, tt, 'nearest', 'extrap'));
    idx = min(max(idx, 1), N);
    return;
end
idx = round(linspace(1, N, numel(t)))';
end

function met = evaluate_reconstruction(full_x, full_y, gt_x, gt_y, gt_pen, max_shift)
met = struct( ...
    'rmse_m', inf, ...
    'mte_m', inf, ...
    'dtw_m', inf, ...
    'start_err_m', inf, ...
    'end_err_m', inf, ...
    'coverage', 0, ...
    'point_errors_m', [], ...
    'aligned_est_x', [], ...
    'aligned_est_y', [], ...
    'aligned_gt_x', [], ...
    'aligned_gt_y', [], ...
    'path_length_m', NaN, ...
    'x_span_m', NaN, ...
    'y_span_m', NaN, ...
    'mean_conf', NaN);

idx_gt = find(gt_pen & isfinite(gt_x) & isfinite(gt_y));
if isempty(idx_gt)
    return;
end

valid_est = find(isfinite(full_x) & isfinite(full_y));
if isempty(valid_est)
    return;
end

N = numel(gt_x);
met.coverage = nnz(isfinite(full_x(idx_gt)) & isfinite(full_y(idx_gt))) / numel(idx_gt);
if numel(valid_est) >= 2
    est_x = interp1(valid_est, full_x(valid_est), 1:N, 'linear', 'extrap').';
    est_y = interp1(valid_est, full_y(valid_est), 1:N, 'linear', 'extrap').';
else
    est_x = repmat(full_x(valid_est(1)), N, 1);
    est_y = repmat(full_y(valid_est(1)), N, 1);
end

best_rmse = inf;
best_shift = 0;
best_keep = [];
for sh = -max_shift:max_shift
    idx_est = idx_gt + sh;
    keep = idx_est >= 1 & idx_est <= N;
    if nnz(keep) < max(8, round(0.2 * numel(idx_gt)))
        continue;
    end
    g_idx = idx_gt(keep);
    e_idx = idx_est(keep);
    err = hypot(est_x(e_idx) - gt_x(g_idx), est_y(e_idx) - gt_y(g_idx));
    rmse = sqrt(mean(err .^ 2));
    if rmse < best_rmse
        best_rmse = rmse;
        best_shift = sh;
        best_keep = keep;
    end
end

if isempty(best_keep)
    return;
end

g_idx = idx_gt(best_keep);
e_idx = g_idx + best_shift;
err = hypot(est_x(e_idx) - gt_x(g_idx), est_y(e_idx) - gt_y(g_idx));

met.rmse_m = sqrt(mean(err .^ 2));
met.mte_m = mean(err);
met.dtw_m = compute_dtw([est_x(e_idx), est_y(e_idx)], [gt_x(g_idx), gt_y(g_idx)]);
met.start_err_m = err(1);
met.end_err_m = err(end);
met.point_errors_m = err;
met.aligned_est_x = est_x(e_idx);
met.aligned_est_y = est_y(e_idx);
met.aligned_gt_x = gt_x(g_idx);
met.aligned_gt_y = gt_y(g_idx);
met.path_length_m = polyline_length(met.aligned_est_x, met.aligned_est_y);
met.x_span_m = span_of(met.aligned_est_x);
met.y_span_m = span_of(met.aligned_est_y);
end

function feat = feature_vector_from_metrics(met)
feat = [ ...
    fallback(met.rmse_m, 1.0), ...
    fallback(met.mte_m, 1.0), ...
    fallback(met.dtw_m, 1.0), ...
    fallback(met.coverage, 0.0), ...
    fallback(met.start_err_m, 1.0), ...
    fallback(met.end_err_m, 1.0), ...
    fallback(met.path_length_m, 0.0), ...
    fallback(met.x_span_m, 0.0), ...
    fallback(met.y_span_m, 0.0), ...
    fallback(mean(met.point_errors_m, 'omitnan'), 1.0), ...
    fallback(std(met.point_errors_m, 0, 'omitnan'), 0.0), ...
    fallback(met.mean_conf, 0.0)];
end

function d = compute_dtw(a, b)
a = sanitize_series(a);
b = sanitize_series(b);
if isempty(a) || isempty(b)
    d = inf;
    return;
end

[ax, ay] = resample_polyline(a(:, 1), a(:, 2), 120);
[bx, by] = resample_polyline(b(:, 1), b(:, 2), 120);
a = [ax, ay];
b = [bx, by];

na = size(a, 1);
nb = size(b, 1);
dp = inf(na + 1, nb + 1);
dp(1, 1) = 0;
for i = 1:na
    for j = 1:nb
        cost = norm(a(i, :) - b(j, :));
        dp(i + 1, j + 1) = cost + min([dp(i, j + 1), dp(i + 1, j), dp(i, j)]);
    end
end
d = dp(end, end) / max(na + nb, 1);
end

function [xr, yr] = resample_polyline(x, y, n_out)
x = x(:);
y = y(:);
keep = isfinite(x) & isfinite(y);
x = x(keep);
y = y(keep);
if isempty(x)
    xr = nan(n_out, 1);
    yr = nan(n_out, 1);
    return;
end
if numel(x) == 1
    xr = repmat(x, n_out, 1);
    yr = repmat(y, n_out, 1);
    return;
end
d = hypot(diff(x), diff(y));
mask = [true; d > eps];
x = x(mask);
y = y(mask);
if numel(x) == 1
    xr = repmat(x, n_out, 1);
    yr = repmat(y, n_out, 1);
    return;
end
s = [0; cumsum(hypot(diff(x), diff(y)))];
sq = linspace(0, s(end), n_out).';
xr = interp1(s, x, sq, 'linear');
yr = interp1(s, y, sq, 'linear');
end

function v = span_of(x)
x = x(isfinite(x));
if isempty(x)
    v = NaN;
else
    v = max(x) - min(x);
end
end

function L = polyline_length(x, y)
x = x(:);
y = y(:);
keep = isfinite(x) & isfinite(y);
x = x(keep);
y = y(keep);
if numel(x) < 2
    L = 0;
else
    L = sum(hypot(diff(x), diff(y)), 'omitnan');
end
end

function sat_list = collect_valid_sats(obs_data)
sat_list = {};
scan_n = min(120, numel(obs_data));
for i = 1:scan_n
    if isempty(obs_data(i).data)
        continue;
    end
    sat_list = [sat_list, fieldnames(obs_data(i).data)']; %#ok<AGROW>
end
sat_list = unique(sat_list, 'stable');
keep = false(size(sat_list));
for i = 1:numel(sat_list)
    keep(i) = ismember(sat_list{i}(1), ['G', 'R', 'C', 'E', 'J']);
end
sat_list = sat_list(keep);
end

function row = make_row(template_name, alg_name, alg_case, num_sats)
row = empty_row();
row.template = string(template_name);
row.algorithm = string(alg_name);
row.rmse_m = alg_case.metrics.rmse_m;
row.mte_m = alg_case.metrics.mte_m;
row.dtw_m = alg_case.metrics.dtw_m;
row.start_err_m = alg_case.metrics.start_err_m;
row.end_err_m = alg_case.metrics.end_err_m;
row.coverage = alg_case.metrics.coverage;
row.visible_satellites = num_sats;
end

function row = empty_row()
row = struct( ...
    'template', "", ...
    'algorithm', "", ...
    'rmse_m', NaN, ...
    'mte_m', NaN, ...
    'dtw_m', NaN, ...
    'start_err_m', NaN, ...
    'end_err_m', NaN, ...
    'coverage', NaN, ...
    'visible_satellites', NaN);
end

function alg_case = empty_alg_case()
alg_case = struct( ...
    'name', "", ...
    'status', "", ...
    'error_message', "", ...
    'x', [], ...
    'y', [], ...
    't', [], ...
    'conf', [], ...
    'plot_x', [], ...
    'plot_y', [], ...
    'full_x', [], ...
    'full_y', [], ...
    'metrics', struct(), ...
    'feature_vector', []);
end

function names = canonicalize_templates(names_in)
if isempty(names_in)
    names = default_templates();
    return;
end
if ischar(names_in) || isstring(names_in)
    names_in = cellstr(names_in);
end
names = cell(size(names_in));
for i = 1:numel(names_in)
    if exist('gesture_template_library', 'file') == 2
        names{i} = gesture_template_library('label', names_in{i});
    else
        names{i} = char(names_in{i});
    end
end
[~, ia] = unique(lower(string(names)), 'stable');
names = names(sort(ia));
end

function names = default_templates()
if exist('gesture_template_library', 'file') == 2
    names = gesture_template_library('all');
else
    names = {'A', 'C', 'M', 'Star', 'L', 'X', 'Z', 'N', 'V', 'Rectangle', 'LeftSwipe', 'RightSwipe'};
end
end

function style_axes(ax, font_size)
set(ax, 'FontName', 'Times New Roman', 'FontSize', font_size, 'LineWidth', 1.0, 'Box', 'on');
end

function annotate_box(ax, txt)
text(ax, 0.98, 0.98, txt, 'Units', 'normalized', ...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', ...
    'BackgroundColor', [1 1 1], 'Margin', 8, 'FontSize', 9);
end

function fit_axis(ax, x1, y1, x2, y2)
all_x = [x1(:); x2(:)];
all_y = [y1(:); y2(:)];
keep = isfinite(all_x) & isfinite(all_y);
all_x = all_x(keep);
all_y = all_y(keep);
if isempty(all_x)
    return;
end
span = max([max(all_x) - min(all_x), max(all_y) - min(all_y), 0.25]);
pad = 0.12 * span;
cx = 0.5 * (max(all_x) + min(all_x));
cy = 0.5 * (max(all_y) + min(all_y));
xlim(ax, [cx - 0.5 * span - pad, cx + 0.5 * span + pad]);
ylim(ax, [cy - 0.5 * span - pad, cy + 0.5 * span + pad]);
end

function save_figure(fig, out_path, resolution, keep_open)
try
    exportgraphics(fig, out_path, 'Resolution', resolution);
catch
    saveas(fig, out_path);
end
if ~keep_open
    close(fig);
end
end

function s = on_off(tf)
if tf
    s = 'on';
else
    s = 'off';
end
end

function ensure_dir(path_in)
if ~exist(path_in, 'dir')
    mkdir(path_in);
end
end

function cfg = merge_cfg(cfg, updates)
if isempty(updates) || ~isstruct(updates)
    return;
end
keys = fieldnames(updates);
for i = 1:numel(keys)
    key = keys{i};
    if isstruct(updates.(key))
        if ~isfield(cfg, key) || ~isstruct(cfg.(key))
            cfg.(key) = updates.(key);
        else
            cfg.(key) = merge_cfg(cfg.(key), updates.(key));
        end
    else
        cfg.(key) = updates.(key);
    end
end
end

function xy = sanitize_series(xy)
if isempty(xy)
    return;
end
if isvector(xy)
    xy = xy(:);
end
xy = xy(all(isfinite(xy), 2), :);
end

function mk = marker_for_algorithm(alg_name)
if alg_name == "Inverse Beam"
    mk = 'o';
else
    mk = 's';
end
end

function v = fallback(v, default_v)
if isempty(v) || ~isfinite(v)
    v = default_v;
end
end
