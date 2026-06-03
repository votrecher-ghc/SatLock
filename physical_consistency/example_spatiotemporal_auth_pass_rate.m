% Example: draw spatiotemporal consistency authentication pass-rate figures.
%
% The authentication rule is:
%   pass if RMS(r_obs - rho_GT) <= tau_ST
%
% tau_ST is calibrated from the authorized spatiotemporal point by the
% selected quantile of per-epoch residual RMS values.

clear;
clc;
close all;

repo_dir = resolve_repo_root_local(fileparts(mfilename('fullpath')));
addpath(genpath(repo_dir));

%% ===================== Edit here =====================
nav_file = fullfile(repo_dir, 'data', '1_8', '2026_1_8.nav');
obs_file = fullfile(repo_dir, 'data', '1_8', 'A_1_8_1.obs');

start_time = datetime(2026, 1, 8, 5, 59, 39.116, 'TimeZone', 'UTC');

% Authorized spatiotemporal point.
lat0_deg = 35.77204196;
lon0_deg = 120.02477210;
height0_m = 92.889;

systems = {'G', 'C', 'E', 'R', 'J'};

% Trial window. Each OBS epoch in this window is treated as one trial.
auth_duration_seconds = 2.0;
auth_max_epochs = 20;
min_epoch_rows = 10;
tau_quantile_pct = 95;

% Spatial pass-rate curve: distance offset step is 5 m.
spatial_offsets_m = 0:5:100;
spatial_direction_deg = [0, 90, 180, 270];

% Temporal pass-rate curve: time offset step is 0.1 s.
time_offsets_s = -5:0.1:5;

% Set false to recompute all pass-rate tables from NAV/OBS.
% Keep true when only updating figure styles/annotations.
reuse_existing_csv = true;

out_dir = fullfile(repo_dir, 'physical_consistency', 'outputs', ...
    'spatiotemporal_uniqueness');
out_pdf_dir = fullfile(out_dir, 'pdf');
out_png_dir = fullfile(out_dir, 'png');
out_fig_dir = fullfile(out_dir, 'fig');
out_csv_dir = fullfile(out_dir, 'csv');

paper_font_name = 'Arial';
paper_tick_font_size = 30;
paper_label_font_size = 30;
paper_axis_line_width = 1.2;
paper_curve_line_width = 3.0;
paper_marker_size = 10;
paper_annotation_font_size = 22;

paper_line_color = [0.19 0.40 0.67];
paper_marker_color = [0.23 0.66 0.88];
paper_threshold_color = [0.62 0.18 0.16];
%% =====================================================

output_dirs = {out_dir, out_pdf_dir, out_png_dir, out_fig_dir, out_csv_dir};
for di = 1:numel(output_dirs)
    if ~exist(output_dirs{di}, 'dir')
        mkdir(output_dirs{di});
    end
end
organize_existing_outputs_local(out_dir, out_png_dir, out_pdf_dir, ...
    out_fig_dir, out_csv_dir);

threshold_csv = fullfile(out_csv_dir, ...
    'spatiotemporal_auth_threshold_calibration.csv');
authorized_epoch_csv = fullfile(out_csv_dir, ...
    'spatiotemporal_authorized_epoch_rms.csv');
spatial_csv = fullfile(out_csv_dir, 'spatial_auth_pass_rate.csv');
temporal_csv = fullfile(out_csv_dir, 'temporal_auth_pass_rate.csv');

if reuse_existing_csv && isfile(threshold_csv) && isfile(spatial_csv) && ...
        isfile(temporal_csv)
    fprintf('Loading cached authentication pass-rate tables...\n');
    threshold_tbl = readtable(threshold_csv);
    spatial_tbl = readtable(spatial_csv);
    temporal_tbl = readtable(temporal_csv);
    tau_st_m = threshold_tbl.tau_st_m(1);
    fprintf('  tau_ST %.3f m loaded from %s\n', tau_st_m, threshold_csv);
else
    fprintf('Parsing NAV/OBS once...\n');
    nav_data = parse_rinex_nav_multi_gnss(nav_file);
    obs_data = parse_rinex_obs(obs_file);

    common_opts = { ...
        'systems', systems, ...
        'elevation_mask_deg', 0, ...
        'min_clock_sats', 4, ...
        'receiver_clock_estimator', 'mean', ...
        'start_time', start_time, ...
        'duration_seconds', auth_duration_seconds, ...
        'max_epochs', auth_max_epochs};

    fprintf('\nCalibrating authorized threshold...\n');
    auth_epoch_tbl = evaluate_candidate_epoch_rms_local(nav_data, obs_data, ...
        lat0_deg, lon0_deg, height0_m, common_opts, 0, min_epoch_rows);

    if height(auth_epoch_tbl) < 3
        error('example_spatiotemporal_auth_pass_rate:NotEnoughAuthorizedTrials', ...
            'Not enough authorized trials for threshold calibration.');
    end

    tau_st_m = percentile_sorted_local(auth_epoch_tbl.rms_m, tau_quantile_pct);
    auth_pass_rate_pct = pass_rate_local(auth_epoch_tbl.rms_m, tau_st_m);
    fprintf('  tau_ST %.3f m (%gth percentile), authorized pass rate %.2f%% (%d trials)\n', ...
        tau_st_m, tau_quantile_pct, auth_pass_rate_pct, height(auth_epoch_tbl));

    threshold_tbl = table(tau_st_m, tau_quantile_pct, ...
        height(auth_epoch_tbl), auth_pass_rate_pct, ...
        median(auth_epoch_tbl.rms_m, 'omitnan'), ...
        min(auth_epoch_tbl.rms_m), max(auth_epoch_tbl.rms_m), ...
        'VariableNames', {'tau_st_m', 'tau_quantile_pct', ...
        'authorized_trial_count', 'authorized_pass_rate_pct', ...
        'authorized_median_rms_m', 'authorized_min_rms_m', ...
        'authorized_max_rms_m'});
    writetable(threshold_tbl, threshold_csv);
    writetable(auth_epoch_tbl, authorized_epoch_csv);

    fprintf('\nComputing spatial authentication pass-rate curve...\n');
    spatial_tbl = compute_spatial_pass_rate_local(nav_data, obs_data, ...
        lat0_deg, lon0_deg, height0_m, common_opts, min_epoch_rows, ...
        spatial_offsets_m, spatial_direction_deg, tau_st_m);
    writetable(spatial_tbl, spatial_csv);

    fprintf('\nComputing temporal authentication pass-rate curve...\n');
    temporal_tbl = compute_temporal_pass_rate_local(nav_data, obs_data, ...
        lat0_deg, lon0_deg, height0_m, common_opts, min_epoch_rows, ...
        time_offsets_s, tau_st_m);
    writetable(temporal_tbl, temporal_csv);
end

[spatial_pass_boundary_m, spatial_reject_start_m] = ...
    derive_spatial_pass_window_local(spatial_tbl);
[temporal_pass_boundary_s, temporal_reject_start_s] = ...
    derive_temporal_pass_window_local(temporal_tbl);

%% Figure 1: same time, different locations.
spatial_fig = figure('Name', 'Spatial authentication pass rate', ...
    'Color', 'w', 'Position', [120, 120, 1280, 860]);
plot(spatial_tbl.spatial_offset_m, spatial_tbl.pass_rate_pct, '-o', ...
    'Color', paper_line_color, ...
    'MarkerFaceColor', paper_marker_color, ...
    'MarkerEdgeColor', paper_line_color, ...
    'LineWidth', paper_curve_line_width, ...
    'MarkerSize', paper_marker_size);
hold on;
if isfinite(spatial_pass_boundary_m)
    plot([spatial_pass_boundary_m spatial_pass_boundary_m], [0 105], '--', ...
        'Color', paper_threshold_color, 'LineWidth', 2.4);
end
hold off;
grid on;
box on;
spatial_axis_min_m = min(spatial_tbl.spatial_offset_m);
spatial_axis_max_m = max(spatial_tbl.spatial_offset_m);
xlim([spatial_axis_min_m - 5, spatial_axis_max_m + 5]);
ylim([0, 105]);
xticks(0:20:spatial_axis_max_m);
yticks(0:20:100);
xlabel('Spatial offset (m)');
ylabel('Authentication pass rate (%)');
ax = gca;
apply_paper_axis_style_local(ax, paper_font_name, paper_tick_font_size, ...
    paper_label_font_size, paper_axis_line_width);
add_threshold_note_local(ax, build_spatial_threshold_note_local( ...
    tau_st_m, spatial_pass_boundary_m, spatial_reject_start_m), ...
    paper_font_name, paper_annotation_font_size);
set(ax, 'Units', 'normalized', 'Position', [0.30, 0.24, 0.64, 0.64]);

spatial_png = fullfile(out_png_dir, 'figure_spatial_auth_pass_rate.png');
spatial_pdf = fullfile(out_pdf_dir, 'figure_spatial_auth_pass_rate.pdf');
spatial_fig_path = fullfile(out_fig_dir, 'figure_spatial_auth_pass_rate.fig');
save_paper_figure_local(spatial_fig, spatial_png, spatial_pdf, spatial_fig_path);

%% Figure 2: same location, different times.
temporal_fig = figure('Name', 'Temporal authentication pass rate', ...
    'Color', 'w', 'Position', [140, 140, 1280, 860]);
plot(temporal_tbl.time_offset_s, temporal_tbl.pass_rate_pct, '-o', ...
    'Color', paper_line_color, ...
    'MarkerFaceColor', paper_marker_color, ...
    'MarkerEdgeColor', paper_line_color, ...
    'LineWidth', paper_curve_line_width, ...
    'MarkerSize', paper_marker_size);
hold on;
if isfinite(temporal_reject_start_s) && temporal_reject_start_s > 0
    plot([-temporal_reject_start_s -temporal_reject_start_s], [0 105], '--', ...
        'Color', paper_threshold_color, 'LineWidth', 2.4);
    plot([temporal_reject_start_s temporal_reject_start_s], [0 105], '--', ...
        'Color', paper_threshold_color, 'LineWidth', 2.4);
elseif isfinite(temporal_pass_boundary_s)
    plot([temporal_pass_boundary_s temporal_pass_boundary_s], [0 105], '--', ...
        'Color', paper_threshold_color, 'LineWidth', 2.4);
end
hold off;
grid on;
box on;
time_axis_min_s = min(temporal_tbl.time_offset_s);
time_axis_max_s = max(temporal_tbl.time_offset_s);
xlim([time_axis_min_s - 0.2, time_axis_max_s + 0.2]);
ylim([0, 105]);
xticks(ceil(time_axis_min_s):1:floor(time_axis_max_s));
yticks(0:20:100);
xlabel('Time offset (s)');
ylabel('Authentication pass rate (%)');
ax = gca;
apply_paper_axis_style_local(ax, paper_font_name, paper_tick_font_size, ...
    paper_label_font_size, paper_axis_line_width);
add_threshold_note_local(ax, build_temporal_threshold_note_local( ...
    tau_st_m, temporal_pass_boundary_s, temporal_reject_start_s), ...
    paper_font_name, paper_annotation_font_size);
set(ax, 'Units', 'normalized', 'Position', [0.30, 0.24, 0.64, 0.64]);

temporal_png = fullfile(out_png_dir, 'figure_temporal_auth_pass_rate.png');
temporal_pdf = fullfile(out_pdf_dir, 'figure_temporal_auth_pass_rate.pdf');
temporal_fig_path = fullfile(out_fig_dir, 'figure_temporal_auth_pass_rate.fig');
save_paper_figure_local(temporal_fig, temporal_png, temporal_pdf, temporal_fig_path);

fprintf('\nSaved authentication pass-rate figures:\n');
fprintf('  %s\n', spatial_png);
fprintf('  %s\n', spatial_pdf);
fprintf('  %s\n', temporal_png);
fprintf('  %s\n', temporal_pdf);

%% ===================== local helpers =====================
function [pass_boundary_m, reject_start_m] = derive_spatial_pass_window_local(spatial_tbl)
offsets = spatial_tbl.spatial_offset_m;
rates = spatial_tbl.pass_rate_pct;
valid = isfinite(offsets) & isfinite(rates);
offsets = offsets(valid);
rates = rates(valid);

pass_offsets = offsets(rates > 0);
if isempty(pass_offsets)
    pass_boundary_m = NaN;
else
    pass_boundary_m = max(pass_offsets);
end

reject_offsets = offsets(offsets > pass_boundary_m & rates == 0);
if isempty(reject_offsets)
    reject_start_m = NaN;
else
    reject_start_m = min(reject_offsets);
end
end

function [pass_boundary_s, reject_start_s] = derive_temporal_pass_window_local(temporal_tbl)
offsets = temporal_tbl.time_offset_s;
rates = temporal_tbl.pass_rate_pct;
valid = isfinite(offsets) & isfinite(rates);
offsets = offsets(valid);
rates = rates(valid);

pass_offsets = abs(offsets(rates > 0));
if isempty(pass_offsets)
    pass_boundary_s = NaN;
else
    pass_boundary_s = max(pass_offsets);
end

reject_offsets = abs(offsets(abs(offsets) > pass_boundary_s & rates == 0));
if isempty(reject_offsets)
    reject_start_s = NaN;
else
    reject_start_s = min(reject_offsets);
end
end

function note = build_spatial_threshold_note_local(tau_st_m, pass_boundary_m, ...
    reject_start_m)
if isfinite(pass_boundary_m)
    pass_line = sprintf('Nonzero pass: d \\leq %.0f m', pass_boundary_m);
else
    pass_line = 'Nonzero pass: none';
end

if isfinite(reject_start_m)
    reject_line = sprintf('0%% pass: d \\geq %.0f m', reject_start_m);
else
    reject_line = '0% pass boundary: not reached';
end

note = sprintf('\\tau_{ST} = %.2f m\n%s\n%s', ...
    tau_st_m, pass_line, reject_line);
end

function note = build_temporal_threshold_note_local(tau_st_m, pass_boundary_s, ...
    reject_start_s)
if isfinite(pass_boundary_s) && pass_boundary_s == 0
    pass_line = 'Nonzero pass: \Delta t = 0 s';
elseif isfinite(pass_boundary_s)
    pass_line = sprintf('Nonzero pass: |\\Delta t| <= %.1f s', pass_boundary_s);
else
    pass_line = 'Nonzero pass: none';
end

if isfinite(reject_start_s)
    reject_line = sprintf('0%% pass: |\\Delta t| \\geq %.1f s', reject_start_s);
elseif isfinite(pass_boundary_s)
    reject_line = '0% pass boundary: not reached';
else
    reject_line = '0% pass boundary: not reached';
end

note = sprintf('\\tau_{ST} = %.2f m\n%s\n%s', ...
    tau_st_m, pass_line, reject_line);
end

function add_threshold_note_local(ax, note, font_name, font_size)
text(ax, 0.97, 0.90, note, ...
    'Units', 'normalized', ...
    'HorizontalAlignment', 'right', ...
    'VerticalAlignment', 'top', ...
    'FontName', font_name, ...
    'FontSize', font_size, ...
    'Interpreter', 'tex', ...
    'BackgroundColor', 'white', ...
    'EdgeColor', [0.25 0.25 0.25], ...
    'Color', [0.05 0.05 0.05], ...
    'Margin', 8);
end

function spatial_tbl = compute_spatial_pass_rate_local(nav_data, obs_data, ...
    lat0_deg, lon0_deg, height0_m, common_opts, min_epoch_rows, ...
    offsets_m, direction_deg, tau_st_m)
n = numel(offsets_m);
pass_rate_pct = NaN(n, 1);
trial_count = zeros(n, 1);
pass_count = zeros(n, 1);
median_rms_m = NaN(n, 1);
p95_rms_m = NaN(n, 1);

for oi = 1:n
    offset_m = offsets_m(oi);
    if offset_m == 0
        dirs = 0;
    else
        dirs = direction_deg(:).';
    end

    all_rms = zeros(0, 1);
    for di = 1:numel(dirs)
        east_m = offset_m * cosd(dirs(di));
        north_m = offset_m * sind(dirs(di));
        [lat_deg, lon_deg, height_m] = offset_geodetic_local( ...
            lat0_deg, lon0_deg, height0_m, east_m, north_m, 0);

        epoch_tbl = evaluate_candidate_epoch_rms_local(nav_data, obs_data, ...
            lat_deg, lon_deg, height_m, common_opts, 0, min_epoch_rows);
        all_rms = [all_rms; epoch_tbl.rms_m]; %#ok<AGROW>
    end

    [pass_rate_pct(oi), pass_count(oi), trial_count(oi)] = ...
        pass_rate_local(all_rms, tau_st_m);
    median_rms_m(oi) = median(all_rms, 'omitnan');
    p95_rms_m(oi) = percentile_sorted_local(all_rms, 95);
    fprintf('  spatial offset %6.1f m -> pass %.2f%% (%d/%d)\n', ...
        offset_m, pass_rate_pct(oi), pass_count(oi), trial_count(oi));
end

spatial_tbl = table(offsets_m(:), pass_rate_pct, pass_count, trial_count, ...
    median_rms_m, p95_rms_m, ...
    'VariableNames', {'spatial_offset_m', 'pass_rate_pct', 'pass_count', ...
    'trial_count', 'median_rms_m', 'p95_rms_m'});
end

function temporal_tbl = compute_temporal_pass_rate_local(nav_data, obs_data, ...
    lat0_deg, lon0_deg, height0_m, common_opts, min_epoch_rows, ...
    time_offsets_s, tau_st_m)
n = numel(time_offsets_s);
pass_rate_pct = NaN(n, 1);
trial_count = zeros(n, 1);
pass_count = zeros(n, 1);
median_rms_m = NaN(n, 1);
p95_rms_m = NaN(n, 1);

for ti = 1:n
    offset_s = time_offsets_s(ti);
    epoch_tbl = evaluate_candidate_epoch_rms_local(nav_data, obs_data, ...
        lat0_deg, lon0_deg, height0_m, common_opts, offset_s, min_epoch_rows);

    [pass_rate_pct(ti), pass_count(ti), trial_count(ti)] = ...
        pass_rate_local(epoch_tbl.rms_m, tau_st_m);
    median_rms_m(ti) = median(epoch_tbl.rms_m, 'omitnan');
    p95_rms_m(ti) = percentile_sorted_local(epoch_tbl.rms_m, 95);
    fprintf('  time offset %+6.1f s -> pass %.2f%% (%d/%d)\n', ...
        offset_s, pass_rate_pct(ti), pass_count(ti), trial_count(ti));
end

temporal_tbl = table(time_offsets_s(:), pass_rate_pct, pass_count, ...
    trial_count, median_rms_m, p95_rms_m, ...
    'VariableNames', {'time_offset_s', 'pass_rate_pct', 'pass_count', ...
    'trial_count', 'median_rms_m', 'p95_rms_m'});
end

function epoch_tbl = evaluate_candidate_epoch_rms_local(nav_data, obs_data, ...
    lat_deg, lon_deg, height_m, common_opts, theory_time_shift_seconds, ...
    min_epoch_rows)
[residual_tbl, ~] = calculate_nav_obs_range_residuals( ...
    nav_data, obs_data, lat_deg, lon_deg, height_m, ...
    common_opts{:}, ...
    'theory_time_shift_seconds', theory_time_shift_seconds);

if isempty(residual_tbl) || height(residual_tbl) == 0
    epoch_tbl = empty_epoch_rms_table_local();
    return;
end

ok_tbl = residual_tbl(string(residual_tbl.status) == "ok", :);
if isempty(ok_tbl) || height(ok_tbl) == 0
    epoch_tbl = empty_epoch_rms_table_local();
    return;
end

epoch_ids = unique(ok_tbl.epoch_idx, 'stable');
n_epoch = numel(epoch_ids);
epoch_idx = zeros(n_epoch, 1);
time = NaT(n_epoch, 1, 'TimeZone', 'UTC');
row_count = zeros(n_epoch, 1);
satellite_count = zeros(n_epoch, 1);
rms_m = NaN(n_epoch, 1);
std_m = NaN(n_epoch, 1);

for ei = 1:n_epoch
    mask = ok_tbl.epoch_idx == epoch_ids(ei);
    vals = ok_tbl.residual_obs_minus_gt_m(mask);
    vals = vals(isfinite(vals));

    epoch_idx(ei) = epoch_ids(ei);
    first_row = find(mask, 1, 'first');
    time(ei) = ok_tbl.time(first_row);
    row_count(ei) = numel(vals);
    satellite_count(ei) = numel(unique(ok_tbl.sat_id(mask)));

    if row_count(ei) > 0
        rms_m(ei) = sqrt(mean(vals.^2, 'omitnan'));
        std_m(ei) = std(vals, 0, 'omitnan');
    end
end

epoch_tbl = table(epoch_idx, time, row_count, satellite_count, rms_m, std_m, ...
    'VariableNames', {'epoch_idx', 'time', 'row_count', ...
    'satellite_count', 'rms_m', 'std_m'});
epoch_tbl = epoch_tbl(epoch_tbl.row_count >= min_epoch_rows & ...
    isfinite(epoch_tbl.rms_m), :);
end

function epoch_tbl = empty_epoch_rms_table_local()
epoch_idx = zeros(0, 1);
time = NaT(0, 1, 'TimeZone', 'UTC');
row_count = zeros(0, 1);
satellite_count = zeros(0, 1);
rms_m = zeros(0, 1);
std_m = zeros(0, 1);
epoch_tbl = table(epoch_idx, time, row_count, satellite_count, rms_m, std_m, ...
    'VariableNames', {'epoch_idx', 'time', 'row_count', ...
    'satellite_count', 'rms_m', 'std_m'});
end

function [rate_pct, pass_count, trial_count] = pass_rate_local(values, tau_st_m)
values = values(:);
values = values(isfinite(values));
trial_count = numel(values);
if trial_count == 0
    rate_pct = NaN;
    pass_count = 0;
    return;
end

pass_count = sum(values <= tau_st_m);
rate_pct = 100 * pass_count / trial_count;
end

function value = percentile_sorted_local(y, pct)
y = sort(y(:));
y = y(isfinite(y));
if isempty(y)
    value = NaN;
    return;
end

n = numel(y);
if n == 1
    value = y(1);
    return;
end

rank = 1 + (n - 1) * pct / 100;
lo = floor(rank);
hi = ceil(rank);
if lo == hi
    value = y(lo);
else
    frac = rank - lo;
    value = y(lo) * (1 - frac) + y(hi) * frac;
end
end

function [lat_deg, lon_deg, height_m] = offset_geodetic_local( ...
    lat0_deg, lon0_deg, height0_m, east_m, north_m, up_m)
lat0 = deg2rad(lat0_deg);
lon0 = deg2rad(lon0_deg);
base_ecef = geodetic_to_ecef_local(lat0_deg, lon0_deg, height0_m);

east_axis = [-sin(lon0); cos(lon0); 0];
north_axis = [-sin(lat0) * cos(lon0); -sin(lat0) * sin(lon0); cos(lat0)];
up_axis = [cos(lat0) * cos(lon0); cos(lat0) * sin(lon0); sin(lat0)];

ecef = base_ecef + [east_axis, north_axis, up_axis] * [east_m; north_m; up_m];
[lat_deg, lon_deg, height_m] = ecef_to_geodetic_local(ecef);
end

function ecef = geodetic_to_ecef_local(lat_deg, lon_deg, height_m)
a = 6378137.0;
f = 1 / 298.257223563;
e2 = f * (2 - f);

lat = deg2rad(lat_deg);
lon = deg2rad(lon_deg);
sin_lat = sin(lat);
cos_lat = cos(lat);
N = a / sqrt(1 - e2 * sin_lat^2);

x = (N + height_m) * cos_lat * cos(lon);
y = (N + height_m) * cos_lat * sin(lon);
z = (N * (1 - e2) + height_m) * sin_lat;
ecef = [x; y; z];
end

function [lat_deg, lon_deg, height_m] = ecef_to_geodetic_local(ecef)
a = 6378137.0;
f = 1 / 298.257223563;
e2 = f * (2 - f);

x = ecef(1);
y = ecef(2);
z = ecef(3);
lon = atan2(y, x);
p = hypot(x, y);
lat = atan2(z, p * (1 - e2));

for i = 1:8
    sin_lat = sin(lat);
    N = a / sqrt(1 - e2 * sin_lat^2);
    height_m = p / cos(lat) - N;
    lat = atan2(z, p * (1 - e2 * N / (N + height_m)));
end

sin_lat = sin(lat);
N = a / sqrt(1 - e2 * sin_lat^2);
height_m = p / cos(lat) - N;
lat_deg = rad2deg(lat);
lon_deg = rad2deg(lon);
end

function apply_paper_axis_style_local(ax, font_name, tick_font_size, ...
    label_font_size, axis_line_width)
set(ax, ...
    'FontName', font_name, ...
    'FontSize', tick_font_size, ...
    'LineWidth', axis_line_width, ...
    'TickDir', 'in', ...
    'Layer', 'top', ...
    'Box', 'on');

if isprop(ax, 'XAxis') && isgraphics(ax.XAxis)
    ax.XAxis.FontName = font_name;
    ax.XAxis.FontSize = tick_font_size;
    if isgraphics(ax.XAxis.Label)
        ax.XAxis.Label.FontName = font_name;
        ax.XAxis.Label.FontSize = label_font_size;
    end
end

if isprop(ax, 'YAxis') && ~isempty(ax.YAxis)
    for k = 1:numel(ax.YAxis)
        if isgraphics(ax.YAxis(k))
            ax.YAxis(k).FontName = font_name;
            ax.YAxis(k).FontSize = tick_font_size;
            if isgraphics(ax.YAxis(k).Label)
                ax.YAxis(k).Label.FontName = font_name;
                ax.YAxis(k).Label.FontSize = label_font_size;
            end
        end
    end
end
end

function organize_existing_outputs_local(out_dir, png_dir, pdf_dir, fig_dir, csv_dir)
move_top_level_outputs_local(out_dir, '*.png', png_dir);
move_top_level_outputs_local(out_dir, '*.pdf', pdf_dir);
move_top_level_outputs_local(out_dir, '*.fig', fig_dir);
move_top_level_outputs_local(out_dir, '*.csv', csv_dir);
end

function move_top_level_outputs_local(out_dir, pattern, target_dir)
files = dir(fullfile(out_dir, pattern));
for k = 1:numel(files)
    if files(k).isdir
        continue;
    end

    src = fullfile(files(k).folder, files(k).name);
    dst = fullfile(target_dir, files(k).name);
    if strcmp(src, dst)
        continue;
    end

    [ok, msg] = movefile(src, dst, 'f');
    if ~ok
        warning('Could not move output file %s to %s: %s', src, target_dir, msg);
    end
end
end

function save_paper_figure_local(fig, png_path, pdf_path, fig_path)
drawnow;
exportgraphics(fig, png_path, 'Resolution', 300, ...
    'BackgroundColor', 'white');
exportgraphics(fig, pdf_path, 'ContentType', 'vector', ...
    'BackgroundColor', 'white');
savefig(fig, fig_path);
end

function repo_dir = resolve_repo_root_local(start_dir)
repo_dir = start_dir;
while true
    if isfolder(fullfile(repo_dir, 'data')) && ...
            isfolder(fullfile(repo_dir, 'nav_parse')) && ...
            isfolder(fullfile(repo_dir, 'obs_parse')) && ...
            isfolder(fullfile(repo_dir, 'calculate_clock_bias_and_positon'))
        return;
    end

    parent_dir = fileparts(repo_dir);
    if strcmp(parent_dir, repo_dir) || isempty(parent_dir)
        error('Cannot locate SatLock repo root from %s.', start_dir);
    end
    repo_dir = parent_dir;
end
end
