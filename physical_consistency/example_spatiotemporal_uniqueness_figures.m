% Example: draw range-domain spatiotemporal uniqueness figures.
%
% Figure 1: same time, different candidate locations.
% Figure 2: same location, different candidate times.
%
% The metric is the multi-satellite residual RMS after estimating receiver
% clock bias per epoch and per GNSS system:
%   e_i = r_obs_i - rho_i_GT(candidate_time, candidate_position)

clear;
clc;
close all;

repo_dir = resolve_repo_root_local(fileparts(mfilename('fullpath')));
addpath(genpath(repo_dir));

%% ===================== Edit here =====================
nav_file = fullfile(repo_dir, 'data', '1_8', '2026_1_8.nav');
obs_file = fullfile(repo_dir, 'data', '1_8', 'A_1_8_1.obs');

start_time = datetime(2026, 1, 8, 5, 59, 39.116, 'TimeZone', 'UTC');

% Reference spatiotemporal point.
lat0_deg = 35.77204196;
lon0_deg = 120.02477210;
height0_m = 92.889;

systems = {'G', 'C', 'E', 'R', 'J'};
max_epochs = 5;
obs_window_seconds = 0.5;

% Spatial heatmap grid around the reference location.
space_offsets_m = -300:60:300;

% Temporal offsets for same-location uniqueness.
time_offsets_s = [-3600, -1800, -900, -600, -300, -120, -60, -30, -10, ...
                  0, 10, 30, 60, 120, 300, 600, 900, 1800, 3600];

% Satellite-wise geometric range-difference distribution for spatial offsets.
range_difference_offsets_m = 0:5:500;

out_dir = fullfile(repo_dir, 'physical_consistency', 'outputs', ...
    'spatiotemporal_uniqueness');
out_pdf_dir = fullfile(out_dir, 'pdf');
out_png_dir = fullfile(out_dir, 'png');
out_fig_dir = fullfile(out_dir, 'fig');
out_csv_dir = fullfile(out_dir, 'csv');

paper_font_name = 'Arial';
paper_tick_font_size = 30;
paper_label_font_size = 30;
paper_legend_font_size = 28;
paper_axis_line_width = 1.2;
paper_curve_line_width = 3.0;
paper_marker_size = 10;
%% =====================================================

output_dirs = {out_dir, out_pdf_dir, out_png_dir, out_fig_dir, out_csv_dir};
for di = 1:numel(output_dirs)
    if ~exist(output_dirs{di}, 'dir')
        mkdir(output_dirs{di});
    end
end
organize_existing_outputs_local(out_dir, out_png_dir, out_pdf_dir, ...
    out_fig_dir, out_csv_dir);

fprintf('Parsing NAV/OBS once...\n');
nav_data = parse_rinex_nav_multi_gnss(nav_file);
obs_data = parse_rinex_obs(obs_file);

common_opts = { ...
    'systems', systems, ...
    'elevation_mask_deg', 0, ...
    'min_clock_sats', 4, ...
    'receiver_clock_estimator', 'mean', ...
    'start_time', start_time, ...
    'duration_seconds', obs_window_seconds, ...
    'max_epochs', max_epochs};

fprintf('\nComputing spatial uniqueness heatmap (%d x %d candidates)...\n', ...
    numel(space_offsets_m), numel(space_offsets_m));
spatial_rms_m = NaN(numel(space_offsets_m), numel(space_offsets_m));
spatial_std_m = NaN(size(spatial_rms_m));
spatial_sat_count = NaN(size(spatial_rms_m));
spatial_row_count = NaN(size(spatial_rms_m));

for ni = 1:numel(space_offsets_m)
    north_m = space_offsets_m(ni);
    for ei = 1:numel(space_offsets_m)
        east_m = space_offsets_m(ei);
        [lat_deg, lon_deg, height_m] = offset_geodetic_local( ...
            lat0_deg, lon0_deg, height0_m, east_m, north_m, 0);

        metric = evaluate_candidate_local(nav_data, obs_data, ...
            lat_deg, lon_deg, height_m, common_opts, 0);

        spatial_rms_m(ni, ei) = metric.rms_m;
        spatial_std_m(ni, ei) = metric.std_m;
        spatial_sat_count(ni, ei) = metric.sat_count;
        spatial_row_count(ni, ei) = metric.row_count;
    end
    fprintf('  spatial row %2d/%2d done\n', ni, numel(space_offsets_m));
end

fprintf('\nComputing temporal uniqueness curve (%d candidates)...\n', ...
    numel(time_offsets_s));
time_rms_m = NaN(numel(time_offsets_s), 1);
time_std_m = NaN(size(time_rms_m));
time_sat_count = NaN(size(time_rms_m));
time_row_count = NaN(size(time_rms_m));

for ti = 1:numel(time_offsets_s)
    metric = evaluate_candidate_local(nav_data, obs_data, ...
        lat0_deg, lon0_deg, height0_m, common_opts, time_offsets_s(ti));
    time_rms_m(ti) = metric.rms_m;
    time_std_m(ti) = metric.std_m;
    time_sat_count(ti) = metric.sat_count;
    time_row_count(ti) = metric.row_count;
    fprintf('  time offset %+7.1f s -> RMS %.3f m\n', ...
        time_offsets_s(ti), metric.rms_m);
end

fprintf('\nComputing satellite-wise geometric range-difference distribution...\n');
range_diff_tbl = compute_range_difference_distribution_local( ...
    nav_data, start_time, lat0_deg, lon0_deg, height0_m, systems, ...
    range_difference_offsets_m);
fprintf('  distribution rows: %d | satellites: %d\n', height(range_diff_tbl), ...
    numel(unique(range_diff_tbl.sat_id)));

%% Save numerical tables.
[east_grid, north_grid] = meshgrid(space_offsets_m, space_offsets_m);
spatial_tbl = table(east_grid(:), north_grid(:), spatial_rms_m(:), ...
    spatial_std_m(:), spatial_sat_count(:), spatial_row_count(:), ...
    'VariableNames', {'east_offset_m', 'north_offset_m', 'residual_rms_m', ...
    'residual_std_m', 'satellite_count', 'row_count'});
writetable(spatial_tbl, fullfile(out_csv_dir, 'spatial_uniqueness_heatmap.csv'));

time_tbl = table(time_offsets_s(:), time_rms_m, time_std_m, ...
    time_sat_count, time_row_count, ...
    'VariableNames', {'time_offset_s', 'residual_rms_m', ...
    'residual_std_m', 'satellite_count', 'row_count'});
writetable(time_tbl, fullfile(out_csv_dir, 'temporal_uniqueness_curve.csv'));

writetable(range_diff_tbl, fullfile(out_csv_dir, ...
    'satellite_range_difference_distribution.csv'));

%% Figure 1a: spatial uniqueness heatmap.
space_heatmap_fig = figure('Name', 'Spatial uniqueness heatmap', ...
    'Color', 'w', 'Position', [80, 80, 1180, 920]);
imagesc(space_offsets_m, space_offsets_m, spatial_rms_m);
ax = gca;
set(ax, 'YDir', 'normal');
axis equal tight;
hold on;
plot(0, 0, 'wx', 'LineWidth', 3.2, 'MarkerSize', 18);
contour(space_offsets_m, space_offsets_m, spatial_rms_m, ...
    'LineColor', [1 1 1] * 0.95, 'LineWidth', 1.15);
hold off;
grid on;
box on;
apply_colormap_local();
cb = colorbar;
cb.Label.String = 'Residual RMS (m)';
xlabel('East offset (m)');
ylabel('North offset (m)');
apply_paper_axis_style_local(ax, paper_font_name, paper_tick_font_size, ...
    paper_label_font_size, paper_axis_line_width);
apply_paper_colorbar_style_local(cb, paper_font_name, paper_tick_font_size, ...
    paper_label_font_size);
set(ax, 'Units', 'normalized', 'Position', [0.25, 0.20, 0.50, 0.70]);
set(cb, 'Units', 'normalized', 'Position', [0.84, 0.20, 0.035, 0.70]);
space_heatmap_png = fullfile(out_png_dir, 'figure_space_uniqueness_heatmap.png');
space_heatmap_pdf = fullfile(out_pdf_dir, 'figure_space_uniqueness_heatmap.pdf');
space_heatmap_fig_path = fullfile(out_fig_dir, 'figure_space_uniqueness_heatmap.fig');
save_paper_figure_local(space_heatmap_fig, space_heatmap_png, ...
    space_heatmap_pdf, space_heatmap_fig_path);

%% Figure 1b: spatial uniqueness cross-section.
space_cross_fig = figure('Name', 'Spatial uniqueness cross sections', ...
    'Color', 'w', 'Position', [100, 100, 1280, 860]);
zero_idx = find(space_offsets_m == 0, 1);
plot(space_offsets_m, spatial_rms_m(zero_idx, :), '-o', ...
    'LineWidth', paper_curve_line_width, 'MarkerSize', paper_marker_size, ...
    'DisplayName', 'East-West');
hold on;
plot(space_offsets_m, spatial_rms_m(:, zero_idx), '-s', ...
    'LineWidth', paper_curve_line_width, 'MarkerSize', paper_marker_size, ...
    'DisplayName', 'North-South');
hold off;
xlim([min(space_offsets_m) - 20, max(space_offsets_m) + 20]);
space_cross_ymax = max([spatial_rms_m(zero_idx, :), ...
    spatial_rms_m(:, zero_idx)'], [], 'omitnan');
ylim([0, space_cross_ymax * 1.10]);
grid on;
box on;
xlabel('Offset (m)');
ylabel('Residual RMS (m)');
lgd = legend('Location', 'northwest');
ax = gca;
apply_paper_axis_style_local(ax, paper_font_name, paper_tick_font_size, ...
    paper_label_font_size, paper_axis_line_width);
apply_paper_legend_style_local(lgd, paper_font_name, paper_legend_font_size);
set(ax, 'Units', 'normalized', 'Position', [0.33, 0.24, 0.60, 0.64]);
space_cross_png = fullfile(out_png_dir, 'figure_space_uniqueness_cross_section.png');
space_cross_pdf = fullfile(out_pdf_dir, 'figure_space_uniqueness_cross_section.pdf');
space_cross_fig_path = fullfile(out_fig_dir, 'figure_space_uniqueness_cross_section.fig');
save_paper_figure_local(space_cross_fig, space_cross_png, ...
    space_cross_pdf, space_cross_fig_path);

%% Figure 2a: temporal uniqueness full range.
time_full_fig = figure('Name', 'Temporal uniqueness full range', ...
    'Color', 'w', 'Position', [120, 120, 1280, 860]);
semilogy(time_offsets_s / 60, time_rms_m, '-o', ...
    'LineWidth', paper_curve_line_width, 'MarkerSize', paper_marker_size, ...
    'MarkerFaceColor', [0.1 0.45 0.85]);
xlim([min(time_offsets_s) / 60 - 3, max(time_offsets_s) / 60 + 3]);
ylim([1, max(time_rms_m, [], 'omitnan') * 1.60]);
grid off;
box on;
xlabel('Time offset (min)');
ylabel('Residual RMS (m)');
ax = gca;
apply_paper_axis_style_local(ax, paper_font_name, paper_tick_font_size, ...
    paper_label_font_size, paper_axis_line_width);
set(ax, 'Units', 'normalized', 'Position', [0.33, 0.24, 0.60, 0.64]);
time_full_png = fullfile(out_png_dir, 'figure_time_uniqueness_full_range.png');
time_full_pdf = fullfile(out_pdf_dir, 'figure_time_uniqueness_full_range.pdf');
time_full_fig_path = fullfile(out_fig_dir, 'figure_time_uniqueness_full_range.fig');
save_paper_figure_local(time_full_fig, time_full_png, ...
    time_full_pdf, time_full_fig_path);

%% Figure 2b: temporal uniqueness near-authorized-time zoom.
time_zoom_fig = figure('Name', 'Temporal uniqueness near-authorized-time zoom', ...
    'Color', 'w', 'Position', [140, 140, 1280, 860]);
zoom_mask = abs(time_offsets_s) <= 60;
time_zoom_line = plot(time_offsets_s(zoom_mask), time_rms_m(zoom_mask), '-o', ...
    'LineWidth', paper_curve_line_width, 'MarkerSize', paper_marker_size, ...
    'MarkerFaceColor', [0.1 0.45 0.85]);
xlim([min(time_offsets_s(zoom_mask)) - 5, max(time_offsets_s(zoom_mask)) + 5]);
zoom_ymax = max(time_rms_m(zoom_mask), [], 'omitnan');
ylim([0, zoom_ymax * 1.15]);
set(time_zoom_line, 'Clipping', 'off');
grid off;
box on;
xlabel('Time offset (s)');
ylabel('Residual RMS (m)');
ax = gca;
apply_paper_axis_style_local(ax, paper_font_name, paper_tick_font_size, ...
    paper_label_font_size, paper_axis_line_width);
set(ax, 'Units', 'normalized', 'Position', [0.33, 0.24, 0.60, 0.64]);
time_zoom_png = fullfile(out_png_dir, 'figure_time_uniqueness_near_time_zoom.png');
time_zoom_pdf = fullfile(out_pdf_dir, 'figure_time_uniqueness_near_time_zoom.pdf');
time_zoom_fig_path = fullfile(out_fig_dir, 'figure_time_uniqueness_near_time_zoom.fig');
save_paper_figure_local(time_zoom_fig, time_zoom_png, ...
    time_zoom_pdf, time_zoom_fig_path);

%% Figure 3: satellite-wise range-difference distribution variants.
ymax_diff = ceil(max(abs(range_diff_tbl.delta_range_m), [], 'omitnan') / 50) * 50;
range_summary_tbl = summarize_range_difference_distribution_local( ...
    range_diff_tbl, range_difference_offsets_m);
writetable(range_summary_tbl, fullfile(out_csv_dir, ...
    'satellite_range_difference_distribution_summary.csv'));

range_diff_png = fullfile(out_png_dir, ...
    'figure_satellite_range_difference_distribution.png');
range_diff_pdf = fullfile(out_pdf_dir, ...
    'figure_satellite_range_difference_distribution.pdf');
range_diff_fig_path = fullfile(out_fig_dir, ...
    'figure_satellite_range_difference_distribution.fig');
plot_range_quantile_band_local(range_summary_tbl, ymax_diff, ...
    range_diff_png, range_diff_pdf, range_diff_fig_path, ...
    paper_font_name, paper_tick_font_size, paper_label_font_size, ...
    paper_axis_line_width, paper_legend_font_size, false);

violin_png = fullfile(out_png_dir, ...
    'figure_satellite_range_difference_violin_envelope.png');
violin_pdf = fullfile(out_pdf_dir, ...
    'figure_satellite_range_difference_violin_envelope.pdf');
violin_fig_path = fullfile(out_fig_dir, ...
    'figure_satellite_range_difference_violin_envelope.fig');
plot_range_violin_envelope_local(range_diff_tbl, range_summary_tbl, ...
    100:100:500, ymax_diff, violin_png, violin_pdf, violin_fig_path, ...
    paper_font_name, paper_tick_font_size, paper_label_font_size, ...
    paper_axis_line_width, paper_legend_font_size);

fan_png = fullfile(out_png_dir, ...
    'figure_satellite_range_difference_fan_chart.png');
fan_pdf = fullfile(out_pdf_dir, ...
    'figure_satellite_range_difference_fan_chart.pdf');
fan_fig_path = fullfile(out_fig_dir, ...
    'figure_satellite_range_difference_fan_chart.fig');
plot_range_quantile_band_local(range_summary_tbl, ymax_diff, ...
    fan_png, fan_pdf, fan_fig_path, ...
    paper_font_name, paper_tick_font_size, paper_label_font_size, ...
    paper_axis_line_width, paper_legend_font_size, true);

fprintf('\nSaved figures:\n');
fprintf('  %s\n', space_heatmap_png);
fprintf('  %s\n', space_heatmap_pdf);
fprintf('  %s\n', space_cross_png);
fprintf('  %s\n', space_cross_pdf);
fprintf('  %s\n', time_full_png);
fprintf('  %s\n', time_full_pdf);
fprintf('  %s\n', time_zoom_png);
fprintf('  %s\n', time_zoom_pdf);
fprintf('  %s\n', range_diff_png);
fprintf('  %s\n', range_diff_pdf);
fprintf('  %s\n', violin_png);
fprintf('  %s\n', violin_pdf);
fprintf('  %s\n', fan_png);
fprintf('  %s\n', fan_pdf);
fprintf('\nAuthorized-point RMS: %.3f m\n', spatial_rms_m(zero_idx, zero_idx));

%% ===================== local helpers =====================
function metric = evaluate_candidate_local(nav_data, obs_data, lat_deg, lon_deg, ...
    height_m, common_opts, theory_time_shift_seconds)
[residual_tbl, ~] = calculate_nav_obs_range_residuals( ...
    nav_data, obs_data, lat_deg, lon_deg, height_m, ...
    common_opts{:}, ...
    'theory_time_shift_seconds', theory_time_shift_seconds);

ok_tbl = residual_tbl(residual_tbl.status == "ok", :);
residual = ok_tbl.residual_obs_minus_gt_m;

metric = struct();
metric.row_count = height(ok_tbl);
metric.sat_count = numel(unique(ok_tbl.sat_id));
metric.mean_m = mean(residual, 'omitnan');
metric.std_m = std(residual, 0, 'omitnan');
metric.rms_m = sqrt(mean(residual.^2, 'omitnan'));
end

function range_diff_tbl = compute_range_difference_distribution_local( ...
    nav_data, target_time, lat0_deg, lon0_deg, height0_m, systems, offsets_m)
ref_tbl = calculate_satellite_ranges_from_nav( ...
    nav_data, target_time, lat0_deg, lon0_deg, height0_m, ...
    'systems', systems, ...
    'elevation_mask_deg', 0, ...
    'visible_only', true);

ref_tbl = ref_tbl(isfinite(ref_tbl.range_m), :);

sat_id = strings(0, 1);
system = strings(0, 1);
prn = zeros(0, 1);
candidate_offset_m = zeros(0, 1);
candidate_east_m = zeros(0, 1);
candidate_north_m = zeros(0, 1);
candidate_lat_deg = zeros(0, 1);
candidate_lon_deg = zeros(0, 1);
candidate_height_m = zeros(0, 1);
azimuth_deg = zeros(0, 1);
elevation_deg = zeros(0, 1);
rho_authorized_m = zeros(0, 1);
rho_candidate_m = zeros(0, 1);
delta_range_m = zeros(0, 1);

for oi = 1:numel(offsets_m)
    offset_m = offsets_m(oi);
    east_m = offset_m;
    north_m = 0;
    [cand_lat_deg, cand_lon_deg, cand_height_m] = offset_geodetic_local( ...
        lat0_deg, lon0_deg, height0_m, east_m, north_m, 0);

    cand_tbl = calculate_satellite_ranges_from_nav( ...
        nav_data, target_time, cand_lat_deg, cand_lon_deg, cand_height_m, ...
        'systems', systems, ...
        'elevation_mask_deg', -90, ...
        'visible_only', false);

    [has_match, cand_idx] = ismember(ref_tbl.sat_id, cand_tbl.sat_id);
    ref_idx = find(has_match);
    cand_idx = cand_idx(has_match);

    valid = isfinite(ref_tbl.range_m(ref_idx)) & isfinite(cand_tbl.range_m(cand_idx));
    ref_idx = ref_idx(valid);
    cand_idx = cand_idx(valid);

    n = numel(ref_idx);
    if n == 0
        continue;
    end

    sat_id = [sat_id; string(ref_tbl.sat_id(ref_idx))]; %#ok<AGROW>
    system = [system; string(ref_tbl.system(ref_idx))]; %#ok<AGROW>
    prn = [prn; ref_tbl.prn(ref_idx)]; %#ok<AGROW>
    candidate_offset_m = [candidate_offset_m; repmat(offset_m, n, 1)]; %#ok<AGROW>
    candidate_east_m = [candidate_east_m; repmat(east_m, n, 1)]; %#ok<AGROW>
    candidate_north_m = [candidate_north_m; repmat(north_m, n, 1)]; %#ok<AGROW>
    candidate_lat_deg = [candidate_lat_deg; repmat(cand_lat_deg, n, 1)]; %#ok<AGROW>
    candidate_lon_deg = [candidate_lon_deg; repmat(cand_lon_deg, n, 1)]; %#ok<AGROW>
    candidate_height_m = [candidate_height_m; repmat(cand_height_m, n, 1)]; %#ok<AGROW>
    azimuth_deg = [azimuth_deg; ref_tbl.azimuth_deg(ref_idx)]; %#ok<AGROW>
    elevation_deg = [elevation_deg; ref_tbl.elevation_deg(ref_idx)]; %#ok<AGROW>
    rho_authorized_m = [rho_authorized_m; ref_tbl.range_m(ref_idx)]; %#ok<AGROW>
    rho_candidate_m = [rho_candidate_m; cand_tbl.range_m(cand_idx)]; %#ok<AGROW>
    delta_range_m = [delta_range_m; cand_tbl.range_m(cand_idx) - ...
        ref_tbl.range_m(ref_idx)]; %#ok<AGROW>
end

range_diff_tbl = table(candidate_offset_m, candidate_east_m, candidate_north_m, ...
    sat_id, system, prn, azimuth_deg, elevation_deg, ...
    candidate_lat_deg, candidate_lon_deg, candidate_height_m, ...
    rho_authorized_m, rho_candidate_m, delta_range_m, ...
    'VariableNames', {'candidate_offset_m', 'candidate_east_m', ...
    'candidate_north_m', 'sat_id', 'system', 'prn', 'azimuth_deg', ...
    'elevation_deg', 'candidate_lat_deg', 'candidate_lon_deg', ...
    'candidate_height_m', 'rho_authorized_m', 'rho_candidate_m', ...
    'delta_range_m'});
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

function summary_tbl = summarize_range_difference_distribution_local( ...
    range_diff_tbl, offsets_m)
n = numel(offsets_m);
p05 = NaN(n, 1);
p25 = NaN(n, 1);
p50 = NaN(n, 1);
p75 = NaN(n, 1);
p95 = NaN(n, 1);
satellite_count = zeros(n, 1);

for oi = 1:n
    mask = range_diff_tbl.candidate_offset_m == offsets_m(oi);
    y = range_diff_tbl.delta_range_m(mask);
    y = y(isfinite(y));
    satellite_count(oi) = numel(y);
    p05(oi) = percentile_sorted_local(y, 5);
    p25(oi) = percentile_sorted_local(y, 25);
    p50(oi) = percentile_sorted_local(y, 50);
    p75(oi) = percentile_sorted_local(y, 75);
    p95(oi) = percentile_sorted_local(y, 95);
end

summary_tbl = table(offsets_m(:), satellite_count, p05, p25, p50, p75, p95, ...
    'VariableNames', {'candidate_offset_m', 'satellite_count', ...
    'p05_delta_range_m', 'p25_delta_range_m', 'median_delta_range_m', ...
    'p75_delta_range_m', 'p95_delta_range_m'});
end

function plot_range_quantile_band_local(summary_tbl, ymax_diff, png_path, ...
    pdf_path, fig_path, font_name, tick_font_size, label_font_size, ...
    axis_line_width, legend_font_size, show_theory_bounds)
fig = figure('Name', 'Satellite-wise range-difference quantile band', ...
    'Color', 'w', 'Position', [160, 160, 1280, 860]);
ax = axes(fig);
hold(ax, 'on');

x = summary_tbl.candidate_offset_m;
x_limits = [min(x) - 35, max(x) + 35];
outer_color = [0.23 0.66 0.88];
inner_color = [0.19 0.40 0.67];
median_color = [0.02 0.08 0.18];
boundary_color = [0.48 0.48 0.48];

outer_band_h = draw_range_band_local(ax, x, summary_tbl.p05_delta_range_m, ...
    summary_tbl.p95_delta_range_m, outer_color, 0.18);
inner_band_h = draw_range_band_local(ax, x, summary_tbl.p25_delta_range_m, ...
    summary_tbl.p75_delta_range_m, inner_color, 0.30);

if show_theory_bounds
    theory_h = plot(ax, x, x, '--', 'Color', boundary_color, 'LineWidth', 1.8);
    plot(ax, x, -x, '--', 'Color', boundary_color, 'LineWidth', 1.8, ...
        'HandleVisibility', 'off');
end

p05_h = plot(ax, x, summary_tbl.p05_delta_range_m, '-', ...
    'Color', outer_color, 'LineWidth', 1.5);
plot(ax, x, summary_tbl.p95_delta_range_m, '-', ...
    'Color', outer_color, 'LineWidth', 1.5, 'HandleVisibility', 'off');
p25_h = plot(ax, x, summary_tbl.p25_delta_range_m, '-', ...
    'Color', inner_color, 'LineWidth', 2.1);
plot(ax, x, summary_tbl.p75_delta_range_m, '-', ...
    'Color', inner_color, 'LineWidth', 2.1, 'HandleVisibility', 'off');
median_h = plot(ax, x, summary_tbl.median_delta_range_m, '-', ...
    'Color', median_color, 'LineWidth', 3.4);
zero_h = plot(ax, x_limits, [0 0], '-', 'Color', [0.20 0.20 0.20], ...
    'LineWidth', 1.1);

hold(ax, 'off');
grid(ax, 'on');
box(ax, 'on');
xlim(ax, x_limits);
ylim(ax, [-ymax_diff - 20, ymax_diff + 20]);
yticks(ax, -ymax_diff:100:ymax_diff);
xlabel(ax, 'Spatial offset (m)');
ylabel(ax, '\Delta\rho_i (m)');
apply_paper_axis_style_local(ax, font_name, tick_font_size, ...
    label_font_size, axis_line_width);
if show_theory_bounds
    legend_handles = [outer_band_h, inner_band_h, median_h, theory_h, zero_h];
    legend_labels = {'5%-95% quantile band', '25%-75% quantile band', ...
        'Median (Q50)', 'Theoretical +/-d bound', 'Zero reference'};
else
    legend_handles = [outer_band_h, inner_band_h, median_h, zero_h];
    legend_labels = {'5%-95% quantile band', '25%-75% quantile band', ...
        'Median (Q50)', 'Zero reference'};
end
lgd = legend(ax, legend_handles, legend_labels, 'Location', 'northwest', ...
    'NumColumns', 2);
apply_paper_legend_style_local(lgd, font_name, legend_font_size - 4);
set(ax, 'Units', 'normalized', 'Position', [0.25, 0.22, 0.70, 0.68]);
save_paper_figure_local(fig, png_path, pdf_path, fig_path, 'tight');
end

function h = draw_range_band_local(ax, x, y_low, y_high, color_rgb, face_alpha)
h = gobjects(1);
x = x(:);
y_low = y_low(:);
y_high = y_high(:);
valid = isfinite(x) & isfinite(y_low) & isfinite(y_high);
if nnz(valid) < 2
    return;
end

x = x(valid);
y_low = y_low(valid);
y_high = y_high(valid);
h = patch(ax, [x; flipud(x)], [y_low; flipud(y_high)], color_rgb, ...
    'EdgeColor', 'none', 'FaceAlpha', face_alpha);
end

function plot_range_violin_envelope_local(range_diff_tbl, summary_tbl, ...
    violin_offsets_m, ymax_diff, png_path, pdf_path, fig_path, font_name, ...
    tick_font_size, label_font_size, axis_line_width, legend_font_size)
fig = figure('Name', 'Satellite-wise range-difference violin envelope', ...
    'Color', 'w', 'Position', [180, 180, 1280, 860]);
ax = axes(fig);
hold(ax, 'on');

x = summary_tbl.candidate_offset_m;
x_limits = [min(x) - 35, max(x) + 35];
envelope_h = plot(ax, x, summary_tbl.p05_delta_range_m, '-', ...
    'Color', [0.48 0.48 0.48], 'LineWidth', 1.8);
plot(ax, x, summary_tbl.p95_delta_range_m, '-', ...
    'Color', [0.48 0.48 0.48], 'LineWidth', 1.8, ...
    'HandleVisibility', 'off');

violin_h = gobjects(1);
for oi = 1:numel(violin_offsets_m)
    offset_m = violin_offsets_m(oi);
    mask = range_diff_tbl.candidate_offset_m == offset_m;
    y = range_diff_tbl.delta_range_m(mask);
    h = draw_violin_local(ax, offset_m, y, 20, [0.23 0.66 0.88]);
    if oi == 1
        violin_h = h;
    elseif isgraphics(h)
        h.HandleVisibility = 'off';
    end
end

median_h = plot(ax, x, summary_tbl.median_delta_range_m, '-', ...
    'Color', [0.02 0.08 0.18], 'LineWidth', 3.2);
zero_h = plot(ax, x_limits, [0 0], '-', 'Color', [0.20 0.20 0.20], ...
    'LineWidth', 1.1);

hold(ax, 'off');
grid(ax, 'on');
box(ax, 'on');
xlim(ax, x_limits);
ylim(ax, [-ymax_diff - 20, ymax_diff + 20]);
yticks(ax, -ymax_diff:100:ymax_diff);
xlabel(ax, 'Spatial offset (m)');
ylabel(ax, '\Delta\rho_i (m)');
apply_paper_axis_style_local(ax, font_name, tick_font_size, ...
    label_font_size, axis_line_width);
lgd = legend(ax, [violin_h, envelope_h, median_h, zero_h], ...
    {'Satellite distribution', '5%/95% envelope', ...
    'Median (Q50)', 'Zero reference'}, ...
    'Location', 'northwest', 'NumColumns', 2);
apply_paper_legend_style_local(lgd, font_name, legend_font_size - 6);
set(ax, 'Units', 'normalized', 'Position', [0.25, 0.22, 0.70, 0.68]);
save_paper_figure_local(fig, png_path, pdf_path, fig_path, 'tight');
end

function h = draw_violin_local(ax, x0, y, width_m, color_rgb)
h = gobjects(1);
y = y(:);
y = y(isfinite(y));
if numel(y) < 3
    return;
end

n_bins = 42;
y_edges = linspace(min(y), max(y), n_bins + 1);
if all(~isfinite(y_edges)) || y_edges(1) == y_edges(end)
    return;
end

y_centers = y_edges(1:end-1) + diff(y_edges) / 2;
counts = histcounts(y, y_edges);
counts = conv(counts(:), gaussian_kernel_local(1.2), 'same');
if max(counts) <= 0
    return;
end

half_width = width_m * counts(:) / max(counts);
h = patch(ax, [x0 - half_width; flipud(x0 + half_width)], ...
    [y_centers(:); flipud(y_centers(:))], color_rgb, ...
    'FaceAlpha', 0.35, 'EdgeColor', [0.25 0.25 0.25], ...
    'LineWidth', 1.2);
end

function kernel = gaussian_kernel_local(sigma)
sigma = max(double(sigma), eps);
radius = max(1, ceil(3 * sigma));
x = -radius:radius;
kernel = exp(-(x .^ 2) ./ (2 * sigma ^ 2));
kernel = kernel ./ sum(kernel);
end

function value = percentile_sorted_local(y_sorted, pct)
y_sorted = sort(y_sorted(:));
y_sorted = y_sorted(isfinite(y_sorted));
n = numel(y_sorted);
if n == 0
    value = NaN;
    return;
end
if n == 1
    value = y_sorted(1);
    return;
end

rank = 1 + (n - 1) * pct / 100;
lo = floor(rank);
hi = ceil(rank);
if lo == hi
    value = y_sorted(lo);
else
    frac = rank - lo;
    value = y_sorted(lo) * (1 - frac) + y_sorted(hi) * frac;
end
end

function apply_colormap_local()
try
    colormap(turbo);
catch
    colormap(parula);
end
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

function apply_paper_colorbar_style_local(cb, font_name, tick_font_size, ...
    label_font_size)
cb.FontName = font_name;
cb.FontSize = tick_font_size;
if isgraphics(cb.Label)
    cb.Label.FontName = font_name;
    cb.Label.FontSize = label_font_size;
end
end

function apply_paper_legend_style_local(lgd, font_name, font_size)
if isempty(lgd) || ~isgraphics(lgd)
    return;
end
set(lgd, ...
    'FontName', font_name, ...
    'FontSize', font_size, ...
    'Box', 'on');
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

function save_paper_figure_local(fig, png_path, pdf_path, fig_path, varargin)
drawnow;

export_mode = "tight";
if ~isempty(varargin)
    export_mode = string(varargin{1});
end

if export_mode == "tight"
    exportgraphics(fig, png_path, 'Resolution', 300, ...
        'BackgroundColor', 'white');
    exportgraphics(fig, pdf_path, 'ContentType', 'vector', ...
        'BackgroundColor', 'white');
    savefig(fig, fig_path);
    return;
end

pdf_margin_in = 0.25;
fig_pos = get(fig, 'Position');
screen_ppi = get(0, 'ScreenPixelsPerInch');
if isempty(screen_ppi) || ~isfinite(screen_ppi) || screen_ppi <= 0
    screen_ppi = 96;
end

paper_width_in = fig_pos(3) / screen_ppi;
paper_height_in = fig_pos(4) / screen_ppi;

set(fig, ...
    'InvertHardcopy', 'off', ...
    'PaperUnits', 'inches', ...
    'PaperSize', [paper_width_in + 2 * pdf_margin_in, ...
                  paper_height_in + 2 * pdf_margin_in], ...
    'PaperPosition', [pdf_margin_in, pdf_margin_in, ...
                      paper_width_in, paper_height_in], ...
    'PaperPositionMode', 'manual');

print(fig, png_path, '-dpng', '-r300');
print(fig, pdf_path, '-dpdf', '-painters');
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
