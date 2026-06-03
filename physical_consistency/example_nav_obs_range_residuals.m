% Example: align pure NAV geometric ranges and OBS corrected pseudoranges.
%
% Goal:
%   rho_gt_i = pure NAV geometric range from satellite i to a fixed point.
%   r_obs_i = raw pseudorange corrected by satellite clock and receiver clock.
%   residual_i = r_obs_i - rho_gt_i.
%
% This example uses one NAV/OBS pair and a short time window. Multi-GNSS is
% enabled by default. Receiver clock bias is estimated separately per system
% to avoid mixing inter-system biases into one clock term.

clear;
clc;
close all;

repo_dir = resolve_repo_root_local(fileparts(mfilename('fullpath')));
addpath(genpath(repo_dir));

%% ===================== Edit here =====================
nav_file = fullfile(repo_dir, 'data', '1_8', '2026_1_8.nav');
obs_file = fullfile(repo_dir, 'data', '1_8', 'A_1_8_1.obs');

start_time = datetime(2026, 1, 8, 5, 59, 39.116, 'TimeZone', 'UTC');
duration_seconds = 5;

% Fixed point used as ground truth receiver location for rho_gt.
% This default is the GPS-only OBS-estimated position from the same dataset.

base_lat_deg = 35.77204196;
base_lon_deg = 120.02477210;
base_height_m = 92.889;

north_offset_m = 30;   % 往北偏移 50 m

lat_deg = base_lat_deg + north_offset_m / 111000;
lon_deg = base_lon_deg;
height_m = base_height_m;

% lat_deg = 35.77204196 + 0.00090127 * 0.5; %00090127对应往北100m + 0.00090127 * 5
% lon_deg = 120.02477210 ;  % 对应向东0.00110592
% height_m = 92.889;

systems = {'G', 'C', 'E', 'R', 'J'};
elevation_mask_deg = 0;
min_clock_sats = 4;
receiver_clock_estimator = 'mean';

plot_enable = true;
write_csv = false;
output_csv = fullfile(repo_dir, 'physical_consistency', 'outputs', ...
    'nav_obs_range_residuals.csv');
%% =====================================================

fprintf('NAV file       : %s\n', nav_file);
fprintf('OBS file       : %s\n', obs_file);
fprintf('Start time UTC : %s\n', datestr(start_time, 'yyyy-mm-dd HH:MM:SS.FFF'));
fprintf('Duration       : %.3f s\n', duration_seconds);
fprintf('Fixed point    : %.8f deg, %.8f deg, %.3f m\n', lat_deg, lon_deg, height_m);
fprintf('Systems        : %s\n', strjoin(cellstr(string(systems)), ', '));

[residual_tbl, clock_tbl] = calculate_nav_obs_range_residuals( ...
    nav_file, obs_file, lat_deg, lon_deg, height_m, ...
    'systems', systems, ...
    'elevation_mask_deg', elevation_mask_deg, ...
    'min_clock_sats', min_clock_sats, ...
    'receiver_clock_estimator', receiver_clock_estimator, ...
    'start_time', start_time, ...
    'duration_seconds', duration_seconds);

ok_tbl = residual_tbl(residual_tbl.status == "ok", :);
summary_tbl = summarize_residuals_by_satellite_local(ok_tbl);
system_summary_tbl = summarize_residuals_by_system_local(ok_tbl);

fprintf('\nReceiver clock estimates:\n');
fprintf('  total rows: %d\n', height(clock_tbl));
clock_preview_idx = preview_indices_local(height(clock_tbl), 10, 5);
disp(clock_tbl(clock_preview_idx, {'epoch_idx', 'time', 'system', 'used_count', ...
    'receiver_clock_bias_m', 'clock_sample_std_m', 'status'}));

fprintf('\nFirst aligned residual rows:\n');
disp(ok_tbl(1:min(20, height(ok_tbl)), {'epoch_idx', 'time', 'sat_id', ...
    'pseudorange_code', 'pseudorange_raw_m', 'rho_gt_m', ...
    'r_obs_corrected_m', 'residual_obs_minus_gt_m', ...
    'receiver_clock_bias_m', 'sat_clock_correction_m'}));

fprintf('\nResidual summary by satellite:\n');
disp(summary_tbl);

fprintf('\nResidual summary by system:\n');
disp(system_summary_tbl);

fprintf('\nOverall residual statistics:\n');
fprintf('  rows        : %d\n', height(ok_tbl));
fprintf('  satellites  : %d\n', numel(unique(ok_tbl.sat_id)));
fprintf('  mean        : %.3f m\n', mean(ok_tbl.residual_obs_minus_gt_m, 'omitnan'));
fprintf('  std         : %.3f m\n', std(ok_tbl.residual_obs_minus_gt_m, 0, 'omitnan'));
fprintf('  rms         : %.3f m\n', sqrt(mean(ok_tbl.residual_obs_minus_gt_m.^2, 'omitnan')));
fprintf('  min / max   : %.3f / %.3f m\n', ...
    min(ok_tbl.residual_obs_minus_gt_m), max(ok_tbl.residual_obs_minus_gt_m));

if write_csv
    out_dir = fileparts(output_csv);
    if ~isempty(out_dir) && ~exist(out_dir, 'dir')
        mkdir(out_dir);
    end
    writetable(residual_tbl, output_csv);
    fprintf('Saved residual table: %s\n', output_csv);
end

if plot_enable && ~isempty(ok_tbl)
    plot_residuals_local(ok_tbl, start_time, lat_deg, lon_deg, height_m, systems);
end

%% ===================== local helpers =====================
function summary_tbl = summarize_residuals_by_satellite_local(ok_tbl)
if isempty(ok_tbl)
    summary_tbl = table();
    return;
end

sat_ids = unique(ok_tbl.sat_id, 'stable');
n = numel(sat_ids);
systems = strings(n, 1);
count = zeros(n, 1);
mean_m = NaN(n, 1);
std_m = NaN(n, 1);
rms_m = NaN(n, 1);
min_m = NaN(n, 1);
max_m = NaN(n, 1);
ptp_m = NaN(n, 1);

for i = 1:n
    mask = ok_tbl.sat_id == sat_ids(i);
    residual = ok_tbl.residual_obs_minus_gt_m(mask);
    systems(i) = extractBefore(sat_ids(i), 2);
    count(i) = numel(residual);
    mean_m(i) = mean(residual, 'omitnan');
    std_m(i) = std(residual, 0, 'omitnan');
    rms_m(i) = sqrt(mean(residual.^2, 'omitnan'));
    min_m(i) = min(residual);
    max_m(i) = max(residual);
    ptp_m(i) = max_m(i) - min_m(i);
end

summary_tbl = table(sat_ids, systems, count, mean_m, std_m, rms_m, min_m, max_m, ptp_m, ...
    'VariableNames', {'sat_id', 'system', 'count', 'mean_m', 'std_m', 'rms_m', ...
    'min_m', 'max_m', 'peak_to_peak_m'});
summary_tbl = sortrows(summary_tbl, 'rms_m', 'descend');
end

function summary_tbl = summarize_residuals_by_system_local(ok_tbl)
if isempty(ok_tbl)
    summary_tbl = table();
    return;
end

systems = unique(ok_tbl.system, 'stable');
n = numel(systems);
satellite_count = zeros(n, 1);
row_count = zeros(n, 1);
mean_m = NaN(n, 1);
std_m = NaN(n, 1);
rms_m = NaN(n, 1);
min_m = NaN(n, 1);
max_m = NaN(n, 1);

for i = 1:n
    mask = ok_tbl.system == systems(i);
    residual = ok_tbl.residual_obs_minus_gt_m(mask);
    satellite_count(i) = numel(unique(ok_tbl.sat_id(mask)));
    row_count(i) = numel(residual);
    mean_m(i) = mean(residual, 'omitnan');
    std_m(i) = std(residual, 0, 'omitnan');
    rms_m(i) = sqrt(mean(residual.^2, 'omitnan'));
    min_m(i) = min(residual);
    max_m(i) = max(residual);
end

summary_tbl = table(systems, satellite_count, row_count, mean_m, std_m, rms_m, min_m, max_m, ...
    'VariableNames', {'system', 'satellite_count', 'row_count', 'mean_m', ...
    'std_m', 'rms_m', 'min_m', 'max_m'});
summary_tbl = sortrows(summary_tbl, 'system');
end

function idx = preview_indices_local(n, head_n, tail_n)
if n <= head_n + tail_n
    idx = 1:n;
else
    idx = unique([1:head_n, (n - tail_n + 1):n], 'stable');
end
end

function plot_residuals_local(ok_tbl, start_time, lat_deg, lon_deg, height_m, systems)
fig_h = figure('Name', 'OBS corrected range minus NAV GT range', ...
    'Color', 'w', ...
    'Position', [80, 80, 1350, 760]);
ax = axes(fig_h);
hold(ax, 'on');

sat_ids = unique(ok_tbl.sat_id, 'stable');
for i = 1:numel(sat_ids)
    sat_id = sat_ids(i);
    mask = ok_tbl.sat_id == sat_id;
    x = seconds(ok_tbl.time(mask) - start_time);
    y = ok_tbl.residual_obs_minus_gt_m(mask);
    plot(ax, x, y, '-o', ...
        'LineWidth', 1.1, ...
        'MarkerSize', 3.5, ...
        'DisplayName', char(sat_id));
end

yline(ax, 0, ':', 'Color', [0.35 0.35 0.35], 'LineWidth', 0.8);
grid(ax, 'on');
box(ax, 'on');
xlabel(ax, 'Seconds after start time');
ylabel(ax, 'r_{obs} - \rho_{GT} (m)');
title(ax, sprintf('Corrected OBS range residuals | systems: %s', ...
    strjoin(cellstr(string(systems)), ',')));
subtitle(ax, sprintf('Fixed point: %.8f deg, %.8f deg, %.3f m', ...
    lat_deg, lon_deg, height_m));
legend(ax, 'show', 'Location', 'eastoutside');
hold(ax, 'off');
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
        error('example_nav_obs_range_residuals:RepoRootNotFound', ...
            'Cannot locate SatLock repo root from %s.', start_dir);
    end
    repo_dir = parent_dir;
end
end
