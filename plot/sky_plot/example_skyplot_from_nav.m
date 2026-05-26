% Example: compare a nav-only sky plot against receiver-observed satellites.
%
% This script draws three figures:
%   1) Nav-only prediction from the broadcast ephemeris.
%   2) Receiver-observed satellites from the nearest OBS epoch.
%   3) Overlay comparison: nav prediction vs receiver observation.

clear;
clc;
close all;

repo_dir = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(genpath(repo_dir));

nav_file = fullfile(repo_dir, 'data', '1_8', '2026_1_8.nav');
obs_file = fullfile(repo_dir, 'data', '1_8', 'A_1_8_1.obs');

% Default interpretation is UTC. If your input is Beijing time, use:
% target_time = datetime(2026,1,8,13,59,39.116,'TimeZone','Asia/Shanghai');
target_time = datetime(2026, 1, 8, 5, 59, 39.116, 'TimeZone', 'UTC');

lat_deg = 31.2304;
lon_deg = 121.4737;
height_m = 10;

systems = {'G', 'C', 'E', 'R', 'J'};
elevation_mask_deg = 0;

%% 1. Parse input data once
fprintf('Parsing NAV: %s\n', nav_file);
nav_data = parse_rinex_nav_multi_gnss(nav_file);

fprintf('Parsing OBS: %s\n', obs_file);
obs_data = parse_rinex_obs(obs_file);

%% 2. Nav-only prediction at the requested time
[nav_sky_tbl, fig_nav] = calculate_and_plot_skyplot_from_nav( ...
    nav_data, target_time, lat_deg, lon_deg, height_m, ...
    'elevation_mask_deg', elevation_mask_deg, ...
    'systems', systems);
fig_nav.Name = 'Nav-only predicted sky plot';

%% 3. Receiver-observed satellites at the nearest OBS epoch
[obs_sky_tbl, obs_epoch_time] = calculate_obs_epoch_sky_local( ...
    obs_data, nav_data, target_time, lat_deg, lon_deg, height_m, ...
    systems, elevation_mask_deg);
fig_obs = plot_receiver_observed_sky_local( ...
    obs_sky_tbl, obs_epoch_time, lat_deg, lon_deg, height_m, elevation_mask_deg);

%% 4. Overlay comparison
fig_overlay = plot_nav_obs_overlay_local( ...
    nav_sky_tbl, obs_sky_tbl, target_time, obs_epoch_time, ...
    lat_deg, lon_deg, height_m, elevation_mask_deg);

%% 5. Print a compact comparison
nav_visible = nav_sky_tbl.visible & isfinite(nav_sky_tbl.azimuth_deg) & isfinite(nav_sky_tbl.elevation_deg);
obs_visible = obs_sky_tbl.visible & isfinite(obs_sky_tbl.azimuth_deg) & isfinite(obs_sky_tbl.elevation_deg);

nav_ids = nav_sky_tbl.sat_id(nav_visible);
obs_ids = obs_sky_tbl.sat_id(obs_visible);
matched_ids = intersect(nav_ids, obs_ids);
nav_only_ids = setdiff(nav_ids, obs_ids);
obs_only_ids = setdiff(obs_ids, nav_ids);

fprintf('\n================ Sky Plot Comparison ================\n');
fprintf('Requested time       : %s UTC\n', datestr(to_utc_local(target_time), 'yyyy-mm-dd HH:MM:SS.FFF'));
fprintf('Nearest OBS epoch    : %s UTC\n', datestr(obs_epoch_time, 'yyyy-mm-dd HH:MM:SS.FFF'));
fprintf('Epoch time difference: %.3f s\n', abs(seconds(obs_epoch_time - to_utc_local(target_time))));
fprintf('Nav visible count    : %d\n', numel(nav_ids));
fprintf('OBS visible count    : %d\n', numel(obs_ids));
fprintf('Matched count        : %d\n', numel(matched_ids));
fprintf('Nav-only satellites  : %s\n', join_ids_local(nav_only_ids));
fprintf('OBS-only satellites  : %s\n', join_ids_local(obs_only_ids));
fprintf('=====================================================\n\n');

disp('Nav-only visible satellites:');
disp(nav_sky_tbl(nav_visible, {'sat_id', 'azimuth_deg', 'elevation_deg', 'status'}));

disp('Receiver-observed visible satellites:');
disp(obs_sky_tbl(obs_visible, {'sat_id', 'azimuth_deg', 'elevation_deg', 'pseudorange_code', 'status'}));

%% Local helpers
function [obs_tbl, epoch_time_utc] = calculate_obs_epoch_sky_local( ...
    obs_data, nav_data, target_time, lat_deg, lon_deg, height_m, systems, elevation_mask_deg)

epoch_times = [obs_data.time];
epoch_times = to_utc_local(epoch_times);
target_time_utc = to_utc_local(target_time);
[~, epoch_idx] = min(abs(seconds(epoch_times - target_time_utc)));

epoch = obs_data(epoch_idx);
epoch_time_utc = epoch_times(epoch_idx);
receiver_ecef = geodetic_to_ecef_local(lat_deg, lon_deg, height_m);
sat_ids = fieldnames(epoch.data);

rows = repmat(empty_obs_sky_row_local(), max(1, numel(sat_ids)), 1);
row_count = 0;

for i = 1:numel(sat_ids)
    sat_id = sat_ids{i};
    sys_char = upper(sat_id(1));
    if ~system_enabled_local(sys_char, systems)
        continue;
    end

    sat_obs = epoch.data.(sat_id);
    [pseudorange, pr_code] = select_pseudorange_local(sat_obs, sys_char);
    if ~isfinite(pseudorange) || pseudorange <= 0
        continue;
    end

    row_count = row_count + 1;
    rows(row_count).sat_id = string(sat_id);
    rows(row_count).system = system_name_local(sys_char);
    rows(row_count).prn = str2double(sat_id(2:end));
    rows(row_count).pseudorange_m = pseudorange;
    rows(row_count).pseudorange_code = string(pr_code);

    try
        [sat_pos_tx, ~, sat_clk_err] = calculate_satellite_state(epoch_time_utc, pseudorange, sat_id, nav_data);
        sat_pos_rx = rotate_satellite_to_receive_frame_local(sat_pos_tx, pseudorange / 299792458.0);
        los_ecef = sat_pos_rx - receiver_ecef;
        [east, north, up] = ecef2enu(los_ecef(1), los_ecef(2), los_ecef(3), lat_deg, lon_deg, 0);
        enu_norm = norm([east, north, up]);

        az_deg = atan2d(east, north);
        if az_deg < 0
            az_deg = az_deg + 360;
        end
        el_deg = asind(up / enu_norm);

        rows(row_count).azimuth_deg = az_deg;
        rows(row_count).elevation_deg = el_deg;
        rows(row_count).range_m = norm(los_ecef);
        rows(row_count).sat_clock_error_s = sat_clk_err;
        rows(row_count).visible = el_deg >= elevation_mask_deg;
        if rows(row_count).visible
            rows(row_count).status = "observed_visible";
        else
            rows(row_count).status = "observed_below_mask";
        end
    catch ME
        rows(row_count).status = "failed: " + string(ME.message);
    end
end

if row_count == 0
    obs_tbl = struct2table(repmat(empty_obs_sky_row_local(), 0, 1));
else
    obs_tbl = sortrows(struct2table(rows(1:row_count)), {'system', 'prn'});
end
end

function fig_h = plot_receiver_observed_sky_local(obs_tbl, epoch_time_utc, lat_deg, lon_deg, height_m, elevation_mask_deg)
fig_h = figure('Name', 'Receiver-observed sky plot', 'Color', 'w');
pax = setup_sky_axes_local(fig_h);
title(pax, sprintf('Receiver OBS sky plot | %s UTC | mask %.1f deg', ...
    datestr(epoch_time_utc, 'yyyy-mm-dd HH:MM:SS.FFF'), elevation_mask_deg));
hold(pax, 'on');

visible_mask = obs_tbl.visible & isfinite(obs_tbl.azimuth_deg) & isfinite(obs_tbl.elevation_deg);
plot_sky_points_local(pax, obs_tbl(visible_mask, :), 'square', true, true);

subtitle(pax, sprintf('Receiver location: %.6f, %.6f, %.1f m', lat_deg, lon_deg, height_m));
hold(pax, 'off');
end

function fig_h = plot_nav_obs_overlay_local(nav_tbl, obs_tbl, target_time, obs_epoch_time, lat_deg, lon_deg, height_m, elevation_mask_deg)
fig_h = figure('Name', 'Nav prediction vs receiver observation', 'Color', 'w');
pax = setup_sky_axes_local(fig_h);
title(pax, sprintf('Nav prediction vs OBS | requested %s UTC | OBS %s UTC', ...
    datestr(to_utc_local(target_time), 'yyyy-mm-dd HH:MM:SS.FFF'), ...
    datestr(obs_epoch_time, 'yyyy-mm-dd HH:MM:SS.FFF')));
subtitle(pax, sprintf('Receiver location: %.6f, %.6f, %.1f m | mask %.1f deg', ...
    lat_deg, lon_deg, height_m, elevation_mask_deg));
hold(pax, 'on');

nav_mask = nav_tbl.visible & isfinite(nav_tbl.azimuth_deg) & isfinite(nav_tbl.elevation_deg);
obs_mask = obs_tbl.visible & isfinite(obs_tbl.azimuth_deg) & isfinite(obs_tbl.elevation_deg);

plot_sky_points_local(pax, nav_tbl(nav_mask, :), 'circle', false, false);
plot_sky_points_local(pax, obs_tbl(obs_mask, :), 'square', true, true);

h_nav = polarplot(pax, NaN, NaN, 'o', 'MarkerSize', 7, ...
    'MarkerFaceColor', 'none', 'MarkerEdgeColor', [0.55 0.55 0.55], 'LineWidth', 1.2);
h_obs = polarplot(pax, NaN, NaN, 's', 'MarkerSize', 7, ...
    'MarkerFaceColor', [0.10 0.45 0.85], 'MarkerEdgeColor', [0.10 0.10 0.10], 'LineWidth', 0.8);
legend(pax, [h_nav h_obs], {'Nav predicted', 'Receiver observed'}, 'Location', 'eastoutside');

hold(pax, 'off');
end

function pax = setup_sky_axes_local(fig_h)
pax = polaraxes(fig_h);
set(pax, ...
    'ThetaZeroLocation', 'top', ...
    'ThetaDir', 'clockwise', ...
    'RAxisLocation', 90, ...
    'RDir', 'reverse', ...
    'FontName', 'Arial');
rlim(pax, [0 90]);
rticks(pax, [0 30 60 90]);
rticklabels(pax, {'90', '60', '30', '0'});
pax.RAxis.Label.String = 'Elevation (deg)';
grid(pax, 'on');
end

function plot_sky_points_local(pax, tbl, marker_kind, use_system_color, label_enable)
if isempty(tbl)
    return;
end
if nargin < 5
    label_enable = true;
end

for i = 1:height(tbl)
    sat_id = char(tbl.sat_id(i));
    sys_char = sat_id(1);
    col = [0.55 0.55 0.55];
    if use_system_color
        col = system_color_local(sys_char);
    end

    theta = deg2rad(tbl.azimuth_deg(i));
    radius = tbl.elevation_deg(i);

    switch marker_kind
        case 'square'
            marker = 's';
            face_col = col;
            edge_col = [0.10 0.10 0.10];
            alpha_face = 0.95;
        otherwise
            marker = 'o';
            face_col = 'none';
            edge_col = col;
            alpha_face = 1.0; %#ok<NASGU>
    end

    polarplot(pax, theta, radius, marker, ...
        'MarkerSize', 7, ...
        'MarkerFaceColor', face_col, ...
        'MarkerEdgeColor', edge_col, ...
        'LineWidth', 1.1);
    if label_enable
        text(pax, theta, radius, ['  ' sat_id], ...
            'FontSize', 8, ...
            'FontName', 'Arial', ...
            'FontWeight', 'bold', ...
            'Color', col);
    end
end
end

function [pseudorange, pr_code] = select_pseudorange_local(sat_obs, sys_char)
pseudorange = NaN;
pr_code = "";
if ~isfield(sat_obs, 'pseudorange') || isempty(sat_obs.pseudorange)
    return;
end

preferred = {};
switch upper(sys_char)
    case {'G', 'E', 'J'}
        preferred = {'C1C', 'C1X', 'C1B', 'C2L', 'C7Q'};
    case 'C'
        preferred = {'C2I', 'C1I', 'C1X', 'C7I'};
    case 'R'
        preferred = {'C1C', 'C1P', 'C2C', 'C2P'};
end

for i = 1:numel(preferred)
    code = preferred{i};
    if isfield(sat_obs.pseudorange, code)
        val = sat_obs.pseudorange.(code);
        if isfinite(val) && val > 0
            pseudorange = val;
            pr_code = string(code);
            return;
        end
    end
end

codes = fieldnames(sat_obs.pseudorange);
for i = 1:numel(codes)
    val = sat_obs.pseudorange.(codes{i});
    if isfinite(val) && val > 0
        pseudorange = val;
        pr_code = string(codes{i});
        return;
    end
end
end

function sat_rx = rotate_satellite_to_receive_frame_local(sat_tx, transit_time_s)
omega_e = 7.2921151467e-5;
theta = omega_e * transit_time_s;
rot_z = [cos(theta), sin(theta), 0; ...
        -sin(theta), cos(theta), 0; ...
         0,          0,          1];
sat_rx = rot_z * sat_tx;
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

function t = to_utc_local(t)
if isempty(t.TimeZone)
    t.TimeZone = 'UTC';
end
t.TimeZone = 'UTC';
end

function tf = system_enabled_local(sys_char, systems)
systems = normalize_system_list_local(systems);
tf = any(strcmpi(sys_char, systems));
end

function systems = normalize_system_list_local(value)
if ischar(value)
    systems = cellstr(upper(value(:)));
elseif isstring(value)
    systems = cellstr(upper(value(:)));
elseif iscell(value)
    systems = cell(size(value));
    for i = 1:numel(value)
        systems{i} = upper(char(string(value{i})));
    end
else
    error('Invalid systems input.');
end
systems = systems(:).';
end

function name = system_name_local(sys_char)
switch upper(char(sys_char))
    case 'G'
        name = "GPS";
    case 'C'
        name = "BeiDou";
    case 'R'
        name = "GLONASS";
    case 'E'
        name = "Galileo";
    case 'J'
        name = "QZSS";
    otherwise
        name = string(sys_char);
end
end

function col = system_color_local(sys_char)
switch upper(char(sys_char))
    case 'G'
        col = [0.0000 0.4470 0.7410];
    case 'C'
        col = [0.8500 0.3250 0.0980];
    case 'E'
        col = [0.4660 0.6740 0.1880];
    case 'R'
        col = [0.4940 0.1840 0.5560];
    case 'J'
        col = [0.9290 0.6940 0.1250];
    otherwise
        col = [0.3500 0.3500 0.3500];
end
end

function row = empty_obs_sky_row_local()
row = struct( ...
    'sat_id', "", ...
    'system', "", ...
    'prn', NaN, ...
    'azimuth_deg', NaN, ...
    'elevation_deg', NaN, ...
    'range_m', NaN, ...
    'pseudorange_m', NaN, ...
    'pseudorange_code', "", ...
    'sat_clock_error_s', NaN, ...
    'visible', false, ...
    'status', "");
end

function txt = join_ids_local(ids)
if isempty(ids)
    txt = "(none)";
else
    txt = char(strjoin(cellstr(string(ids(:).')), ', '));
end
end
