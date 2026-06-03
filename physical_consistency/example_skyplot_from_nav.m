% Simple sky plot entry.
%
% Three supported workflows:
%   A) Give NAV + OBS, no manual position:
%      The script estimates receiver lat/lon/height from OBS, draws the
%      NAV sky plot at that estimated position, then compares with the OBS
%      sky plot from the nearest observation epoch.
%
%   B) Give NAV + OBS + manual position:
%      The script draws the NAV sky plot at your manual position, and compares
%      it with the receiver OBS sky plot computed from the OBS-estimated
%      receiver position.
%
%   C) Give NAV + manual position, no OBS:
%      The script draws only the NAV sky plot at the given time and position.

clear;
clc;
close all;

repo_dir = resolve_repo_root_local(fileparts(mfilename('fullpath')));
addpath(genpath(repo_dir));

%% ===================== Edit here =====================
nav_file = fullfile(repo_dir, 'data', '1_8', '2026_1_8.nav');

% Use UTC by default. For Beijing time, use TimeZone = 'Asia/Shanghai'.
target_time = datetime(2026, 1, 8, 5, 59, 39.116, 'TimeZone', 'UTC');

% Optional receiver observation file.
% Set obs_file = "" if you only want the NAV sky plot at a manual position.
obs_file = fullfile(repo_dir, 'data', '1_8', 'A_1_8_1.obs');
% obs_file = "";

% Manual position. Used when use_manual_position = true, or when obs_file = "".
use_manual_position = true;
manual_lat_deg = 31.2304;
manual_lon_deg = 121.4737;
manual_height_m = 10;

systems = {'G', 'C', 'E', 'R', 'J'};
elevation_mask_deg = 0;

% Used only for estimating receiver position from OBS. The current position
% solver has one receiver clock term, so single-system positioning is more
% stable than mixing all constellations directly.
position_system_priority = {'G', 'E', 'C'};
%% =====================================================

target_time = to_utc_local(target_time);

fprintf('NAV file : %s\n', nav_file);
fprintf('Time UTC : %s\n', datestr(target_time, 'yyyy-mm-dd HH:MM:SS.FFF'));

nav_data = parse_rinex_nav_multi_gnss(nav_file);

obs_data = [];
obs_position = struct();
obs_epoch_time = NaT('TimeZone', 'UTC');
has_obs = strlength(string(obs_file)) > 0;

if has_obs
    if ~isfile(obs_file)
        error('example_skyplot_from_nav:MissingObsFile', ...
            'OBS file does not exist: %s', obs_file);
    end

    fprintf('OBS file : %s\n', obs_file);
    obs_data = parse_rinex_obs(obs_file);

    [obs_position, obs_epoch_time] = estimate_receiver_position_from_obs_local( ...
        obs_data, nav_data, target_time, position_system_priority);

    fprintf('\nOBS-estimated receiver position:\n');
    fprintf('  lat    : %.8f deg\n', obs_position.lat_deg);
    fprintf('  lon    : %.8f deg\n', obs_position.lon_deg);
    fprintf('  height : %.3f m\n', obs_position.height_m);
    fprintf('  epochs used for position estimate: %d\n', obs_position.success_count);
    fprintf('  positioning system: %s\n', obs_position.solution_system);
end

if use_manual_position
    nav_position = struct( ...
        'source', "manual", ...
        'lat_deg', manual_lat_deg, ...
        'lon_deg', manual_lon_deg, ...
        'height_m', manual_height_m);
elseif has_obs
    nav_position = struct( ...
        'source', "obs_estimated", ...
        'lat_deg', obs_position.lat_deg, ...
        'lon_deg', obs_position.lon_deg, ...
        'height_m', obs_position.height_m);
else
    error('example_skyplot_from_nav:MissingPosition', ...
        'No OBS file was provided, so use_manual_position must be true.');
end

fprintf('\nNAV sky plot position (%s):\n', nav_position.source);
fprintf('  lat    : %.8f deg\n', nav_position.lat_deg);
fprintf('  lon    : %.8f deg\n', nav_position.lon_deg);
fprintf('  height : %.3f m\n', nav_position.height_m);

%% 1. NAV sky plot at selected position
[nav_sky_tbl, fig_nav] = calculate_and_plot_skyplot_from_nav( ...
    nav_data, target_time, nav_position.lat_deg, nav_position.lon_deg, nav_position.height_m, ...
    'systems', systems, ...
    'elevation_mask_deg', elevation_mask_deg);
fig_nav.Name = 'NAV sky plot';

nav_visible = nav_sky_tbl.visible & isfinite(nav_sky_tbl.azimuth_deg) & isfinite(nav_sky_tbl.elevation_deg);

fprintf('\nNAV visible satellites: %d\n', sum(nav_visible));
disp(nav_sky_tbl(nav_visible, {'sat_id', 'azimuth_deg', 'elevation_deg', 'status'}));

%% 2. Optional OBS sky plot and NAV-vs-OBS comparison
if has_obs
    [obs_sky_tbl, obs_epoch_time] = calculate_obs_epoch_sky_local( ...
        obs_data, nav_data, target_time, ...
        obs_position.lat_deg, obs_position.lon_deg, obs_position.height_m, ...
        systems, elevation_mask_deg);

    fig_obs = plot_receiver_observed_sky_local( ...
        obs_sky_tbl, obs_epoch_time, ...
        obs_position.lat_deg, obs_position.lon_deg, obs_position.height_m, elevation_mask_deg); %#ok<NASGU>

    fig_compare = plot_nav_obs_overlay_local( ...
        nav_sky_tbl, obs_sky_tbl, target_time, obs_epoch_time, ...
        nav_position, obs_position, elevation_mask_deg); %#ok<NASGU>

    obs_visible = obs_sky_tbl.visible & isfinite(obs_sky_tbl.azimuth_deg) & isfinite(obs_sky_tbl.elevation_deg);
    nav_ids = nav_sky_tbl.sat_id(nav_visible);
    obs_ids = obs_sky_tbl.sat_id(obs_visible);

    matched_ids = intersect(nav_ids, obs_ids);
    nav_only_ids = setdiff(nav_ids, obs_ids);
    obs_only_ids = setdiff(obs_ids, nav_ids);

    fprintf('\n================ NAV vs OBS ================\n');
    fprintf('Nearest OBS epoch    : %s UTC\n', datestr(obs_epoch_time, 'yyyy-mm-dd HH:MM:SS.FFF'));
    fprintf('Epoch time difference: %.3f s\n', abs(seconds(obs_epoch_time - target_time)));
    fprintf('NAV position source  : %s\n', nav_position.source);
    fprintf('OBS position source  : obs_estimated\n');
    fprintf('NAV visible count    : %d\n', numel(nav_ids));
    fprintf('OBS visible count    : %d\n', numel(obs_ids));
    fprintf('Matched count        : %d\n', numel(matched_ids));
    fprintf('NAV-only satellites  : %s\n', join_ids_local(nav_only_ids));
    fprintf('OBS-only satellites  : %s\n', join_ids_local(obs_only_ids));
    fprintf('============================================\n\n');

    fprintf('Receiver OBS visible satellites:\n');
    disp(obs_sky_tbl(obs_visible, {'sat_id', 'azimuth_deg', 'elevation_deg', 'pseudorange_code', 'status'}));
end

%% ===================== local helpers =====================
function [pos, obs_epoch_time] = estimate_receiver_position_from_obs_local(obs_data, nav_data, target_time, system_priority)
epoch_times = to_utc_local([obs_data.time]);
[~, nearest_idx] = min(abs(seconds(epoch_times - target_time)));
obs_epoch_time = epoch_times(nearest_idx);

half_window = 5;
candidate_idx = max(1, nearest_idx - half_window):min(numel(obs_data), nearest_idx + half_window);

ecef_list = [];
success_count = 0;
solution_system = "";
for s = 1:numel(system_priority)
    sys_char = char(string(system_priority{s}));
    ecef_list = [];
    success_count = 0;

    for k = candidate_idx
        try
            obs_one = keep_epoch_system_only_local(obs_data, k, sys_char);
            [ecef_pos, ~, ~] = calculate_receiver_position(obs_one, nav_data, k);
            if numel(ecef_pos) == 3 && all(isfinite(ecef_pos(:)))
                success_count = success_count + 1;
                ecef_list(success_count, :) = ecef_pos(:).'; %#ok<AGROW>
            end
        catch
            continue;
        end
    end

    if ~isempty(ecef_list)
        solution_system = string(sys_char);
        break;
    end
end

if isempty(ecef_list)
    error('example_skyplot_from_nav:ObsPositionFailed', ...
        'Failed to estimate receiver position from OBS near %s UTC using systems: %s.', ...
        datestr(target_time, 'yyyy-mm-dd HH:MM:SS.FFF'), strjoin(cellstr(string(system_priority)), ', '));
end

ecef_mean = mean(ecef_list, 1);
[lat_deg, lon_deg, height_m] = ecef2geodetic(ecef_mean(1), ecef_mean(2), ecef_mean(3));

pos = struct();
pos.source = "obs_estimated";
pos.lat_deg = lat_deg;
pos.lon_deg = lon_deg;
pos.height_m = height_m;
pos.ecef_m = ecef_mean(:);
pos.success_count = success_count;
pos.solution_system = solution_system;
pos.nearest_epoch_index = nearest_idx;
pos.nearest_epoch_time = obs_epoch_time;
end

function obs_one = keep_epoch_system_only_local(obs_data, epoch_idx, sys_char)
obs_one = obs_data;
sat_ids = fieldnames(obs_one(epoch_idx).data);
for i = 1:numel(sat_ids)
    sat_id = sat_ids{i};
    if ~startsWith(sat_id, sys_char, 'IgnoreCase', true)
        obs_one(epoch_idx).data = rmfield(obs_one(epoch_idx).data, sat_id);
    end
end
end

function [obs_tbl, epoch_time_utc] = calculate_obs_epoch_sky_local( ...
    obs_data, nav_data, target_time, lat_deg, lon_deg, height_m, systems, elevation_mask_deg)

epoch_times = to_utc_local([obs_data.time]);
[~, epoch_idx] = min(abs(seconds(epoch_times - target_time)));

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
fig_h = figure('Name', 'Receiver OBS sky plot', 'Color', 'w');
pax = setup_sky_axes_local(fig_h);
title(pax, sprintf('Receiver OBS sky plot | %s UTC | mask %.1f deg', ...
    datestr(epoch_time_utc, 'yyyy-mm-dd HH:MM:SS.FFF'), elevation_mask_deg));
subtitle(pax, sprintf('OBS-estimated receiver: %.6f, %.6f, %.1f m', lat_deg, lon_deg, height_m));
hold(pax, 'on');

visible_mask = obs_tbl.visible & isfinite(obs_tbl.azimuth_deg) & isfinite(obs_tbl.elevation_deg);
plot_sky_points_local(pax, obs_tbl(visible_mask, :), 's', true);

hold(pax, 'off');
end

function fig_h = plot_nav_obs_overlay_local(nav_tbl, obs_tbl, target_time, obs_epoch_time, nav_position, obs_position, elevation_mask_deg)
fig_h = figure('Name', 'NAV prediction vs receiver OBS', 'Color', 'w');
pax = setup_sky_axes_local(fig_h);
title(pax, sprintf('NAV prediction vs OBS | NAV %s UTC | OBS %s UTC', ...
    datestr(target_time, 'yyyy-mm-dd HH:MM:SS.FFF'), ...
    datestr(obs_epoch_time, 'yyyy-mm-dd HH:MM:SS.FFF')));
subtitle(pax, sprintf('NAV pos (%s): %.6f, %.6f, %.1f m | OBS pos: %.6f, %.6f, %.1f m | mask %.1f deg', ...
    nav_position.source, nav_position.lat_deg, nav_position.lon_deg, nav_position.height_m, ...
    obs_position.lat_deg, obs_position.lon_deg, obs_position.height_m, elevation_mask_deg));
hold(pax, 'on');

nav_mask = nav_tbl.visible & isfinite(nav_tbl.azimuth_deg) & isfinite(nav_tbl.elevation_deg);
obs_mask = obs_tbl.visible & isfinite(obs_tbl.azimuth_deg) & isfinite(obs_tbl.elevation_deg);

plot_sky_points_local(pax, nav_tbl(nav_mask, :), 'o', false);
plot_sky_points_local(pax, obs_tbl(obs_mask, :), 's', true);

h_nav = polarplot(pax, NaN, NaN, 'o', 'MarkerSize', 10, ...
    'MarkerFaceColor', 'none', 'MarkerEdgeColor', [0.55 0.55 0.55], 'LineWidth', 1.2);
h_obs = polarplot(pax, NaN, NaN, 's', 'MarkerSize', 10, ...
    'MarkerFaceColor', [0.10 0.45 0.85], 'MarkerEdgeColor', [0.10 0.10 0.10], 'LineWidth', 0.8);
legend(pax, [h_nav h_obs], {'NAV predicted', 'Receiver OBS'}, 'Location', 'eastoutside');

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

function plot_sky_points_local(pax, tbl, marker, use_system_color)
for i = 1:height(tbl)
    sat_id = char(tbl.sat_id(i));
    sys_char = sat_id(1);
    if use_system_color
        col = system_color_local(sys_char);
        face_col = col;
        label_enable = true;
    else
        col = [0.55 0.55 0.55];
        face_col = 'none';
        label_enable = true;
    end

    theta = deg2rad(tbl.azimuth_deg(i));
    radius = tbl.elevation_deg(i);
    polarplot(pax, theta, radius, marker, ...
        'MarkerSize', 10, ...
        'MarkerFaceColor', face_col, ...
        'MarkerEdgeColor', col, ...
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

switch upper(sys_char)
    case {'G', 'E', 'J'}
        preferred = {'C1C', 'C1X', 'C1B', 'C2L', 'C7Q'};
    case 'C'
        preferred = {'C2I', 'C1I', 'C1X', 'C7I'};
    case 'R'
        preferred = {'C1C', 'C1P', 'C2C', 'C2P'};
    otherwise
        preferred = {};
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
systems = cellstr(upper(string(systems)));
tf = any(strcmpi(sys_char, systems));
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

function repo_dir = resolve_repo_root_local(start_dir)
repo_dir = start_dir;
while true
    if isfolder(fullfile(repo_dir, 'data')) && ...
            isfolder(fullfile(repo_dir, 'nav_parse')) && ...
            isfolder(fullfile(repo_dir, 'calculate_clock_bias_and_positon'))
        return;
    end

    parent_dir = fileparts(repo_dir);
    if strcmp(parent_dir, repo_dir) || isempty(parent_dir)
        error('example_skyplot_from_nav:RepoRootNotFound', ...
            'Cannot locate SatLock repo root from %s.', start_dir);
    end
    repo_dir = parent_dir;
end
end
