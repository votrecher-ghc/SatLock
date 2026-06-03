function [sky_tbl, fig_h] = calculate_and_plot_skyplot_from_nav(nav_input, target_time, lat_deg, lon_deg, height_m, varargin)
% CALCULATE_AND_PLOT_SKYPLOT_FROM_NAV
% Build a one-epoch GNSS sky plot from a RINEX navigation file and a fixed
% receiver location. This function does not require receiver observations,
% pseudorange, carrier phase, or SNR data.
%
% Example:
%   nav_file = 'D:\Matproject\SatLock\data\1_8\2026_1_8.nav';
%   t_utc = datetime(2026,1,8,4,0,0,'TimeZone','UTC');
%   [sky_tbl, fig_h] = calculate_and_plot_skyplot_from_nav( ...
%       nav_file, t_utc, 31.2304, 121.4737, 10);
%
% Name-value options:
%   'time_zone'                 Time zone for timezone-free input times.
%                               Default: 'UTC'.
%   'systems'                   Cell/string array of system letters.
%                               Default: {'G','C','E','R','J'}.
%   'elevation_mask_deg'         Elevation cutoff for plotting. Default: 0.
%   'plot_enable'               Draw the sky plot. Default: true.
%   'apply_sagnac'              Apply Earth-rotation correction. Default: true.
%   'max_broadcast_age_s'        Max ephemeris age for G/C/E/J. Default: 14400.
%   'max_glonass_age_s'          Max ephemeris age for R. Default: 3600.
%   'max_iterations'             Signal travel-time iterations. Default: 6.
%   'range_tolerance_m'          Iteration convergence tolerance. Default: 1e-3.
%   'save_path'                 Optional output PNG/PDF path.
%   'resolution'                Export resolution. Default: 300.

if nargin < 5 || isempty(height_m)
    height_m = 0;
end

repo_dir = resolve_repo_root_local(fileparts(mfilename('fullpath')));
ensure_project_paths_local(repo_dir);

cfg = default_cfg_local();
cfg = merge_options_local(cfg, varargin{:});

t_utc = normalize_time_local(target_time, cfg.time_zone);
nav_data = load_nav_data_local(nav_input);
receiver_ecef = geodetic_to_ecef_local(lat_deg, lon_deg, height_m);

system_order = {'G', 'C', 'R', 'E', 'J'};
rows = repmat(empty_sky_row_local(), 63 * numel(system_order), 1);
row_count = 0;

for sys_idx = 1:min(size(nav_data, 2), numel(system_order))
    sys_char = system_order{sys_idx};
    if ~system_enabled_local(sys_char, cfg.systems)
        continue;
    end

    for prn = 1:size(nav_data, 1)
        if isempty(nav_data{prn, sys_idx})
            continue;
        end

        row_count = row_count + 1;
        if row_count > numel(rows)
            rows(row_count + 64) = empty_sky_row_local(); %#ok<AGROW>
        end

        sat_id = sprintf('%s%02d', sys_char, prn);
        rows(row_count).sat_id = string(sat_id);
        rows(row_count).system = system_name_local(sys_char);
        rows(row_count).prn = prn;

        eph_age_s = ephemeris_age_seconds_local(nav_data{prn, sys_idx}, sys_char, t_utc);
        rows(row_count).ephemeris_age_s = eph_age_s;

        max_age_s = cfg.max_broadcast_age_s;
        if strcmp(sys_char, 'R')
            max_age_s = cfg.max_glonass_age_s;
        end

        if isfinite(max_age_s) && isfinite(eph_age_s) && eph_age_s > max_age_s
            rows(row_count).status = "ephemeris_too_old";
            continue;
        end

        try
            sat_res = compute_satellite_az_el_local( ...
                t_utc, receiver_ecef, lat_deg, lon_deg, sat_id, nav_data, cfg);

            rows(row_count).azimuth_deg = sat_res.azimuth_deg;
            rows(row_count).elevation_deg = sat_res.elevation_deg;
            rows(row_count).range_m = sat_res.range_m;
            rows(row_count).transit_time_s = sat_res.transit_time_s;
            rows(row_count).sat_clock_error_s = sat_res.sat_clock_error_s;
            rows(row_count).visible = sat_res.elevation_deg >= cfg.elevation_mask_deg;
            if rows(row_count).visible
                rows(row_count).status = "visible";
            else
                rows(row_count).status = "below_elevation_mask";
            end
        catch ME
            rows(row_count).status = "failed: " + string(ME.message);
        end
    end
end

if row_count == 0
    sky_tbl = empty_sky_table_local();
else
    sky_tbl = struct2table(rows(1:row_count));
    sky_tbl = sortrows(sky_tbl, {'system', 'prn'});
end

fig_h = gobjects(0);
if cfg.plot_enable
    fig_h = plot_sky_table_local(sky_tbl, t_utc, lat_deg, lon_deg, height_m, cfg);
    if strlength(string(cfg.save_path)) > 0
        save_figure_local(fig_h, char(string(cfg.save_path)), cfg.resolution);
    end
end
end

function cfg = default_cfg_local()
cfg = struct();
cfg.time_zone = 'UTC';
cfg.systems = {'G', 'C', 'E', 'R', 'J'};
cfg.elevation_mask_deg = 0;
cfg.plot_enable = true;
cfg.apply_sagnac = true;
cfg.max_broadcast_age_s = 4 * 3600;
cfg.max_glonass_age_s = 3600;
cfg.max_iterations = 6;
cfg.range_tolerance_m = 1e-3;
cfg.initial_range_m = 26560000;
cfg.save_path = "";
cfg.resolution = 300;
end

function cfg = merge_options_local(cfg, varargin)
if isempty(varargin)
    return;
end

if numel(varargin) == 1 && isstruct(varargin{1})
    updates = varargin{1};
    keys = fieldnames(updates);
    for i = 1:numel(keys)
        cfg = set_option_local(cfg, keys{i}, updates.(keys{i}));
    end
    return;
end

if mod(numel(varargin), 2) ~= 0
    error('calculate_and_plot_skyplot_from_nav:InvalidOptions', ...
        'Options must be provided as name-value pairs or one struct.');
end

for i = 1:2:numel(varargin)
    cfg = set_option_local(cfg, varargin{i}, varargin{i + 1});
end
end

function cfg = set_option_local(cfg, key, value)
key_norm = lower(regexprep(char(string(key)), '[^a-zA-Z0-9]', ''));
switch key_norm
    case 'timezone'
        cfg.time_zone = char(string(value));
    case 'systems'
        cfg.systems = normalize_system_list_local(value);
    case {'elevationmaskdeg', 'minelevationdeg', 'maskdeg'}
        cfg.elevation_mask_deg = double(value);
    case {'plotenable', 'plot'}
        cfg.plot_enable = logical(value);
    case 'applysagnac'
        cfg.apply_sagnac = logical(value);
    case {'maxbroadcastages', 'maxephemerisages', 'maxephages'}
        cfg.max_broadcast_age_s = double(value);
    case 'maxglonassages'
        cfg.max_glonass_age_s = double(value);
    case 'maxiterations'
        cfg.max_iterations = double(value);
    case 'rangetolerancem'
        cfg.range_tolerance_m = double(value);
    case 'initialrangem'
        cfg.initial_range_m = double(value);
    case 'savepath'
        cfg.save_path = string(value);
    case 'resolution'
        cfg.resolution = double(value);
    otherwise
        error('calculate_and_plot_skyplot_from_nav:UnknownOption', ...
            'Unknown option: %s', char(string(key)));
end
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
    error('calculate_and_plot_skyplot_from_nav:InvalidSystems', ...
        'systems must be a char, string array, or cell array.');
end
systems = systems(:).';
end

function t_utc = normalize_time_local(target_time, default_time_zone)
if isa(target_time, 'datetime')
    t_utc = target_time;
elseif isnumeric(target_time) && isscalar(target_time)
    t_utc = datetime(target_time, 'ConvertFrom', 'datenum');
elseif isnumeric(target_time) && numel(target_time) >= 6
    vec = double(target_time(:).');
    t_utc = datetime(vec(1), vec(2), vec(3), vec(4), vec(5), vec(6));
else
    txt = char(string(target_time));
    try
        t_utc = datetime(txt);
    catch
        t_utc = datetime(txt, 'InputFormat', 'yyyy-MM-dd HH:mm:ss');
    end
end

if isempty(t_utc.TimeZone)
    t_utc.TimeZone = char(string(default_time_zone));
end
t_utc.TimeZone = 'UTC';
end

function nav_data = load_nav_data_local(nav_input)
if iscell(nav_input)
    nav_data = nav_input;
    return;
end

nav_path = char(string(nav_input));
if ~exist(nav_path, 'file')
    error('calculate_and_plot_skyplot_from_nav:MissingNavFile', ...
        'Navigation file not found: %s', nav_path);
end

nav_data = parse_rinex_nav_multi_gnss(nav_path);
end

function sat_res = compute_satellite_az_el_local(t_utc, receiver_ecef, lat_deg, lon_deg, sat_id, nav_data, cfg)
c = 299792458.0;
rho = cfg.initial_range_m;
sat_clk_err = NaN;
sat_rx = [NaN; NaN; NaN];

for iter = 1:max(1, round(cfg.max_iterations))
    [sat_tx, ~, sat_clk_err] = calculate_satellite_state(t_utc, rho, sat_id, nav_data);
    tau = rho / c;

    if cfg.apply_sagnac
        sat_rx = rotate_satellite_to_receive_frame_local(sat_tx, tau);
    else
        sat_rx = sat_tx;
    end

    rho_new = norm(sat_rx - receiver_ecef);
    if abs(rho_new - rho) <= cfg.range_tolerance_m
        rho = rho_new;
        break;
    end
    rho = rho_new;
end

los_ecef = sat_rx - receiver_ecef;
[east, north, up] = ecef2enu(los_ecef(1), los_ecef(2), los_ecef(3), lat_deg, lon_deg, 0);
enu_norm = norm([east, north, up]);

az_deg = atan2d(east, north);
if az_deg < 0
    az_deg = az_deg + 360;
end
el_deg = asind(up / enu_norm);

sat_res = struct();
sat_res.azimuth_deg = az_deg;
sat_res.elevation_deg = el_deg;
sat_res.range_m = rho;
sat_res.transit_time_s = rho / c;
sat_res.sat_clock_error_s = sat_clk_err;
end

function sat_rx = rotate_satellite_to_receive_frame_local(sat_tx, transit_time_s)
omega_e = 7.2921151467e-5;
theta = omega_e * transit_time_s;
rot_z = [cos(theta), sin(theta), 0; ...
        -sin(theta), cos(theta), 0; ...
         0,          0,          1];
sat_rx = rot_z * sat_tx;
end

function age_s = ephemeris_age_seconds_local(eph_sets, sys_char, t_utc)
age_s = NaN;
if isempty(eph_sets)
    return;
end

switch sys_char
    case 'R'
        t_obs_num = datenum(t_utc);
        best_age = inf;
        for i = 1:numel(eph_sets)
            toc = eph_sets(i).Toc;
            t_ref = datenum([toc.Year, toc.Month, toc.Day, toc.Hour, toc.Minute, toc.Second]);
            this_age = abs((t_obs_num - t_ref) * 86400);
            if this_age < best_age
                best_age = this_age;
            end
        end
        age_s = best_age;

    otherwise
        [~, ~, ~, hh, mm, ss] = datevec(t_utc);
        t_sow = (weekday(t_utc) - 1) * 86400 + hh * 3600 + mm * 60 + ss;
        best_age = inf;
        for i = 1:numel(eph_sets)
            if ~isfield(eph_sets(i), 'Toe') || ~isfinite(eph_sets(i).Toe)
                continue;
            end
            this_age = abs(t_sow - eph_sets(i).Toe);
            if this_age > 302400
                this_age = 604800 - this_age;
            end
            if this_age < best_age
                best_age = this_age;
            end
        end
        age_s = best_age;
end
end

function fig_h = plot_sky_table_local(sky_tbl, t_utc, lat_deg, lon_deg, height_m, cfg)
fig_h = figure('Name', 'Sky plot from navigation ephemeris', 'Color', 'w');
pax = polaraxes(fig_h);
set(pax, ...
    'ThetaZeroLocation', 'top', ...
    'ThetaDir', 'clockwise', ...
    'RAxisLocation', 90, ...
    'RDir', 'reverse');
rlim(pax, [0 90]);
rticks(pax, [0 30 60 90]);
rticklabels(pax, {'90', '60', '30', '0'});
pax.RAxis.Label.String = 'Elevation (deg)';
title(pax, sprintf('Sky plot | %s UTC | %.6f, %.6f, %.1f m', ...
    datestr(t_utc, 'yyyy-mm-dd HH:MM:SS'), lat_deg, lon_deg, height_m));
hold(pax, 'on');
grid(pax, 'on');

if isempty(sky_tbl)
    hold(pax, 'off');
    warning('calculate_and_plot_skyplot_from_nav:NoSatellites', ...
        'No satellites were found in the navigation data.');
    return;
end

visible_mask = sky_tbl.visible & isfinite(sky_tbl.azimuth_deg) & isfinite(sky_tbl.elevation_deg);
visible_tbl = sky_tbl(visible_mask, :);

if isempty(visible_tbl)
    hold(pax, 'off');
    warning('calculate_and_plot_skyplot_from_nav:NoVisibleSatellites', ...
        'No satellites passed the elevation mask %.2f deg.', cfg.elevation_mask_deg);
    return;
end

legend_handles = gobjects(0);
legend_labels = {};
seen_systems = {};

for i = 1:height(visible_tbl)
    sys_char = char(extractBefore(visible_tbl.sat_id(i), 2));
    col = system_color_local(sys_char);
    theta = deg2rad(visible_tbl.azimuth_deg(i));
    radius = visible_tbl.elevation_deg(i);

    h = polarplot(pax, theta, radius, 'o', ...
        'MarkerSize', 10, ...
        'MarkerFaceColor', col, ...
        'MarkerEdgeColor', [0.15 0.15 0.15], ...
        'LineWidth', 0.8);
    text(pax, theta, radius, ['  ' char(visible_tbl.sat_id(i))], ...
        'FontSize', 9, ...
        'FontName', 'Arial', ...
        'Color', col, ...
        'FontWeight', 'bold');

    if ~any(strcmp(seen_systems, sys_char))
        seen_systems{end + 1} = sys_char; %#ok<AGROW>
        legend_handles(end + 1) = h; %#ok<AGROW>
        legend_labels{end + 1} = char(system_name_local(sys_char)); %#ok<AGROW>
    end
end

legend(pax, legend_handles, legend_labels, 'Location', 'eastoutside');
hold(pax, 'off');
end

function tf = system_enabled_local(sys_char, systems)
systems = normalize_system_list_local(systems);
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

function row = empty_sky_row_local()
row = struct( ...
    'sat_id', "", ...
    'system', "", ...
    'prn', NaN, ...
    'azimuth_deg', NaN, ...
    'elevation_deg', NaN, ...
    'range_m', NaN, ...
    'transit_time_s', NaN, ...
    'sat_clock_error_s', NaN, ...
    'ephemeris_age_s', NaN, ...
    'visible', false, ...
    'status', "");
end

function tbl = empty_sky_table_local()
tbl = struct2table(repmat(empty_sky_row_local(), 0, 1));
end

function save_figure_local(fig_h, out_path, resolution)
out_dir = fileparts(out_path);
if ~isempty(out_dir) && ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

try
    exportgraphics(fig_h, out_path, 'Resolution', resolution);
catch
    saveas(fig_h, out_path);
end
end

function ensure_project_paths_local(repo_dir)
addpath(fullfile(repo_dir, 'nav_parse'));
addpath(fullfile(repo_dir, 'calculate_clock_bias_and_positon'));
addpath(fullfile(repo_dir, 'plot', 'sky_plot'));
addpath(fullfile(repo_dir, 'physical_consistency'));
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
        error('calculate_and_plot_skyplot_from_nav:RepoRootNotFound', ...
            'Cannot locate SatLock repo root from %s.', start_dir);
    end
    repo_dir = parent_dir;
end
end
