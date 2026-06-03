function [range_tbl, meta] = calculate_satellite_ranges_from_nav(nav_input, target_time, lat_deg, lon_deg, height_m, varargin)
% CALCULATE_SATELLITE_RANGES_FROM_NAV
% Compute pure NAV-predicted satellite-to-point geometric ranges.
%
% This function does not use receiver OBS data, pseudorange, carrier phase,
% SNR, or any receiver measurement. It only needs broadcast navigation data,
% a target time, and a fixed geodetic point.
%
% Example:
%   nav_file = 'D:\Matproject\SatLock\data\1_8\2026_1_8.nav';
%   t_utc = datetime(2026, 1, 8, 5, 59, 39.116, 'TimeZone', 'UTC');
%   range_tbl = calculate_satellite_ranges_from_nav( ...
%       nav_file, t_utc, 35.77204197, 120.02477211, 92.889);
%
% Inputs:
%   nav_input    RINEX nav file path, or parsed nav_data returned by
%                parse_rinex_nav_multi_gnss.
%   target_time  datetime, datenum, date vector, or parseable time string.
%   lat_deg      Geodetic latitude in degrees.
%   lon_deg      Geodetic longitude in degrees.
%   height_m     Ellipsoidal height in meters. Empty defaults to 0.
%
% Extra options:
%   'visible_only'    Return only satellites above elevation_mask_deg.
%                     Default: false.
%   'include_failed'  Keep rows whose range could not be computed. Default:
%                     false.
%
% Other name-value options are forwarded to calculate_and_plot_skyplot_from_nav,
% including:
%   'time_zone', 'systems', 'elevation_mask_deg', 'apply_sagnac',
%   'max_broadcast_age_s', 'max_glonass_age_s', 'max_iterations',
%   'range_tolerance_m', and 'initial_range_m'.
%
% Output range_tbl columns:
%   sat_id, system, prn, azimuth_deg, elevation_deg, range_m, range_km,
%   transit_time_s, sat_clock_error_s, sat_clock_correction_m,
%   ephemeris_age_s, visible, status.
%
% Notes:
%   range_m is a NAV-predicted geometric range. It is computed with iterative
%   signal travel time and optional Earth-rotation (Sagnac) correction. It is
%   not an observed pseudorange and does not include receiver clock bias,
%   atmospheric delay, hardware delay, or multipath.

if nargin < 5 || isempty(height_m)
    height_m = 0;
end

validateattributes(lat_deg, {'numeric'}, {'scalar', 'real', 'finite', '>=', -90, '<=', 90}, ...
    mfilename, 'lat_deg');
validateattributes(lon_deg, {'numeric'}, {'scalar', 'real', 'finite'}, ...
    mfilename, 'lon_deg');
validateattributes(height_m, {'numeric'}, {'scalar', 'real', 'finite'}, ...
    mfilename, 'height_m');

option_pairs = options_to_pairs_local(varargin{:});
[visible_only, option_pairs] = pop_boolean_option_local( ...
    option_pairs, {'visibleonly', 'onlyvisible'}, false);
[include_failed, option_pairs] = pop_boolean_option_local( ...
    option_pairs, {'includefailed', 'keepfailed'}, false);

[sky_tbl, ~] = calculate_and_plot_skyplot_from_nav( ...
    nav_input, target_time, lat_deg, lon_deg, height_m, ...
    option_pairs{:}, 'plot_enable', false);

if isempty(sky_tbl)
    range_tbl = empty_range_table_local();
else
    keep_mask = true(height(sky_tbl), 1);

    if ~include_failed
        keep_mask = keep_mask & isfinite(sky_tbl.range_m);
    end

    if visible_only
        keep_mask = keep_mask & sky_tbl.visible;
    end

    range_tbl = sky_tbl(keep_mask, :);
    if ~isempty(range_tbl)
        c = 299792458.0;
        range_tbl.range_km = range_tbl.range_m / 1000;
        range_tbl.sat_clock_correction_m = range_tbl.sat_clock_error_s * c;
        range_tbl = reorder_columns_local(range_tbl);
    else
        range_tbl = empty_range_table_local();
    end
end

meta = struct();
meta.lat_deg = lat_deg;
meta.lon_deg = lon_deg;
meta.height_m = height_m;
meta.visible_only = visible_only;
meta.include_failed = include_failed;
meta.row_count = height(range_tbl);
meta.visible_count = sum(range_tbl.visible);
end

function option_pairs = options_to_pairs_local(varargin)
if isempty(varargin)
    option_pairs = {};
    return;
end

if numel(varargin) == 1 && isstruct(varargin{1})
    updates = varargin{1};
    keys = fieldnames(updates);
    option_pairs = cell(1, 2 * numel(keys));
    for i = 1:numel(keys)
        option_pairs{2 * i - 1} = keys{i};
        option_pairs{2 * i} = updates.(keys{i});
    end
    return;
end

if mod(numel(varargin), 2) ~= 0
    error('calculate_satellite_ranges_from_nav:InvalidOptions', ...
        'Options must be name-value pairs, or one scalar struct.');
end

option_pairs = varargin;
end

function [value, option_pairs] = pop_boolean_option_local(option_pairs, aliases, default_value)
value = default_value;
if isempty(option_pairs)
    return;
end

keep = true(1, numel(option_pairs));
for i = 1:2:numel(option_pairs)
    key_norm = lower(regexprep(char(string(option_pairs{i})), '[^a-zA-Z0-9]', ''));
    if any(strcmp(key_norm, aliases))
        value = logical(option_pairs{i + 1});
        keep(i:i + 1) = false;
    end
end

option_pairs = option_pairs(keep);
end

function tbl = reorder_columns_local(tbl)
preferred = {'sat_id', 'system', 'prn', ...
    'azimuth_deg', 'elevation_deg', ...
    'range_m', 'range_km', 'transit_time_s', ...
    'sat_clock_error_s', 'sat_clock_correction_m', ...
    'ephemeris_age_s', 'visible', 'status'};

existing = tbl.Properties.VariableNames;
ordered = preferred(ismember(preferred, existing));
remaining = existing(~ismember(existing, ordered));
tbl = tbl(:, [ordered, remaining]);
end

function tbl = empty_range_table_local()
tbl = table( ...
    string.empty(0, 1), ...
    string.empty(0, 1), ...
    zeros(0, 1), ...
    zeros(0, 1), ...
    zeros(0, 1), ...
    zeros(0, 1), ...
    zeros(0, 1), ...
    zeros(0, 1), ...
    zeros(0, 1), ...
    zeros(0, 1), ...
    zeros(0, 1), ...
    false(0, 1), ...
    string.empty(0, 1), ...
    'VariableNames', {'sat_id', 'system', 'prn', ...
    'azimuth_deg', 'elevation_deg', ...
    'range_m', 'range_km', 'transit_time_s', ...
    'sat_clock_error_s', 'sat_clock_correction_m', ...
    'ephemeris_age_s', 'visible', 'status'});
end
