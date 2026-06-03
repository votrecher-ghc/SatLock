function [residual_tbl, clock_tbl] = calculate_nav_obs_range_residuals(nav_input, obs_input, lat_deg, lon_deg, height_m, varargin)
% CALCULATE_NAV_OBS_RANGE_RESIDUALS
% Align NAV-predicted geometric ranges and OBS pseudorange-derived ranges.
%
% For each selected OBS epoch and satellite:
%   1. Compute NAV-only geometric range rho_gt_m at the fixed geodetic point.
%   2. Read raw pseudorange P_i from OBS.
%   3. Compute satellite clock error dt_s from broadcast NAV.
%   4. Estimate receiver clock bias dt_r from the fixed-position equation:
%        P_i = rho_i + c * (dt_r - dt_s) + residual_i
%   5. Build corrected observed range:
%        r_obs_i = P_i - c * dt_r + c * dt_s
%   6. Output:
%        residual_i = r_obs_i - rho_gt_i
%
% Important:
%   Receiver clock bias is estimated per epoch and per GNSS system by default.
%   This avoids mixing GPS/BDS/Galileo/GLONASS/QZSS inter-system biases into a
%   single clock term. Set systems = {'G'} if you want a GPS-only diagnostic.

if nargin < 5 || isempty(height_m)
    height_m = 0;
end

repo_dir = resolve_repo_root_local(fileparts(mfilename('fullpath')));
ensure_project_paths_local(repo_dir);

cfg = default_cfg_local();
cfg = merge_options_local(cfg, varargin{:});

nav_data = load_nav_data_local(nav_input);
obs_data = load_obs_data_local(obs_input);
epoch_indices = select_epoch_indices_local(obs_data, cfg);
receiver_ecef = geodetic_to_ecef_local(lat_deg, lon_deg, height_m);

c = 299792458.0;
rows = repmat(empty_residual_row_local(), 1024, 1);
row_count = 0;
clock_rows = repmat(empty_clock_row_local(), 256, 1);
clock_count = 0;

for idx_pos = 1:numel(epoch_indices)
    epoch_idx = epoch_indices(idx_pos);
    epoch_time = to_utc_local(obs_data(epoch_idx).time, cfg.time_zone);

    obs_candidates = build_epoch_candidates_local( ...
        obs_data(epoch_idx), epoch_time, epoch_idx, receiver_ecef, ...
        lat_deg, lon_deg, nav_data, cfg);

    for sys_i = 1:numel(cfg.systems)
        sys_char = upper(char(string(cfg.systems{sys_i})));
        sys_mask = strcmp({obs_candidates.system}, sys_char);
        sys_candidates = obs_candidates(sys_mask);

        valid_clock_mask = [sys_candidates.clock_sample_valid];
        clock_samples_m = [sys_candidates(valid_clock_mask).clock_sample_m];

        clock_count = clock_count + 1;
        if clock_count > numel(clock_rows)
            clock_rows(end + 1:end + 128, 1) = repmat(empty_clock_row_local(), 128, 1); %#ok<AGROW>
        end

        clock_rows(clock_count).epoch_idx = epoch_idx;
        clock_rows(clock_count).time = epoch_time;
        clock_rows(clock_count).system = string(sys_char);
        clock_rows(clock_count).candidate_count = numel(sys_candidates);
        clock_rows(clock_count).used_count = numel(clock_samples_m);

        if numel(clock_samples_m) < cfg.min_clock_sats
            clock_rows(clock_count).status = "not_enough_sats";
            receiver_clock_bias_m = NaN;
        else
            receiver_clock_bias_m = estimate_clock_bias_local(clock_samples_m, cfg.receiver_clock_estimator);
            clock_rows(clock_count).receiver_clock_bias_m = receiver_clock_bias_m;
            clock_rows(clock_count).receiver_clock_bias_s = receiver_clock_bias_m / c;
            clock_rows(clock_count).clock_sample_std_m = std(clock_samples_m, 0);
            clock_rows(clock_count).status = "ok";
        end

        for k = 1:numel(sys_candidates)
            cand = sys_candidates(k);
            row_count = row_count + 1;
            if row_count > numel(rows)
                rows(end + 1:end + 512, 1) = repmat(empty_residual_row_local(), 512, 1); %#ok<AGROW>
            end

            rows(row_count).epoch_idx = cand.epoch_idx;
            rows(row_count).time = cand.time;
            rows(row_count).sat_id = cand.sat_id;
            rows(row_count).system = string(cand.system);
            rows(row_count).pseudorange_code = cand.pseudorange_code;
            rows(row_count).pseudorange_raw_m = cand.pseudorange_raw_m;
            rows(row_count).rho_gt_m = cand.rho_gt_m;
            rows(row_count).azimuth_deg = cand.azimuth_deg;
            rows(row_count).elevation_deg = cand.elevation_deg;
            rows(row_count).visible = cand.visible;
            rows(row_count).sat_clock_error_s = cand.sat_clock_error_s;
            rows(row_count).sat_clock_correction_m = c * cand.sat_clock_error_s;
            rows(row_count).receiver_clock_bias_m = receiver_clock_bias_m;
            rows(row_count).receiver_clock_bias_s = receiver_clock_bias_m / c;
            rows(row_count).clock_sample_m = cand.clock_sample_m;
            rows(row_count).clock_sample_valid = cand.clock_sample_valid;

            if cand.clock_sample_valid && isfinite(receiver_clock_bias_m)
                r_obs = cand.pseudorange_raw_m - receiver_clock_bias_m + c * cand.sat_clock_error_s;
                rows(row_count).r_obs_corrected_m = r_obs;
                rows(row_count).residual_obs_minus_gt_m = r_obs - cand.rho_gt_m;
                rows(row_count).status = "ok";
            elseif ~cand.clock_sample_valid
                rows(row_count).status = cand.status;
            else
                rows(row_count).status = "receiver_clock_unavailable";
            end
        end
    end
end

if row_count == 0
    residual_tbl = empty_residual_table_local();
else
    residual_tbl = struct2table(rows(1:row_count));
end

if clock_count == 0
    clock_tbl = empty_clock_table_local();
else
    clock_tbl = struct2table(clock_rows(1:clock_count));
end
end

function candidates = build_epoch_candidates_local( ...
    epoch, epoch_time, epoch_idx, receiver_ecef, lat_deg, lon_deg, nav_data, cfg)
c = 299792458.0;
sat_ids = fieldnames(epoch.data);
candidates = repmat(empty_candidate_row_local(), max(1, numel(sat_ids)), 1);
candidate_count = 0;
theory_time = epoch_time + seconds(cfg.theory_time_shift_seconds);

for i = 1:numel(sat_ids)
    sat_id = sat_ids{i};
    sys_char = upper(sat_id(1));
    if ~system_enabled_local(sys_char, cfg.systems)
        continue;
    end

    [pseudorange, pr_code] = select_pseudorange_local(epoch.data.(sat_id), sys_char);
    if ~isfinite(pseudorange) || pseudorange <= 0
        continue;
    end

    candidate_count = candidate_count + 1;
    if candidate_count > numel(candidates)
        candidates(end + 1:end + 64, 1) = repmat(empty_candidate_row_local(), 64, 1); %#ok<AGROW>
    end

    candidates(candidate_count).epoch_idx = epoch_idx;
    candidates(candidate_count).time = epoch_time;
    candidates(candidate_count).sat_id = string(sat_id);
    candidates(candidate_count).system = sys_char;
    candidates(candidate_count).pseudorange_code = string(pr_code);
    candidates(candidate_count).pseudorange_raw_m = pseudorange;

    try
        gt = compute_gt_range_for_satellite_local( ...
            theory_time, receiver_ecef, lat_deg, lon_deg, sat_id, nav_data, cfg);
    catch ME
        candidates(candidate_count).status = "rho_gt_failed: " + string(ME.message);
        continue;
    end

    candidates(candidate_count).rho_gt_m = gt.range_m;
    candidates(candidate_count).azimuth_deg = gt.azimuth_deg;
    candidates(candidate_count).elevation_deg = gt.elevation_deg;
    candidates(candidate_count).visible = gt.elevation_deg >= cfg.elevation_mask_deg;

    try
        [~, ~, sat_clk_err] = calculate_satellite_state(theory_time, pseudorange, sat_id, nav_data);
        candidates(candidate_count).sat_clock_error_s = sat_clk_err;
        candidates(candidate_count).clock_sample_m = ...
            pseudorange - candidates(candidate_count).rho_gt_m + c * sat_clk_err;
        candidates(candidate_count).clock_sample_valid = ...
            isfinite(candidates(candidate_count).clock_sample_m);
        candidates(candidate_count).status = "ok";
    catch ME
        candidates(candidate_count).status = "sat_clock_failed: " + string(ME.message);
    end
end

candidates = candidates(1:candidate_count);
end

function value = estimate_clock_bias_local(samples_m, estimator)
switch lower(char(string(estimator)))
    case 'mean'
        value = mean(samples_m);
    case 'median'
        value = median(samples_m);
    otherwise
        error('calculate_nav_obs_range_residuals:InvalidClockEstimator', ...
            'receiver_clock_estimator must be ''mean'' or ''median''.');
end
end

function gt = compute_gt_range_for_satellite_local(t_utc, receiver_ecef, lat_deg, lon_deg, sat_id, nav_data, cfg)
c = 299792458.0;
rho = cfg.initial_range_m;
sat_rx = [NaN; NaN; NaN];
sat_clk_err = NaN;

for iter = 1:max(1, round(cfg.max_iterations))
    [sat_tx, ~, sat_clk_err] = calculate_satellite_state(t_utc, rho, sat_id, nav_data);
    transit_time_s = rho / c;

    if cfg.apply_sagnac
        sat_rx = rotate_satellite_to_receive_frame_local(sat_tx, transit_time_s);
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

gt = struct();
gt.range_m = rho;
gt.azimuth_deg = az_deg;
gt.elevation_deg = el_deg;
gt.sat_clock_error_s = sat_clk_err;
end

function sat_rx = rotate_satellite_to_receive_frame_local(sat_tx, transit_time_s)
omega_e = 7.2921151467e-5;
theta = omega_e * transit_time_s;
rot_z = [cos(theta), sin(theta), 0; ...
        -sin(theta), cos(theta), 0; ...
         0,          0,          1];
sat_rx = rot_z * sat_tx;
end

function epoch_indices = select_epoch_indices_local(obs_data, cfg)
all_times = to_utc_local([obs_data.time], cfg.time_zone);
mask = true(1, numel(obs_data));

if ~isempty(cfg.start_time)
    start_time = to_utc_local(cfg.start_time, cfg.time_zone);
    mask = mask & (all_times >= start_time);
end

if ~isempty(cfg.end_time)
    end_time = to_utc_local(cfg.end_time, cfg.time_zone);
    mask = mask & (all_times <= end_time);
elseif ~isempty(cfg.start_time) && isfinite(cfg.duration_seconds)
    end_time = to_utc_local(cfg.start_time, cfg.time_zone) + seconds(cfg.duration_seconds);
    mask = mask & (all_times <= end_time);
end

epoch_indices = find(mask);
if ~isempty(cfg.epoch_indices)
    epoch_indices = intersect(epoch_indices, cfg.epoch_indices(:).');
end

if isfinite(cfg.max_epochs) && numel(epoch_indices) > cfg.max_epochs
    epoch_indices = epoch_indices(1:cfg.max_epochs);
end
end

function cfg = default_cfg_local()
cfg = struct();
cfg.time_zone = 'UTC';
cfg.systems = {'G', 'C', 'E', 'R', 'J'};
cfg.elevation_mask_deg = 0;
cfg.min_clock_sats = 4;
cfg.receiver_clock_estimator = 'mean';
cfg.apply_sagnac = true;
cfg.max_iterations = 6;
cfg.range_tolerance_m = 1e-3;
cfg.initial_range_m = 26560000;
cfg.start_time = [];
cfg.end_time = [];
cfg.duration_seconds = NaN;
cfg.epoch_indices = [];
cfg.max_epochs = Inf;
cfg.theory_time_shift_seconds = 0;
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
    error('calculate_nav_obs_range_residuals:InvalidOptions', ...
        'Options must be name-value pairs or one scalar struct.');
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
    case 'minclocksats'
        cfg.min_clock_sats = double(value);
    case {'receiverclockestimator', 'clockestimator'}
        cfg.receiver_clock_estimator = char(string(value));
    case 'applysagnac'
        cfg.apply_sagnac = logical(value);
    case 'maxiterations'
        cfg.max_iterations = double(value);
    case 'rangetolerancem'
        cfg.range_tolerance_m = double(value);
    case 'initialrangem'
        cfg.initial_range_m = double(value);
    case 'starttime'
        cfg.start_time = value;
    case 'endtime'
        cfg.end_time = value;
    case {'durationseconds', 'windows'}
        cfg.duration_seconds = double(value);
    case 'epochindices'
        cfg.epoch_indices = double(value);
    case 'maxepochs'
        cfg.max_epochs = double(value);
    case {'theorytimeshiftseconds', 'timeshiftseconds'}
        cfg.theory_time_shift_seconds = double(value);
    otherwise
        error('calculate_nav_obs_range_residuals:UnknownOption', ...
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
    error('calculate_nav_obs_range_residuals:InvalidSystems', ...
        'systems must be a char, string array, or cell array.');
end
systems = systems(:).';
end

function nav_data = load_nav_data_local(nav_input)
if iscell(nav_input)
    nav_data = nav_input;
    return;
end

nav_path = char(string(nav_input));
if ~exist(nav_path, 'file')
    error('calculate_nav_obs_range_residuals:MissingNavFile', ...
        'Navigation file not found: %s', nav_path);
end
nav_data = parse_rinex_nav_multi_gnss(nav_path);
end

function obs_data = load_obs_data_local(obs_input)
if isstruct(obs_input)
    obs_data = obs_input;
    return;
end

obs_path = char(string(obs_input));
if ~exist(obs_path, 'file')
    error('calculate_nav_obs_range_residuals:MissingObsFile', ...
        'Observation file not found: %s', obs_path);
end
obs_data = parse_rinex_obs(obs_path);
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

function tf = system_enabled_local(sys_char, systems)
systems = cellstr(upper(string(systems)));
tf = any(strcmpi(sys_char, systems));
end

function t = to_utc_local(t, default_time_zone)
if isa(t, 'datetime')
    if isempty(t.TimeZone)
        t.TimeZone = default_time_zone;
    end
    t.TimeZone = 'UTC';
    return;
end

if isnumeric(t) && isscalar(t)
    t = datetime(t, 'ConvertFrom', 'datenum', 'TimeZone', default_time_zone);
elseif isnumeric(t) && numel(t) >= 6
    vec = double(t(:).');
    t = datetime(vec(1), vec(2), vec(3), vec(4), vec(5), vec(6), ...
        'TimeZone', default_time_zone);
else
    t = datetime(char(string(t)), 'TimeZone', default_time_zone);
end
t.TimeZone = 'UTC';
end

function row = empty_candidate_row_local()
row = struct( ...
    'epoch_idx', NaN, ...
    'time', NaT('TimeZone', 'UTC'), ...
    'sat_id', "", ...
    'system', "", ...
    'pseudorange_code', "", ...
    'pseudorange_raw_m', NaN, ...
    'rho_gt_m', NaN, ...
    'azimuth_deg', NaN, ...
    'elevation_deg', NaN, ...
    'visible', false, ...
    'sat_clock_error_s', NaN, ...
    'clock_sample_m', NaN, ...
    'clock_sample_valid', false, ...
    'status', "");
end

function row = empty_residual_row_local()
row = struct( ...
    'epoch_idx', NaN, ...
    'time', NaT('TimeZone', 'UTC'), ...
    'sat_id', "", ...
    'system', "", ...
    'pseudorange_code', "", ...
    'pseudorange_raw_m', NaN, ...
    'rho_gt_m', NaN, ...
    'r_obs_corrected_m', NaN, ...
    'residual_obs_minus_gt_m', NaN, ...
    'receiver_clock_bias_s', NaN, ...
    'receiver_clock_bias_m', NaN, ...
    'sat_clock_error_s', NaN, ...
    'sat_clock_correction_m', NaN, ...
    'clock_sample_m', NaN, ...
    'clock_sample_valid', false, ...
    'azimuth_deg', NaN, ...
    'elevation_deg', NaN, ...
    'visible', false, ...
    'status', "");
end

function row = empty_clock_row_local()
row = struct( ...
    'epoch_idx', NaN, ...
    'time', NaT('TimeZone', 'UTC'), ...
    'system', "", ...
    'candidate_count', NaN, ...
    'used_count', NaN, ...
    'receiver_clock_bias_s', NaN, ...
    'receiver_clock_bias_m', NaN, ...
    'clock_sample_std_m', NaN, ...
    'status', "");
end

function tbl = empty_residual_table_local()
tbl = struct2table(repmat(empty_residual_row_local(), 0, 1));
end

function tbl = empty_clock_table_local()
tbl = struct2table(repmat(empty_clock_row_local(), 0, 1));
end

function ensure_project_paths_local(repo_dir)
addpath(fullfile(repo_dir, 'nav_parse'));
addpath(fullfile(repo_dir, 'obs_parse'));
addpath(fullfile(repo_dir, 'calculate_clock_bias_and_positon'));
addpath(fullfile(repo_dir, 'physical_consistency'));
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
        error('calculate_nav_obs_range_residuals:RepoRootNotFound', ...
            'Cannot locate SatLock repo root from %s.', start_dir);
    end
    repo_dir = parent_dir;
end
end
