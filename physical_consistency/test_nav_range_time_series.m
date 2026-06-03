% Test pure NAV-predicted satellite range changes over a short time window.
%
% This script uses only:
%   1. RINEX navigation data
%   2. A start time
%   3. A fixed geodetic point: latitude, longitude, height
%
% It does not use OBS, pseudorange, carrier phase, SNR, or receiver data.

clear;
clc;
close all;

repo_dir = resolve_repo_root_local(fileparts(mfilename('fullpath')));
addpath(genpath(repo_dir));

%% ===================== Edit here =====================
nav_file = fullfile(repo_dir, 'data', '1_8', '2026_1_8.nav');

start_time = datetime(2026, 1, 8, 5, 59, 39.116, 'TimeZone', 'UTC');
duration_minutes = 1;
time_step_seconds = 5;

lat_deg = 35.77204197;
lon_deg = 120.02477211;
height_m = 92.889;

systems = {'G', 'C', 'E', 'R', 'J'};
elevation_mask_deg = 0;

% false: plot every satellite with a valid NAV-computed range.
% true : plot only satellites visible at least once in this time window.
plot_visible_only = false;

% One satellite per subplot. Satellites are split across pages/figures.
satellites_per_figure = 24;
subplot_columns = 4;

% Plot only the small change relative to each satellite's first valid range.
% Use 'm' for meter-level change, or 'km' for kilometer-level change.
range_change_unit = 'm';

summary_top_n = 20;

save_plot = false;
output_png = fullfile(repo_dir, 'physical_consistency', 'outputs', ...
    'nav_range_time_series.png');
%% =====================================================

fprintf('NAV file          : %s\n', nav_file);
fprintf('Start time UTC    : %s\n', datestr(start_time, 'yyyy-mm-dd HH:MM:SS.FFF'));
fprintf('Duration          : %.3f min\n', duration_minutes);
fprintf('Time step         : %.3f s\n', time_step_seconds);
fprintf('Location          : %.8f deg, %.8f deg, %.3f m\n', lat_deg, lon_deg, height_m);

time_offsets_s = 0:time_step_seconds:(duration_minutes * 60);
time_vec = start_time + seconds(time_offsets_s);

fprintf('\nParsing NAV once...\n');
nav_data = parse_rinex_nav_multi_gnss(nav_file);

range_tables = cell(numel(time_vec), 1);
all_sat_ids = string.empty(0, 1);

fprintf('Computing pure NAV ranges for %d epochs...\n', numel(time_vec));
for k = 1:numel(time_vec)
    range_tbl = calculate_satellite_ranges_from_nav( ...
        nav_data, time_vec(k), lat_deg, lon_deg, height_m, ...
        'systems', systems, ...
        'elevation_mask_deg', elevation_mask_deg, ...
        'visible_only', false);

    range_tables{k} = range_tbl;
    all_sat_ids = union(all_sat_ids, range_tbl.sat_id);

    fprintf('  %2d/%2d  %s UTC  satellites: %d, visible: %d\n', ...
        k, numel(time_vec), datestr(time_vec(k), 'HH:MM:SS.FFF'), ...
        height(range_tbl), sum(range_tbl.visible));
end

all_sat_ids = sort_sat_ids_local(all_sat_ids);
range_km = NaN(numel(time_vec), numel(all_sat_ids));
visible_mat = false(numel(time_vec), numel(all_sat_ids));

for k = 1:numel(time_vec)
    range_tbl = range_tables{k};
    [tf, loc] = ismember(range_tbl.sat_id, all_sat_ids);
    loc = loc(tf);
    range_km(k, loc) = range_tbl.range_km(tf);
    visible_mat(k, loc) = range_tbl.visible(tf);
end

valid_count = sum(isfinite(range_km), 1);
plot_mask = valid_count >= 2;
if plot_visible_only
    plot_mask = plot_mask & any(visible_mat, 1);
end

plot_sat_ids = all_sat_ids(plot_mask);
plot_range_km = range_km(:, plot_mask);
plot_visible_mat = visible_mat(:, plot_mask);

if isempty(plot_sat_ids)
    error('test_nav_range_time_series:NoSatellitesToPlot', ...
        'No satellites have enough valid range samples to plot.');
end

baseline_km = first_valid_value_local(plot_range_km);
[plot_range_change, y_unit] = range_change_from_baseline_local( ...
    plot_range_km, baseline_km, range_change_unit);
range_change_summary = build_range_change_summary_local( ...
    plot_sat_ids, plot_range_km, plot_visible_mat);
[~, summary_order] = sort(range_change_summary.peak_to_peak_change_m, 'descend');
range_change_summary = range_change_summary(summary_order, :);

fprintf('\nTotal satellites with valid ranges : %d\n', numel(all_sat_ids));
fprintf('Satellites plotted                : %d\n', numel(plot_sat_ids));
fprintf('Plot mode                         : %s\n', plot_mode_text_local(plot_visible_only));
fprintf('Range offset                       : first valid range of each satellite\n');
fprintf('Range change unit                  : %s\n', y_unit);

fprintf('\nTop %d satellites by %.1f-minute peak-to-peak range change:\n', ...
    min(summary_top_n, height(range_change_summary)), duration_minutes);
disp(range_change_summary(1:min(summary_top_n, height(range_change_summary)), ...
    {'sat_id', 'baseline_km', 'end_change_m', 'min_change_m', ...
    'max_change_m', 'peak_to_peak_change_m', 'visible_any'}));

minute_axis = minutes(time_vec - start_time);
satellites_per_figure = max(1, round(satellites_per_figure));
subplot_columns = max(1, round(subplot_columns));
page_count = ceil(numel(plot_sat_ids) / satellites_per_figure);
fig_handles = gobjects(page_count, 1);

fprintf('Subplot pages                      : %d\n', page_count);

for page_idx = 1:page_count
    sat_start = (page_idx - 1) * satellites_per_figure + 1;
    sat_end = min(page_idx * satellites_per_figure, numel(plot_sat_ids));
    sat_idx = sat_start:sat_end;
    n_tiles = numel(sat_idx);
    n_cols = min(subplot_columns, n_tiles);
    n_rows = ceil(n_tiles / n_cols);

    fig_h = figure('Name', sprintf('Pure NAV range changes page %d', page_idx), ...
        'Color', 'w', ...
        'Position', [60, 40, 1600, 920]);
    fig_handles(page_idx) = fig_h;

    tlo = tiledlayout(fig_h, n_rows, n_cols, ...
        'TileSpacing', 'compact', ...
        'Padding', 'compact');

    title(tlo, sprintf(['Pure NAV satellite-to-point range changes | %s UTC ' ...
        'to +%.1f min | page %d/%d'], ...
        datestr(start_time, 'yyyy-mm-dd HH:MM:SS.FFF'), duration_minutes, ...
        page_idx, page_count));
    subtitle(tlo, sprintf(['Point: %.8f deg, %.8f deg, %.3f m | ' ...
        'each subplot subtracts its first valid range'], ...
        lat_deg, lon_deg, height_m));
    xlabel(tlo, 'Minutes after start time');
    ylabel(tlo, sprintf('Range change (%s)', y_unit));

    for local_i = 1:n_tiles
        col_idx = sat_idx(local_i);
        sat_id = char(plot_sat_ids(col_idx));
        col = system_color_local(sat_id(1));

        ax = nexttile(tlo);
        plot(ax, minute_axis, plot_range_change(:, col_idx), ...
            'LineWidth', 1.1, ...
            'Color', col);
        hold(ax, 'on');
        yline(ax, 0, ':', 'Color', [0.45 0.45 0.45], 'LineWidth', 0.7);
        hold(ax, 'off');

        grid(ax, 'on');
        box(ax, 'on');
        xlim(ax, [minute_axis(1), minute_axis(end)]);
        ylim(ax, padded_ylim_local(plot_range_change(:, col_idx)));
        title(ax, sprintf('%s | offset %.3f km', sat_id, baseline_km(col_idx)), ...
            'FontSize', 8, ...
            'FontWeight', 'bold', ...
            'Color', col);
        ax.FontSize = 7;
    end

    if save_plot
        out_path = page_output_path_local(output_png, page_idx, page_count);
        out_dir = fileparts(out_path);
        if ~isempty(out_dir) && ~exist(out_dir, 'dir')
            mkdir(out_dir);
        end
        exportgraphics(fig_h, out_path, 'Resolution', 300);
        fprintf('Saved plot page %d: %s\n', page_idx, out_path);
    end
end

%% ===================== local helpers =====================
function baseline_km = first_valid_value_local(range_km)
baseline_km = NaN(1, size(range_km, 2));
for i = 1:size(range_km, 2)
    first_idx = find(isfinite(range_km(:, i)), 1, 'first');
    if ~isempty(first_idx)
        baseline_km(i) = range_km(first_idx, i);
    end
end
end

function [range_change, y_unit] = range_change_from_baseline_local(range_km, baseline_km, unit_name)
unit_key = lower(char(string(unit_name)));
range_change = range_km - baseline_km;

switch unit_key
    case {'m', 'meter', 'meters'}
        range_change = range_change * 1000;
        y_unit = 'm';
    case {'km', 'kilometer', 'kilometers'}
        y_unit = 'km';
    otherwise
        error('test_nav_range_time_series:InvalidRangeChangeUnit', ...
            'range_change_unit must be ''m'' or ''km''.');
end
end

function summary_tbl = build_range_change_summary_local(sat_ids, range_km, visible_mat)
sat_ids = string(sat_ids(:));
n_sat = numel(sat_ids);

baseline_km = NaN(n_sat, 1);
end_change_m = NaN(n_sat, 1);
min_change_m = NaN(n_sat, 1);
max_change_m = NaN(n_sat, 1);
peak_to_peak_change_m = NaN(n_sat, 1);
visible_any = false(n_sat, 1);

for i = 1:n_sat
    valid_idx = find(isfinite(range_km(:, i)));
    if isempty(valid_idx)
        continue;
    end

    baseline_km(i) = range_km(valid_idx(1), i);
    change_m = (range_km(valid_idx, i) - baseline_km(i)) * 1000;
    end_change_m(i) = change_m(end);
    min_change_m(i) = min(change_m);
    max_change_m(i) = max(change_m);
    peak_to_peak_change_m(i) = max_change_m(i) - min_change_m(i);
    visible_any(i) = any(visible_mat(:, i));
end

summary_tbl = table(sat_ids, baseline_km, end_change_m, min_change_m, ...
    max_change_m, peak_to_peak_change_m, visible_any, ...
    'VariableNames', {'sat_id', 'baseline_km', 'end_change_m', ...
    'min_change_m', 'max_change_m', 'peak_to_peak_change_m', ...
    'visible_any'});
end

function limits = padded_ylim_local(values)
values = values(isfinite(values));
if isempty(values)
    limits = [-1, 1];
    return;
end

v_min = min(values);
v_max = max(values);
span = v_max - v_min;
if span <= 0
    pad = max(1, 0.05 * max(1, abs(v_min)));
else
    pad = 0.12 * span;
end
limits = [v_min - pad, v_max + pad];
end

function out_path = page_output_path_local(output_png, page_idx, page_count)
[out_dir, base_name, ext] = fileparts(output_png);
if isempty(ext)
    ext = '.png';
end

if page_count <= 1
    out_name = [base_name ext];
else
    out_name = sprintf('%s_page_%02d%s', base_name, page_idx, ext);
end
out_path = fullfile(out_dir, out_name);
end

function sat_ids = sort_sat_ids_local(sat_ids)
sat_ids = string(sat_ids(:));
if isempty(sat_ids)
    return;
end

sys_order = containers.Map({'G', 'C', 'E', 'R', 'J'}, [1, 2, 3, 4, 5]);
sort_key = zeros(numel(sat_ids), 2);
for i = 1:numel(sat_ids)
    sat_id = char(sat_ids(i));
    sys_char = upper(sat_id(1));
    if isKey(sys_order, sys_char)
        sort_key(i, 1) = sys_order(sys_char);
    else
        sort_key(i, 1) = 99;
    end
    sort_key(i, 2) = str2double(sat_id(2:end));
end

[~, idx] = sortrows(sort_key, [1, 2]);
sat_ids = sat_ids(idx);
end

function txt = plot_mode_text_local(plot_visible_only)
if plot_visible_only
    txt = 'visible satellites only';
else
    txt = 'all satellites with valid NAV ranges';
end
end

function add_system_legend_local(pax, sat_ids)
seen = {};
handles = gobjects(0);
labels = {};

for i = 1:numel(sat_ids)
    sat_id = char(sat_ids(i));
    sys_char = upper(sat_id(1));
    if any(strcmp(seen, sys_char))
        continue;
    end

    seen{end + 1} = sys_char; %#ok<AGROW>
    h = plot(pax, NaN, NaN, '-', ...
        'Color', system_color_local(sys_char), ...
        'LineWidth', 1.8);
    handles(end + 1) = h; %#ok<AGROW>
    labels{end + 1} = system_label_local(sys_char); %#ok<AGROW>
end

if ~isempty(handles)
    legend(pax, handles, labels, 'Location', 'eastoutside');
end
end

function label = system_label_local(sys_char)
switch upper(char(sys_char))
    case 'G'
        label = 'GPS';
    case 'C'
        label = 'BeiDou';
    case 'E'
        label = 'Galileo';
    case 'R'
        label = 'GLONASS';
    case 'J'
        label = 'QZSS';
    otherwise
        label = upper(char(sys_char));
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
        error('test_nav_range_time_series:RepoRootNotFound', ...
            'Cannot locate SatLock repo root from %s.', start_dir);
    end
    repo_dir = parent_dir;
end
end
