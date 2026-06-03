% Example: compare NAV-predicted sky plots at one location over time.
%
% This script uses only broadcast NAV data and a fixed geodetic point. It
% samples the sky plot from start_time to start_time + compare_duration_hours,
% then draws one overlay sky plot showing satellite motion over that interval.

clear;
clc;
close all;

repo_dir = resolve_repo_root_local(fileparts(mfilename('fullpath')));
addpath(genpath(repo_dir));

%% ===================== Edit here =====================
nav_file = fullfile(repo_dir, 'data', '1_8', '2026_1_8.nav');

start_time = datetime(2026, 1, 8, 5, 59, 39.116, 'TimeZone', 'UTC');
compare_duration_hours = 0.5;
time_step_minutes = 15;

lat_deg = 35.77204196;
lon_deg = 120.02477210;
height_m = 92.889;

systems = {'G', 'C', 'E', 'R', 'J'};
elevation_mask_deg = 0;

show_satellite_labels = false;
show_title = false;

out_dir = fullfile(repo_dir, 'physical_consistency', 'outputs', ...
    'skyplot_time_comparison');
out_pdf_dir = fullfile(out_dir, 'pdf');
out_png_dir = fullfile(out_dir, 'png');
out_fig_dir = fullfile(out_dir, 'fig');
out_csv_dir = fullfile(out_dir, 'csv');

paper_font_name = 'Arial';
paper_tick_font_size = 18;
paper_label_font_size = 20;
paper_legend_font_size = 15;
paper_note_font_size = 14;
paper_axis_line_width = 1.1;
%% =====================================================

output_dirs = {out_dir, out_pdf_dir, out_png_dir, out_fig_dir, out_csv_dir};
for di = 1:numel(output_dirs)
    if ~exist(output_dirs{di}, 'dir')
        mkdir(output_dirs{di});
    end
end
organize_existing_outputs_local(out_dir, out_png_dir, out_pdf_dir, ...
    out_fig_dir, out_csv_dir);

duration_minutes = compare_duration_hours * 60;
time_offsets_min = 0:time_step_minutes:duration_minutes;
if time_offsets_min(end) < duration_minutes
    time_offsets_min(end + 1) = duration_minutes; %#ok<SAGROW>
end
time_vec = start_time + minutes(time_offsets_min);

duration_tag = duration_tag_local(compare_duration_hours);
figure_base = sprintf('figure_skyplot_time_comparison_%s', duration_tag);
png_path = fullfile(out_png_dir, [figure_base '.png']);
pdf_path = fullfile(out_pdf_dir, [figure_base '.pdf']);
fig_path = fullfile(out_fig_dir, [figure_base '.fig']);
all_epochs_csv = fullfile(out_csv_dir, ...
    sprintf('skyplot_time_comparison_%s_all_epochs.csv', duration_tag));
endpoint_csv = fullfile(out_csv_dir, ...
    sprintf('skyplot_time_comparison_%s_endpoint_delta.csv', duration_tag));
summary_csv = fullfile(out_csv_dir, ...
    sprintf('skyplot_time_comparison_%s_summary.csv', duration_tag));

fprintf('NAV file : %s\n', nav_file);
fprintf('Location : %.8f deg, %.8f deg, %.3f m\n', ...
    lat_deg, lon_deg, height_m);
fprintf('Window   : %s UTC to %s UTC (%.2f h)\n', ...
    datestr(start_time, 'yyyy-mm-dd HH:MM:SS.FFF'), ...
    datestr(time_vec(end), 'yyyy-mm-dd HH:MM:SS.FFF'), ...
    compare_duration_hours);

fprintf('\nParsing NAV once...\n');
nav_data = parse_rinex_nav_multi_gnss(nav_file);

fprintf('Computing sky plots for %d epochs...\n', numel(time_vec));
all_sky_tbl = table();
for ti = 1:numel(time_vec)
    [sky_tbl, ~] = calculate_and_plot_skyplot_from_nav( ...
        nav_data, time_vec(ti), lat_deg, lon_deg, height_m, ...
        'systems', systems, ...
        'elevation_mask_deg', elevation_mask_deg, ...
        'plot_enable', false);

    sky_tbl.epoch_time = repmat(time_vec(ti), height(sky_tbl), 1);
    sky_tbl.time_offset_min = repmat(time_offsets_min(ti), height(sky_tbl), 1);
    sky_tbl = movevars(sky_tbl, {'epoch_time', 'time_offset_min'}, 'Before', 1);
    all_sky_tbl = [all_sky_tbl; sky_tbl]; %#ok<AGROW>

    visible_count = sum(sky_tbl.visible & isfinite(sky_tbl.azimuth_deg) & ...
        isfinite(sky_tbl.elevation_deg));
    fprintf('  %+7.1f min -> %d visible satellites\n', ...
        time_offsets_min(ti), visible_count);
end

visible_tbl = all_sky_tbl(all_sky_tbl.visible & ...
    isfinite(all_sky_tbl.azimuth_deg) & ...
    isfinite(all_sky_tbl.elevation_deg), :);

start_tbl = visible_tbl(visible_tbl.time_offset_min == time_offsets_min(1), :);
end_tbl = visible_tbl(visible_tbl.time_offset_min == time_offsets_min(end), :);
endpoint_tbl = build_endpoint_comparison_local(start_tbl, end_tbl);
summary_tbl = build_summary_table_local(endpoint_tbl, start_time, time_vec(end), ...
    compare_duration_hours, time_step_minutes, lat_deg, lon_deg, height_m, ...
    elevation_mask_deg);

writetable(all_sky_tbl, all_epochs_csv);
writetable(endpoint_tbl, endpoint_csv);
writetable(summary_tbl, summary_csv);

fig_h = plot_time_comparison_skyplot_local(visible_tbl, endpoint_tbl, ...
    time_vec, time_offsets_min, lat_deg, lon_deg, height_m, ...
    elevation_mask_deg, compare_duration_hours, show_satellite_labels, ...
    show_title, paper_font_name, paper_tick_font_size, paper_label_font_size, ...
    paper_legend_font_size, paper_note_font_size, paper_axis_line_width);

save_paper_figure_local(fig_h, png_path, pdf_path, fig_path);

fprintf('\nSaved skyplot time-comparison outputs:\n');
fprintf('  %s\n', png_path);
fprintf('  %s\n', pdf_path);
fprintf('  %s\n', fig_path);
fprintf('  %s\n', endpoint_csv);
fprintf('  %s\n', summary_csv);

common_mask = endpoint_tbl.visibility_class == "common";
fprintf('\nEndpoint comparison:\n');
fprintf('  t0 visible      : %d\n', height(start_tbl));
fprintf('  t0+%.1fh visible: %d\n', compare_duration_hours, height(end_tbl));
fprintf('  common visible  : %d\n', sum(common_mask));
if any(common_mask)
    fprintf('  median sky angular shift: %.2f deg\n', ...
        median(endpoint_tbl.angular_shift_deg(common_mask), 'omitnan'));
    fprintf('  max sky angular shift   : %.2f deg\n', ...
        max(endpoint_tbl.angular_shift_deg(common_mask), [], 'omitnan'));
end

%% ===================== local helpers =====================
function endpoint_tbl = build_endpoint_comparison_local(start_tbl, end_tbl)
start_ids = string(start_tbl.sat_id);
end_ids = string(end_tbl.sat_id);
all_ids = unique([start_ids; end_ids], 'stable');
n = numel(all_ids);

sat_id = strings(n, 1);
system = strings(n, 1);
prn = NaN(n, 1);
azimuth_start_deg = NaN(n, 1);
elevation_start_deg = NaN(n, 1);
azimuth_end_deg = NaN(n, 1);
elevation_end_deg = NaN(n, 1);
delta_azimuth_deg = NaN(n, 1);
delta_elevation_deg = NaN(n, 1);
angular_shift_deg = NaN(n, 1);
visibility_class = strings(n, 1);

for i = 1:n
    id = all_ids(i);
    sat_id(i) = id;
    si = find(start_ids == id, 1);
    ei = find(end_ids == id, 1);

    if ~isempty(si)
        system(i) = string(start_tbl.system(si));
        prn(i) = start_tbl.prn(si);
        azimuth_start_deg(i) = start_tbl.azimuth_deg(si);
        elevation_start_deg(i) = start_tbl.elevation_deg(si);
    end

    if ~isempty(ei)
        if strlength(system(i)) == 0
            system(i) = string(end_tbl.system(ei));
            prn(i) = end_tbl.prn(ei);
        end
        azimuth_end_deg(i) = end_tbl.azimuth_deg(ei);
        elevation_end_deg(i) = end_tbl.elevation_deg(ei);
    end

    if ~isempty(si) && ~isempty(ei)
        visibility_class(i) = "common";
        delta_azimuth_deg(i) = wrap_to_180_local( ...
            azimuth_end_deg(i) - azimuth_start_deg(i));
        delta_elevation_deg(i) = elevation_end_deg(i) - elevation_start_deg(i);
        angular_shift_deg(i) = sky_angular_separation_deg_local( ...
            azimuth_start_deg(i), elevation_start_deg(i), ...
            azimuth_end_deg(i), elevation_end_deg(i));
    elseif ~isempty(si)
        visibility_class(i) = "start_only";
    else
        visibility_class(i) = "end_only";
    end
end

endpoint_tbl = table(sat_id, system, prn, ...
    azimuth_start_deg, elevation_start_deg, ...
    azimuth_end_deg, elevation_end_deg, ...
    delta_azimuth_deg, delta_elevation_deg, angular_shift_deg, ...
    visibility_class);
endpoint_tbl = sortrows(endpoint_tbl, {'visibility_class', 'sat_id'});
end

function summary_tbl = build_summary_table_local(endpoint_tbl, start_time, ...
    end_time, duration_hours, time_step_minutes, lat_deg, lon_deg, height_m, ...
    elevation_mask_deg)
common_mask = endpoint_tbl.visibility_class == "common";
start_only_mask = endpoint_tbl.visibility_class == "start_only";
end_only_mask = endpoint_tbl.visibility_class == "end_only";

common_shift = endpoint_tbl.angular_shift_deg(common_mask);
summary_tbl = table( ...
    start_time, end_time, duration_hours, time_step_minutes, ...
    lat_deg, lon_deg, height_m, elevation_mask_deg, ...
    sum(common_mask), sum(start_only_mask), sum(end_only_mask), ...
    median(common_shift, 'omitnan'), ...
    mean(common_shift, 'omitnan'), ...
    max(common_shift, [], 'omitnan'), ...
    'VariableNames', {'start_time', 'end_time', 'duration_hours', ...
    'time_step_minutes', 'lat_deg', 'lon_deg', 'height_m', ...
    'elevation_mask_deg', 'common_visible_count', ...
    'start_only_count', 'end_only_count', ...
    'median_angular_shift_deg', 'mean_angular_shift_deg', ...
    'max_angular_shift_deg'});
end

function fig_h = plot_time_comparison_skyplot_local(visible_tbl, endpoint_tbl, ...
    time_vec, time_offsets_min, lat_deg, lon_deg, height_m, elevation_mask_deg, ...
    duration_hours, show_satellite_labels, show_title, font_name, ...
    tick_font_size, label_font_size, legend_font_size, note_font_size, ...
    axis_line_width)
fig_h = figure('Name', 'Sky plot time comparison', ...
    'Color', 'w', 'Position', [120, 90, 1280, 980]);
pax = polaraxes(fig_h);
set(pax, ...
    'ThetaZeroLocation', 'top', ...
    'ThetaDir', 'clockwise', ...
    'RAxisLocation', 90, ...
    'FontName', font_name, ...
    'FontSize', tick_font_size, ...
    'LineWidth', axis_line_width);
rlim(pax, [0 90]);
rticks(pax, [0 30 60 90]);
rticklabels(pax, {'90', '60', '30', '0'});
pax.RAxis.Label.String = '';
pax.RAxis.Label.FontName = font_name;
pax.RAxis.Label.FontSize = label_font_size;
grid(pax, 'on');
hold(pax, 'on');

if show_title
    title(pax, sprintf('Sky plot comparison over %.1f h', duration_hours), ...
        'FontName', font_name, 'FontSize', label_font_size);
end

track_color = [0.62 0.62 0.62];
start_color = [0.19 0.40 0.67];
end_color = [0.85 0.33 0.10];
start_only_color = [0.00 0.45 0.74];
end_only_color = [0.93 0.69 0.13];
label_color = [0.10 0.10 0.10];

sat_ids = unique(string(visible_tbl.sat_id), 'stable');
for i = 1:numel(sat_ids)
    sat_mask = string(visible_tbl.sat_id) == sat_ids(i);
    sat_tbl = sortrows(visible_tbl(sat_mask, :), 'time_offset_min');
    if height(sat_tbl) < 2
        continue;
    end

    theta = deg2rad(unwrap_degrees_local(sat_tbl.azimuth_deg));
    radius = 90 - sat_tbl.elevation_deg;
    polarplot(pax, theta, radius, '-', ...
        'Color', track_color, 'LineWidth', 1.0, ...
        'HandleVisibility', 'off');
end

common_mask = endpoint_tbl.visibility_class == "common";
start_only_mask = endpoint_tbl.visibility_class == "start_only";
end_only_mask = endpoint_tbl.visibility_class == "end_only";

legend_handles = gobjects(0);
legend_labels = {};

h_track = polarplot(pax, NaN, NaN, '-', ...
    'Color', track_color, 'LineWidth', 1.4);
legend_handles(end + 1) = h_track; %#ok<AGROW>
legend_labels{end + 1} = 'satellite trajectory'; %#ok<AGROW>

if any(common_mask | start_only_mask)
    start_rows = endpoint_tbl(common_mask | start_only_mask, :);
    h_start = polarplot(pax, deg2rad(start_rows.azimuth_start_deg), ...
        90 - start_rows.elevation_start_deg, 'o', ...
        'MarkerSize', 8, ...
        'MarkerFaceColor', 'w', ...
        'MarkerEdgeColor', start_color, ...
        'LineWidth', 1.8);
    legend_handles(end + 1) = h_start; %#ok<AGROW>
    legend_labels{end + 1} = sprintf('t0: %s UTC', ...
        datestr(time_vec(1), 'HH:MM')); %#ok<AGROW>
end

if any(common_mask | end_only_mask)
    end_rows = endpoint_tbl(common_mask | end_only_mask, :);
    h_end = polarplot(pax, deg2rad(end_rows.azimuth_end_deg), ...
        90 - end_rows.elevation_end_deg, 's', ...
        'MarkerSize', 8, ...
        'MarkerFaceColor', end_color, ...
        'MarkerEdgeColor', [0.15 0.15 0.15], ...
        'LineWidth', 1.0);
    legend_handles(end + 1) = h_end; %#ok<AGROW>
    legend_labels{end + 1} = sprintf('t0 + %.1f h: %s UTC', ...
        duration_hours, datestr(time_vec(end), 'HH:MM')); %#ok<AGROW>
end

if any(start_only_mask)
    rows = endpoint_tbl(start_only_mask, :);
    h_start_only = polarplot(pax, deg2rad(rows.azimuth_start_deg), ...
        90 - rows.elevation_start_deg, 'x', ...
        'MarkerSize', 9, ...
        'MarkerEdgeColor', start_only_color, ...
        'LineWidth', 1.8);
    legend_handles(end + 1) = h_start_only; %#ok<AGROW>
    legend_labels{end + 1} = 'visible only at t0'; %#ok<AGROW>
end

if any(end_only_mask)
    rows = endpoint_tbl(end_only_mask, :);
    h_end_only = polarplot(pax, deg2rad(rows.azimuth_end_deg), ...
        90 - rows.elevation_end_deg, 'd', ...
        'MarkerSize', 8, ...
        'MarkerFaceColor', end_only_color, ...
        'MarkerEdgeColor', [0.15 0.15 0.15], ...
        'LineWidth', 1.0);
    legend_handles(end + 1) = h_end_only; %#ok<AGROW>
    legend_labels{end + 1} = 'visible only at t0+duration'; %#ok<AGROW>
end

if show_satellite_labels
    label_endpoint_satellites_local(pax, endpoint_tbl(common_mask, :), ...
        "end", label_color);
    label_endpoint_satellites_local(pax, endpoint_tbl(start_only_mask, :), ...
        "start", start_only_color);
    label_endpoint_satellites_local(pax, endpoint_tbl(end_only_mask, :), ...
        "end", end_only_color);
end

lgd = legend(pax, legend_handles, legend_labels, ...
    'Location', 'eastoutside');
set(lgd, 'FontName', font_name, 'FontSize', legend_font_size, 'Box', 'on');

summary_note = skyplot_summary_note_local(endpoint_tbl, time_vec, ...
    time_offsets_min, lat_deg, lon_deg, height_m, elevation_mask_deg);
annotation(fig_h, 'textbox', [0.68 0.09 0.28 0.20], ...
    'String', summary_note, ...
    'FontName', font_name, ...
    'FontSize', note_font_size, ...
    'Interpreter', 'tex', ...
    'FitBoxToText', 'on', ...
    'BackgroundColor', 'white', ...
    'EdgeColor', [0.25 0.25 0.25], ...
    'Margin', 7);

set(pax, 'Units', 'normalized', 'Position', [0.08 0.11 0.64 0.78]);
hold(pax, 'off');
end

function label_endpoint_satellites_local(pax, rows, endpoint_name, color_rgb)
if isempty(rows)
    return;
end

for i = 1:height(rows)
    if endpoint_name == "start"
        az_deg = rows.azimuth_start_deg(i);
        el_deg = rows.elevation_start_deg(i);
    else
        az_deg = rows.azimuth_end_deg(i);
        el_deg = rows.elevation_end_deg(i);
    end

    if ~isfinite(az_deg) || ~isfinite(el_deg)
        continue;
    end

    radius = min(max(90 - el_deg + 2.0, 0), 90);
    text(pax, deg2rad(az_deg), radius, [' ' char(rows.sat_id(i))], ...
        'FontName', 'Arial', ...
        'FontSize', 7, ...
        'FontWeight', 'bold', ...
        'Color', color_rgb, ...
        'HorizontalAlignment', 'left', ...
        'VerticalAlignment', 'middle');
end
end

function note = skyplot_summary_note_local(endpoint_tbl, time_vec, ...
    time_offsets_min, lat_deg, lon_deg, height_m, elevation_mask_deg)
common_mask = endpoint_tbl.visibility_class == "common";
start_only_mask = endpoint_tbl.visibility_class == "start_only";
end_only_mask = endpoint_tbl.visibility_class == "end_only";
shift = endpoint_tbl.angular_shift_deg(common_mask);

if isempty(shift)
    median_shift = NaN;
    max_shift = NaN;
else
    median_shift = median(shift, 'omitnan');
    max_shift = max(shift, [], 'omitnan');
end

note = sprintf(['Location: %.6f, %.6f, %.1f m\n', ...
    'Time: %s to %s UTC\n', ...
    'Samples: %.0f-%.0f min, mask %.0f deg\n', ...
    'Common/start-only/end-only: %d/%d/%d\n', ...
    'Median/max angular shift: %.1f/%.1f deg'], ...
    lat_deg, lon_deg, height_m, ...
    datestr(time_vec(1), 'HH:MM'), datestr(time_vec(end), 'HH:MM'), ...
    min(time_offsets_min), max(time_offsets_min), elevation_mask_deg, ...
    sum(common_mask), sum(start_only_mask), sum(end_only_mask), ...
    median_shift, max_shift);
end

function theta_unwrapped = unwrap_degrees_local(theta_deg)
theta_unwrapped = theta_deg(:);
for i = 2:numel(theta_unwrapped)
    step = theta_unwrapped(i) - theta_unwrapped(i - 1);
    if step > 180
        theta_unwrapped(i:end) = theta_unwrapped(i:end) - 360;
    elseif step < -180
        theta_unwrapped(i:end) = theta_unwrapped(i:end) + 360;
    end
end
end

function d = wrap_to_180_local(x)
d = mod(x + 180, 360) - 180;
end

function angle_deg = sky_angular_separation_deg_local(az1_deg, el1_deg, ...
    az2_deg, el2_deg)
v1 = sky_unit_vector_local(az1_deg, el1_deg);
v2 = sky_unit_vector_local(az2_deg, el2_deg);
dot_value = max(-1, min(1, dot(v1, v2)));
angle_deg = acosd(dot_value);
end

function v = sky_unit_vector_local(az_deg, el_deg)
az = deg2rad(az_deg);
el = deg2rad(el_deg);
v = [cos(el) * sin(az); cos(el) * cos(az); sin(el)];
end

function tag = duration_tag_local(duration_hours)
tag = sprintf('%gh', duration_hours);
tag = strrep(tag, '.', 'p');
tag = regexprep(tag, '[^a-zA-Z0-9_]', '_');
end

function save_paper_figure_local(fig, png_path, pdf_path, fig_path)
drawnow;
exportgraphics(fig, png_path, 'Resolution', 300, ...
    'BackgroundColor', 'white');
exportgraphics(fig, pdf_path, 'ContentType', 'vector', ...
    'BackgroundColor', 'white');
savefig(fig, fig_path);
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
        error('example_skyplot_time_comparison:RepoRootNotFound', ...
            'Cannot locate SatLock repo root from %s.', start_dir);
    end
    repo_dir = parent_dir;
end
end
