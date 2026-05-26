function [manifest_tbl, out_dir, demo] = generate_system_design_figures(user_cfg)
% GENERATE_SYSTEM_DESIGN_FIGURES
% Build principle-oriented figures for Section III. SYSTEM DESIGN in
% StarDial. This workflow is intentionally isolated from the main
% gesture_test -> paper_figures pipeline.

if nargin < 1 || isempty(user_cfg)
    user_cfg = struct();
end

repo_dir = fileparts(fileparts(mfilename('fullpath')));
old_dir = pwd;
cleanup_obj = onCleanup(@() cd(old_dir)); %#ok<NASGU>
cd(repo_dir);

add_demo_paths_local(repo_dir);

cfg = default_cfg_local(repo_dir);
cfg = merge_cfg_local(cfg, user_cfg);

[out_dir, fmt_dirs] = prepare_output_dirs_local(cfg.output_root);
demo = build_demo_dataset_local(cfg);

manifest_rows = repmat(struct( ...
    'figure_id', "", ...
    'png_path', "", ...
    'pdf_path', "", ...
    'fig_path', ""), 0, 1);

fig_a = plot_signal_preprocessing_local(demo, cfg);
manifest_rows(end + 1, 1) = save_figure_bundle_local(fig_a, "A_signal_preprocessing", fmt_dirs, cfg); %#ok<AGROW>

fig_b1 = plot_diffraction_feature_extraction_local(demo, cfg);
manifest_rows(end + 1, 1) = save_figure_bundle_local(fig_b1, "B1_diffraction_waveforms", fmt_dirs, cfg); %#ok<AGROW>

fig_b2 = plot_diffraction_event_features_local(demo, cfg);
manifest_rows(end + 1, 1) = save_figure_bundle_local(fig_b2, "B2_diffraction_feature_values", fmt_dirs, cfg); %#ok<AGROW>

fig_b3 = plot_diffraction_saliency_scores_local(demo, cfg);
manifest_rows(end + 1, 1) = save_figure_bundle_local(fig_b3, "B3_diffraction_saliency_scores", fmt_dirs, cfg); %#ok<AGROW>

fig_c1 = plot_geometric_inversion_local(demo, cfg);
manifest_rows(end + 1, 1) = save_figure_bundle_local(fig_c1, "C1_gesture_plane_geometric_inversion", fmt_dirs, cfg); %#ok<AGROW>

fig_c2 = plot_geometric_constraint_rows_local(demo, cfg);
manifest_rows(end + 1, 1) = save_figure_bundle_local(fig_c2, "C2_geometric_constraint_rows", fmt_dirs, cfg); %#ok<AGROW>

manifest_tbl = struct2table(manifest_rows);
postprocess_pdf_groups_local(manifest_tbl);
writetable(manifest_tbl, fullfile(out_dir, 'system_design_manifest.csv'));
save(fullfile(out_dir, 'system_design_demo_summary.mat'), 'demo', 'cfg', 'manifest_tbl');

fprintf('\nSystem-design figures exported to:\n%s\n', out_dir);
disp(manifest_tbl);
end

function cfg = default_cfg_local(repo_dir)
cfg = struct();
cfg.obs_filepath = fullfile(repo_dir, 'data', '1_8', 'A_1_8_1.obs');
cfg.nav_filepath = fullfile(repo_dir, 'data', '1_8', '2026_1_8.nav');
cfg.output_root = fullfile(repo_dir, 'system_polt', 'results');
cfg.figure_resolution = 320;
cfg.sample_rate_hz = 25;
cfg.window_pad_sec = 1.2;
cfg.max_display_segments = 2;
cfg.max_top_sats = 3;
cfg.feature_window_pad_sec = 0.8;
cfg.feature_demo_display_pad_sec = 0.7;
cfg.feature_demo_max_segments = 4;
cfg.feature_demo_gap_fill = 4;
cfg.feature_demo_min_pts = 10;
cfg.feature_demo_continuity_pts = 11;
cfg.feature_demo_depth_floor = 0.8;
cfg.feature_demo_osc_floor = 0.25;
cfg.geometry_active_satellites = 5;
cfg.geometry_silent_satellites = 3;
cfg.geometry_sample_frames = 12;
cfg.geometry_min_elevation_deg = 5;
cfg.geometry_max_projection_radius = 0.90;
cfg.geometry_contact_label_frames = 3;
cfg.geometry_ref_epoch_count = 8;
cfg.geometry_interaction_half_span = 0.25;
cfg.geometry_gesture_height_m = 0.30;
cfg.geometry_effective_radius_m = 0.07;
cfg.geometry_silent_max_radius = 0.32;
cfg.geometry_path_link_samples = 20;

cfg.recovery_cfg = struct();
cfg.recovery_cfg.debug = struct('verbose', false, 'plot', false);
cfg.recovery_cfg.model = struct('max_hand_radius', 0.40);
cfg.recovery_cfg.grid = struct('x_min', -0.35, 'x_max', 0.35, 'y_min', -0.35, 'y_max', 0.35, 'step', 0.015);
cfg.recovery_cfg.track = struct( ...
    'lambda_smooth', 12.0, ...
    'final_smooth_pts', 2, ...
    'max_jump_m', 0.18, ...
    'use_active_interpolation', true, ...
    'use_process_window_output', true, ...
    'output_pad_frames', 14, ...
    'use_draw_mask_output', false, ...
    'use_draw_energy_gate', false, ...
    'enforce_piecewise_linear', true, ...
    'polyline_min_segments', 1, ...
    'polyline_rdp_eps', 0.022, ...
    'polyline_corner_angle_deg', 26, ...
    'polyline_max_fit_err', 0.19, ...
    'polyline_len_ratio_min', 0.60, ...
    'polyline_len_ratio_max', 1.25, ...
    'template_snap_enable', false, ...
    'endpoint_lock_enable', true, ...
    'endpoint_lock_blend', 0.72, ...
    'endpoint_lock_len_pts', 10, ...
    'axis_regularize_enable', true, ...
    'axis_regularize_min_major_span', 0.16, ...
    'axis_regularize_max_minor_span', 0.08, ...
    'axis_regularize_min_aspect', 3.2, ...
    'axis_regularize_monotonicity_min', 0.78, ...
    'axis_regularize_path_ratio_max', 1.85, ...
    'axis_regularize_max_turns', 2, ...
    'axis_regularize_target_span', 0.42, ...
    'axis_regularize_max_span', 0.50, ...
    'axis_regularize_blend', 0.82, ...
    'axis_regularize_minor_keep', 0.15, ...
    'shape_guided_enable', false);
end

function add_demo_paths_local(repo_dir)
paths = { ...
    fullfile(repo_dir, 'obs_parse'), ...
    fullfile(repo_dir, 'nav_parse'), ...
    fullfile(repo_dir, 'calculate_clock_bias_and_positon'), ...
    fullfile(repo_dir, 'gesture_analysis', 'preprocess_feature_extraction'), ...
    fullfile(repo_dir, 'gesture_analysis', 'continue'), ...
    fullfile(repo_dir, 'gesture_analysis', 'utils'), ...
    fullfile(repo_dir, 'gesture_analysis', 'templates')};

for i = 1:numel(paths)
    if exist(paths{i}, 'dir') == 7
        addpath(paths{i});
    end
end
end

function [out_dir, fmt_dirs] = prepare_output_dirs_local(output_root)
stamp = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
out_dir = fullfile(output_root, ['system_design_figures_', stamp]);
fmt_dirs = struct();
fmt_dirs.png = fullfile(out_dir, 'png');
fmt_dirs.pdf = fullfile(out_dir, 'pdf');
fmt_dirs.fig = fullfile(out_dir, 'fig');

if exist(out_dir, 'dir') ~= 7
    mkdir(out_dir);
end
if exist(fmt_dirs.png, 'dir') ~= 7
    mkdir(fmt_dirs.png);
end
if exist(fmt_dirs.pdf, 'dir') ~= 7
    mkdir(fmt_dirs.pdf);
end
if exist(fmt_dirs.fig, 'dir') ~= 7
    mkdir(fmt_dirs.fig);
end
end

function demo = build_demo_dataset_local(cfg)
demo = build_schematic_demo_dataset_local(cfg);
end

function demo = build_schematic_demo_dataset_local(cfg)
rng(26);

sample_idx = (1:round(12 * cfg.sample_rate_hz)).';
time_s = (sample_idx - 1) / cfg.sample_rate_hz;
display_windows = [2.0 3.4; 5.0 6.2; 8.1 9.5];
display_segments_rel = build_segments_from_windows_local(time_s, display_windows);

[sat1_raw, sat1_clean, sat1_base] = build_signal_channel_schematic_local( ...
    time_s, 45.0, display_windows, ...
    [2.6 1.7 2.3], [1.8 2.2 1.6], [2.2 2.8 2.0], 0.22, 0.75);
sat1_raw = build_sat_a_reference_raw_local(sample_idx, sat1_base);
[sat2_raw, sat2_clean, sat2_base] = build_signal_channel_schematic_local( ...
    time_s, 43.6, display_windows, ...
    [1.9 2.4 1.8], [2.4 1.9 2.5], [2.4 2.0 2.6], 0.20, 1.35);
[sat3_raw, sat3_clean, sat3_base] = build_signal_channel_schematic_local( ...
    time_s, 46.1, display_windows, ...
    [2.2 2.0 2.7], [1.7 2.6 1.9], [1.8 2.3 2.7], 0.18, 2.10);

feature_demo = build_schematic_feature_demo_local(cfg);
geometry_demo = build_schematic_geometry_demo_local(cfg);

combined_dev = zeros(size(time_s));
clean_stack = [sat1_clean, sat2_clean, sat3_clean];
baseline_stack = [sat1_base, sat2_base, sat3_base];
for i = 1:size(clean_stack, 2)
    baseline_i = baseline_stack(i);
    diff_i = max(0, baseline_i - clean_stack(:, i));
    combined_dev = combined_dev + (diff_i / max(max(diff_i), 1e-6));
end
combined_dev = movmean_local(combined_dev, 11);
combined_dev = combined_dev / max(max(combined_dev), 1e-6);
combined_dev = 0.06 + 1.18 * combined_dev;
gvi_threshold = 0.48;

demo = struct();
demo.sample_name = "schematic illustration";
demo.obs_filepath = "";
demo.nav_filepath = "";
demo.t_window = seconds(time_s);
demo.time_s = time_s;
demo.sample_idx = sample_idx;
demo.sample_rate_hz = cfg.sample_rate_hz;
demo.segment_idx = 2;
demo.display_segments_rel = display_segments_rel;
demo.gvi_window = combined_dev;
demo.gvi_threshold = gvi_threshold;
demo.window_mask = true(size(time_s));
demo.top_sat_ids = {'Sat-A', 'Sat-B', 'Sat-C'};
demo.top_sat_raw = [sat1_raw, sat2_raw, sat3_raw];
demo.top_sat_clean = [sat1_clean, sat2_clean, sat3_clean];
demo.top_sat_baselines = round([sat1_base, sat2_base, sat3_base]);
demo.top_sat_scores = [0.94, 0.88, 0.91];
demo.feature_demo = feature_demo;
demo.geometry_demo = geometry_demo;
end

function [raw_wave, clean_wave, baseline] = build_signal_channel_schematic_local( ...
    t_s, baseline, windows, amps, freqs, widths, noise_scale, phase_shift)
clean_wave = baseline + 0.20 * sin(2 * pi * 0.07 * t_s + 0.45 * phase_shift);
raw_wave = baseline ...
    + 0.05 * sin(2 * pi * 0.04 * t_s + 0.28 * phase_shift) ...
    + 0.04 * cos(2 * pi * 0.08 * t_s + 0.18 * phase_shift);
for k = 1:size(windows, 1)
    center = mean(windows(k, :));
    sigma = 0.14 * widths(k);
    env = exp(-0.5 * ((t_s - center) / max(sigma, 1e-3)).^2);
    osc = 0.30 + 0.70 * (0.5 + 0.5 * cos(2 * pi * freqs(k) * (t_s - center) + phase_shift));
    dip = amps(k) * env .* osc;
    clean_wave = clean_wave - dip;
    raw_style = 1 + mod(k + round(phase_shift * 3), 3);
    raw_profile = diffraction_profile_local((t_s - center) / max(widths(k), 1e-3), raw_style);
    raw_wave = raw_wave - amps(k) * (0.92 + 0.15 * noise_scale) * raw_profile;
end
clean_wave = movmean_local(clean_wave, 7);
raw_wave = quantize_signal_local(raw_wave, 1.0);
end

function raw_wave = build_sat_a_reference_raw_local(sample_idx, baseline)
raw_wave = baseline * ones(size(sample_idx));

raw_wave(21:30) = baseline - 1;
raw_wave(43:46) = baseline + 1;

raw_wave(51:54) = baseline - 1;
raw_wave(55:57) = baseline - 2;
raw_wave(58:60) = baseline - 3;
raw_wave(61:63) = baseline - 4;
raw_wave(64:66) = baseline - 2;
raw_wave(67:68) = baseline - 1;
raw_wave(69:70) = baseline - 3;
raw_wave(71:73) = baseline - 2;
raw_wave(74:76) = baseline - 1;
raw_wave(77:90) = baseline;
raw_wave(91:110) = baseline + 1;

raw_wave(111:118) = baseline;
raw_wave(119:121) = baseline - 1;
raw_wave(122:124) = baseline;
raw_wave(125:129) = baseline - 1;
raw_wave(130:133) = baseline - 2;
raw_wave(134:136) = baseline - 1;
raw_wave(137:141) = baseline - 3;
raw_wave(142:145) = baseline - 2;
raw_wave(146:148) = baseline - 1;
raw_wave(149:156) = baseline;
raw_wave(157:168) = baseline + 1;

raw_wave(169:186) = baseline;
raw_wave(187:202) = baseline - 1;
raw_wave(203:210) = baseline;
raw_wave(211:214) = baseline - 1;
raw_wave(215:217) = baseline - 2;
raw_wave(218:220) = baseline - 4;
raw_wave(221:223) = baseline - 1;
raw_wave(224:226) = baseline - 2;
raw_wave(227:229) = baseline - 1;
raw_wave(230:233) = baseline - 4;
raw_wave(234:237) = baseline - 3;
raw_wave(238:244) = baseline - 1;
raw_wave(245:254) = baseline;
raw_wave(255:266) = baseline + 1;
raw_wave(267:280) = baseline;
raw_wave(281:end) = baseline - 1;
end

function raw_wave = build_feature_raw_from_a_style_local( ...
    t_s, baseline, event_windows, raw_drop_levels, event_style_ids, w_mix, width_scales)
sample_idx = (1:numel(t_s)).';
sample_dt = median(diff(t_s));
if ~isfinite(sample_dt) || sample_dt <= 0
    sample_rate_hz = 25;
else
    sample_rate_hz = round(1 / sample_dt);
end

if nargin < 5 || isempty(event_style_ids)
    event_style_ids = 1:numel(raw_drop_levels);
end
if nargin < 6 || isempty(w_mix)
    w_mix = 0.55 * ones(size(raw_drop_levels));
end
if nargin < 7 || isempty(width_scales)
    width_scales = ones(size(raw_drop_levels));
end

event_segments = build_segments_from_windows_local(t_s, event_windows);
template_bank = build_a_style_event_templates_local(baseline, sample_rate_hz);
raw_wave = baseline + build_quiet_integer_noise_local(sample_idx, event_segments);

for i = 1:numel(event_segments)
    template_i = template_bank(1 + mod(i - 1, numel(template_bank)));
    seg_i = expand_segment_local(event_segments(i), numel(sample_idx), 6);
    n_i = seg_i.end_idx - seg_i.start_idx + 1;
    profile_i = build_blended_event_profile_local( ...
        template_i.raw_drop, n_i, event_style_ids(i), w_mix(i), width_scales(i), 0);
    drop_i = quantize_signal_local(raw_drop_levels(i) * profile_i, 1.0);
    drop_i = limit_adjacent_delta_local(drop_i, 2.0);
    raw_wave(seg_i.start_idx:seg_i.end_idx) = baseline - drop_i;
end

raw_wave = quantize_signal_local(raw_wave, 1.0);
raw_wave = min(baseline, max(35, raw_wave));
end

function clean_wave = build_feature_clean_from_a_style_local( ...
    t_s, baseline, event_windows, clean_drop_levels, event_style_ids, w_mix, width_scales)
sample_idx = (1:numel(t_s)).';
sample_dt = median(diff(t_s));
if ~isfinite(sample_dt) || sample_dt <= 0
    sample_rate_hz = 25;
else
    sample_rate_hz = round(1 / sample_dt);
end

if nargin < 5 || isempty(event_style_ids)
    event_style_ids = 1:numel(clean_drop_levels);
end
if nargin < 6 || isempty(w_mix)
    w_mix = 0.45 * ones(size(clean_drop_levels));
end
if nargin < 7 || isempty(width_scales)
    width_scales = ones(size(clean_drop_levels));
end

event_segments = build_segments_from_windows_local(t_s, event_windows);
template_bank = build_a_style_event_templates_local(baseline, sample_rate_hz);
clean_wave = baseline ...
    + 0.08 * sin(2 * pi * 0.040 * t_s + 0.30) ...
    + 0.05 * cos(2 * pi * 0.018 * t_s + 0.72);

for i = 1:numel(event_segments)
    template_i = template_bank(1 + mod(i - 1, numel(template_bank)));
    seg_i = expand_segment_local(event_segments(i), numel(sample_idx), 6);
    n_i = seg_i.end_idx - seg_i.start_idx + 1;
    profile_i = build_blended_event_profile_local( ...
        template_i.clean_drop, n_i, event_style_ids(i), w_mix(i), width_scales(i), 5);
    clean_wave(seg_i.start_idx:seg_i.end_idx) = clean_wave(seg_i.start_idx:seg_i.end_idx) ...
        - clean_drop_levels(i) * profile_i;
end

clean_wave = movmean_local(clean_wave, 3);
end

function profile_out = build_blended_event_profile_local(template_profile, n_out, style_id, w_mix, width_scale, smooth_win)
profile_a = resample_profile_local(template_profile, n_out);
profile_a = profile_a / max(max(profile_a), 1e-6);

u = linspace(-1, 1, n_out).' / max(width_scale, 1e-3);
profile_w = diffraction_w_profile_local(u, style_id);
if nargin >= 6 && smooth_win > 1
    profile_w = movmean_local(profile_w, smooth_win);
end
profile_w = profile_w / max(max(profile_w), 1e-6);

w_mix = min(1, max(0, w_mix));
profile_out = (1 - w_mix) * profile_a + w_mix * profile_w;
profile_out = profile_out / max(max(profile_out), 1e-6);
end

function profile = diffraction_w_profile_local(u, style_id)
u_nodes = [-0.82, -0.66, -0.50, -0.34, -0.22, -0.12, -0.02, 0.10, 0.22, 0.34, 0.48, 0.64, 0.82];
switch mod(style_id - 1, 3) + 1
    case 1
        v_nodes = [0.00, 0.04, 0.16, 0.46, 0.90, 1.00, 0.68, 0.60, 0.88, 0.76, 0.30, 0.08, 0.00];
    case 2
        v_nodes = [0.00, 0.03, 0.12, 0.40, 0.82, 0.98, 0.66, 0.58, 0.92, 0.82, 0.36, 0.10, 0.00];
    otherwise
        v_nodes = [0.00, 0.05, 0.18, 0.50, 0.94, 1.00, 0.70, 0.62, 0.90, 0.80, 0.34, 0.09, 0.00];
end

profile = interp1(u_nodes, v_nodes, u, 'pchip', 0);
profile = max(profile, 0);
end

function segs = build_segments_from_windows_local(t_s, windows)
segs = repmat(struct('start_idx', 1, 'end_idx', 1), size(windows, 1), 1);
for i = 1:size(windows, 1)
    [~, start_idx] = min(abs(t_s - windows(i, 1)));
    [~, end_idx] = min(abs(t_s - windows(i, 2)));
    segs(i).start_idx = min(start_idx, end_idx);
    segs(i).end_idx = max(start_idx, end_idx);
end
end

function feature_demo = build_schematic_feature_demo_local(cfg)
rng(31);
sample_idx = (1:round(20.0 * cfg.sample_rate_hz)).';
t_s = (sample_idx - 1) / cfg.sample_rate_hz;
baseline = 45.0;

event_windows = [ ...
    2.00 3.40; ...
    5.00 6.20; ...
    8.10 9.50; ...
   11.20 12.40; ...
   14.10 15.50; ...
   17.20 18.60];
event_labels = {'e_1', 'e_2', 'e_3', 'e_4', 'e_5', 'e_6'};
clean_drop_levels = [0.92, 1.72, 4.45, 1.28, 3.18, 2.46];
raw_drop_levels = [2, 5, 12, 3, 8, 6];
event_style_ids = [1, 2, 3, 1, 2, 3];
raw_w_mix = [0.42, 0.66, 1.00, 0.52, 0.84, 0.76];
clean_w_mix = [0.30, 0.54, 0.98, 0.38, 0.74, 0.64];
event_width_scales = [0.86, 0.96, 1.24, 0.90, 1.06, 1.00];

event_segments = build_segments_from_windows_local(t_s, event_windows);
num_events = numel(event_segments);
event_center_s = zeros(1, num_events);
clean_wave = build_feature_clean_from_a_style_local( ...
    t_s, baseline, event_windows, clean_drop_levels, ...
    event_style_ids, clean_w_mix, event_width_scales);
raw_wave = build_feature_raw_from_a_style_local( ...
    t_s, baseline, event_windows, raw_drop_levels, ...
    event_style_ids, raw_w_mix, event_width_scales);

for i = 1:num_events
    event_center_s(i) = mean(event_windows(i, :));
end

event_center_idx = round(0.5 * ([event_segments.start_idx] + [event_segments.end_idx]));
[event_depth, event_osc, event_cont] = measure_diffraction_features_local( ...
    raw_wave, clean_wave, baseline, event_segments);

amp_norm = normalize_vector_local(event_depth);
osc_norm = normalize_vector_local(event_osc);
cont_norm = normalize_vector_local(event_cont);
saliency_norm = 0.50 * amp_norm + 0.30 * osc_norm + 0.20 * cont_norm;
saliency_norm = normalize_vector_local(saliency_norm);

feature_demo = struct();
feature_demo.t_s = t_s;
feature_demo.sample_idx = sample_idx;
feature_demo.raw_wave = raw_wave;
feature_demo.clean_wave = clean_wave;
feature_demo.baseline = baseline;
feature_demo.event_segments = event_segments;
feature_demo.event_labels = event_labels;
feature_demo.event_center_s = event_center_s;
feature_demo.event_center_idx = event_center_idx;
feature_demo.amp_norm = amp_norm;
feature_demo.osc_norm = osc_norm;
feature_demo.cont_norm = cont_norm;
feature_demo.saliency_norm = saliency_norm;
feature_demo.feature_mat = [amp_norm(:), osc_norm(:), cont_norm(:)];
end

function y = normalize_vector_local(x)
x = double(x(:));
x_min = min(x);
x_max = max(x);
if ~isfinite(x_min) || ~isfinite(x_max) || (x_max - x_min) < 1e-6
    y = ones(size(x));
else
    y = (x - x_min) / (x_max - x_min);
end
y = 0.18 + 0.82 * y;
y = y(:).';
end

function y_q = quantize_signal_local(y, step)
if nargin < 2 || ~isfinite(step) || step <= 0
    y_q = y;
    return;
end
y_q = step * round(y / step);
end

function template_bank = build_a_style_event_templates_local(baseline, sample_rate_hz)
sample_idx_ref = (1:round(12 * sample_rate_hz)).';
t_ref = (sample_idx_ref - 1) / sample_rate_hz;
windows_ref = [2.0 3.4; 5.0 6.2; 8.1 9.5];
[~, clean_ref, ~] = build_signal_channel_schematic_local( ...
    t_ref, baseline, windows_ref, ...
    [2.6 1.7 2.3], [1.8 2.2 1.6], [2.2 2.8 2.0], 0.22, 0.75);
raw_ref = build_sat_a_reference_raw_local(sample_idx_ref, baseline);
segments_ref = build_segments_from_windows_local(t_ref, windows_ref);

template_bank = repmat(struct('raw_drop', [], 'clean_drop', []), numel(segments_ref), 1);
for i = 1:numel(segments_ref)
    seg_i = expand_segment_local(segments_ref(i), numel(sample_idx_ref), 6);
    raw_drop_i = max(0, baseline - raw_ref(seg_i.start_idx:seg_i.end_idx));
    clean_drop_i = max(0, baseline - clean_ref(seg_i.start_idx:seg_i.end_idx));
    template_bank(i).raw_drop = raw_drop_i / max(max(raw_drop_i), 1e-6);
    template_bank(i).clean_drop = clean_drop_i / max(max(clean_drop_i), 1e-6);
end
end

function seg_out = expand_segment_local(seg_in, n_total, pad_pts)
seg_out = struct( ...
    'start_idx', max(1, seg_in.start_idx - pad_pts), ...
    'end_idx', min(n_total, seg_in.end_idx + pad_pts));
end

function profile_out = resample_profile_local(profile_in, n_out)
profile_in = profile_in(:);
if n_out <= 1 || numel(profile_in) <= 1
    profile_out = profile_in;
    return;
end

u_in = linspace(0, 1, numel(profile_in)).';
u_out = linspace(0, 1, n_out).';
profile_out = interp1(u_in, profile_in, u_out, 'pchip');
profile_out = max(profile_out, 0);
end

function noise = build_quiet_integer_noise_local(sample_idx, event_segments)
sample_idx = sample_idx(:);
ctrl_idx = unique([1; (1:7:numel(sample_idx)).'; numel(sample_idx)]);
ctrl_wave = 0.55 + 0.78 * sin(0.17 * ctrl_idx) + 0.52 * cos(0.05 * ctrl_idx + 0.80);
ctrl_vals = -min(1, max(0, round(ctrl_wave)));
noise = interp1(ctrl_idx, ctrl_vals, (1:numel(sample_idx)).', 'previous', 'extrap');
noise = round(noise);

event_mask = false(size(sample_idx));
for i = 1:numel(event_segments)
    seg_i = expand_segment_local(event_segments(i), numel(sample_idx), 4);
    event_mask(seg_i.start_idx:seg_i.end_idx) = true;
end
noise(event_mask) = 0;
end

function [depth_vals, osc_vals, cont_vals] = measure_diffraction_features_local(raw_wave, clean_wave, baseline, event_segments)
n_events = numel(event_segments);
depth_vals = zeros(1, n_events);
osc_vals = zeros(1, n_events);
cont_vals = zeros(1, n_events);

for i = 1:n_events
    seg_i = expand_segment_local(event_segments(i), numel(raw_wave), 6);
    idx_i = seg_i.start_idx:seg_i.end_idx;
    raw_drop_i = max(0, baseline - raw_wave(idx_i));
    clean_drop_i = max(0, baseline - clean_wave(idx_i));
    peak_i = max(raw_drop_i);
    depth_vals(i) = max(clean_drop_i);
    osc_vals(i) = 0.72 * mean(raw_drop_i) + 0.28 * peak_i;
    if peak_i <= 0
        cont_vals(i) = 0;
    else
        width_ratio_i = nnz(raw_drop_i >= max(1.0, 0.35 * peak_i)) / numel(raw_drop_i);
        cont_vals(i) = peak_i * width_ratio_i;
    end
end
end

function [x_dense, y_dense] = densify_integer_series_local(x, y, factor)
x = x(:);
y = y(:);
factor = max(1, round(factor));
if numel(x) <= 1 || factor == 1
    x_dense = x;
    y_dense = y;
    return;
end

x_dense = linspace(x(1), x(end), factor * (numel(x) - 1) + 1).';
y_dense = interp1(x, y, x_dense, 'pchip');
y_dense = quantize_signal_local(y_dense, 1.0);
y_dense = min(max(y_dense, min(y)), max(y));
end

function [x_dense, y_dense] = densify_curve_local(x, y, factor)
x = x(:);
y = y(:);
factor = max(1, round(factor));
if numel(x) <= 1 || factor == 1
    x_dense = x;
    y_dense = y;
    return;
end

x_dense = linspace(x(1), x(end), factor * (numel(x) - 1) + 1).';
y_dense = interp1(x, y, x_dense, 'pchip');
end

function y_out = limit_adjacent_delta_local(y_in, max_delta)
y_out = y_in(:);
if isempty(y_out)
    return;
end

max_delta = max(0, double(max_delta));
for pass = 1:2
    for i = 2:numel(y_out)
        delta_i = y_out(i) - y_out(i - 1);
        if delta_i > max_delta
            y_out(i) = y_out(i - 1) + max_delta;
        elseif delta_i < -max_delta
            y_out(i) = y_out(i - 1) - max_delta;
        end
    end
    for i = (numel(y_out) - 1):-1:1
        delta_i = y_out(i) - y_out(i + 1);
        if delta_i > max_delta
            y_out(i) = y_out(i + 1) + max_delta;
        elseif delta_i < -max_delta
            y_out(i) = y_out(i + 1) - max_delta;
        end
    end
end
end

function apply_axis_typography_local(ax, tick_font_size, label_font_size)
if nargin < 2 || isempty(tick_font_size)
    tick_font_size = 14;
end
if nargin < 3 || isempty(label_font_size)
    label_font_size = 16;
end

set(ax, ...
    'FontName', plot_font_name_local(), ...
    'FontSize', tick_font_size, ...
    'Box', 'on', ...
    'LineWidth', 1.0, ...
    'TickDir', 'in', ...
    'Layer', 'top');
if isprop(ax, 'Title') && isgraphics(ax.Title)
    set(ax.Title, 'FontName', plot_font_name_local(), 'FontSize', label_font_size);
end
if isprop(ax, 'XLabel') && isgraphics(ax.XLabel)
    set(ax.XLabel, 'FontName', plot_font_name_local(), 'FontSize', label_font_size);
end
if isprop(ax, 'YLabel') && isgraphics(ax.YLabel)
    set(ax.YLabel, 'FontName', plot_font_name_local(), 'FontSize', label_font_size);
end
end

function font_name = plot_font_name_local()
font_name = 'Arial';
end

function profile = diffraction_profile_local(u, style_id)
switch mod(style_id - 1, 3) + 1
    case 1
        u_nodes = [-0.85, -0.64, -0.46, -0.31, -0.20, -0.10, 0.00, 0.10, 0.22, 0.38, 0.58, 0.82];
        v_nodes = [ 0.00,  0.02,  0.08,  0.18,  0.60,  0.24, 1.00, 0.30, 0.58, 0.16, 0.05, 0.00];
    case 2
        u_nodes = [-0.85, -0.66, -0.50, -0.36, -0.22, -0.08, 0.06, 0.18, 0.30, 0.46, 0.64, 0.82];
        v_nodes = [ 0.00,  0.03,  0.10,  0.46,  0.18,  0.78, 0.22, 0.66, 0.26, 0.12, 0.04, 0.00];
    otherwise
        u_nodes = [-0.85, -0.68, -0.52, -0.38, -0.24, -0.10, 0.02, 0.16, 0.30, 0.44, 0.60, 0.82];
        v_nodes = [ 0.00,  0.02,  0.06,  0.20,  0.52,  0.16, 0.88, 0.38, 0.18, 0.74, 0.14, 0.00];
end

profile = interp1(u_nodes, v_nodes, u, 'pchip', 0);
profile = max(profile, 0);
end

function set_sample_axis_ticks_local(ax, x, sample_rate_hz)
x = x(:);
if isempty(x)
    return;
end

x_start = min(x);
x_end = max(x);
span = x_end - x_start;

if span <= 140
    tick_step = sample_rate_hz;
elseif span <= 320
    tick_step = 2 * sample_rate_hz;
else
    tick_step = 4 * sample_rate_hz;
end

xt = unique([x_start; (ceil(x_start / tick_step) * tick_step:tick_step:x_end).']);
if isempty(xt) || (x_end - xt(end)) > 0.35 * tick_step
    xt = unique([xt; x_end]);
end
set(ax, 'XTick', xt);
end

function labels = build_sample_event_ticklabels_local(x, event_labels)
x = x(:);
labels = cell(numel(x), 1);
for i = 1:numel(x)
    labels{i} = sprintf('%d\\newline%s', x(i), event_labels{i});
end
end

function geometry_demo = build_schematic_geometry_demo_local(cfg)
contact_points = [ ...
    -0.18  0.22; ...
     0.16  0.22; ...
    -0.03  0.02; ...
     0.20 -0.18; ...
    -0.20 -0.18];
num_active = size(contact_points, 1);
num_trace_frames = 14;
active_traces = nan(num_trace_frames, num_active, 2);
trace_dirs = [ ...
     0.10  0.00; ...
     0.10  0.00; ...
    -0.11 -0.10; ...
     0.11  0.00; ...
     0.00  0.00];

for i = 1:(num_active - 1)
    alpha = linspace(-1, 1, num_trace_frames).';
    local_curve = contact_points(i, :) + alpha .* trace_dirs(i, :);
    active_traces(:, i, :) = local_curve;
end

corner_ctrl = [ ...
    -0.24 -0.13; ...
    -0.20 -0.18; ...
    -0.12 -0.18];
corner_trace = interpolate_polyline_linear_local(corner_ctrl, num_trace_frames);
active_traces(:, num_active, :) = corner_trace;

traj_ctrl = [ ...
    -0.26  0.23; ...
     0.24  0.23; ...
    -0.20 -0.18; ...
     0.24 -0.18];
traj_pts = interpolate_polyline_linear_local(traj_ctrl, 180);

geometry_demo = struct();
geometry_demo.active_sat_ids = {'Sat-1', 'Sat-2', 'Sat-3', 'Sat-4', 'Sat-5'};
geometry_demo.silent_sat_ids = {'Silent-A', 'Silent-B', 'Silent-C'};
geometry_demo.active_traces = active_traces;
geometry_demo.active_mid_points = contact_points;
geometry_demo.silent_points = [ ...
    -0.05 -0.29; ...
     0.14  0.31; ...
     0.27  0.14];
geometry_demo.traj_x = traj_pts(:, 1);
geometry_demo.traj_y = traj_pts(:, 2);
geometry_demo.interaction_half_span = 0.34;
geometry_demo.effective_radius_m = 0.075;
geometry_demo.window_sample_idx = (1:round(6.0 * cfg.sample_rate_hz)).';
geometry_demo.window_time_s = (geometry_demo.window_sample_idx - 1) / cfg.sample_rate_hz;
geometry_demo.event_segments = { ...
    struct('t_start', 0.7, 't_end', 1.6), ...
    struct('t_start', 1.8, 't_end', 2.7), ...
    struct('t_start', 3.0, 't_end', 3.9), ...
    struct('t_start', 4.4, 't_end', 5.1), ...
    struct('t_start', 3.9, 't_end', 4.5)};
geometry_demo.active_match_time_s = [1.1; 2.2; 3.4; 4.8; 4.2];
geometry_demo.silent_match_time_s = [1.9; 3.8; 5.4];
geometry_demo.event_segments_sample = { ...
    struct('start_idx', round(0.7 * cfg.sample_rate_hz) + 1, 'end_idx', round(1.6 * cfg.sample_rate_hz) + 1), ...
    struct('start_idx', round(1.8 * cfg.sample_rate_hz) + 1, 'end_idx', round(2.7 * cfg.sample_rate_hz) + 1), ...
    struct('start_idx', round(3.0 * cfg.sample_rate_hz) + 1, 'end_idx', round(3.9 * cfg.sample_rate_hz) + 1), ...
    struct('start_idx', round(4.4 * cfg.sample_rate_hz) + 1, 'end_idx', round(5.1 * cfg.sample_rate_hz) + 1), ...
    struct('start_idx', round(3.9 * cfg.sample_rate_hz) + 1, 'end_idx', round(4.5 * cfg.sample_rate_hz) + 1)};
geometry_demo.active_match_sample_idx = round(geometry_demo.active_match_time_s * cfg.sample_rate_hz) + 1;
geometry_demo.silent_match_sample_idx = round(geometry_demo.silent_match_time_s * cfg.sample_rate_hz) + 1;
end

function pts = interpolate_polyline_local(ctrl_pts, n_total)
if size(ctrl_pts, 1) < 2
    pts = ctrl_pts;
    return;
end

seg_len = vecnorm(diff(ctrl_pts, 1, 1), 2, 2);
cum_len = [0; cumsum(seg_len)];
sample_pos = linspace(0, cum_len(end), n_total).';
pts = nan(n_total, 2);
for i = 1:2
    pts(:, i) = interp1(cum_len, ctrl_pts(:, i), sample_pos, 'pchip');
end
end

function pts = interpolate_polyline_linear_local(ctrl_pts, n_total)
if size(ctrl_pts, 1) < 2
    pts = ctrl_pts;
    return;
end

seg_len = vecnorm(diff(ctrl_pts, 1, 1), 2, 2);
cum_len = [0; cumsum(seg_len)];
sample_pos = linspace(0, cum_len(end), n_total).';
pts = nan(n_total, 2);
for i = 1:2
    pts(:, i) = interp1(cum_len, ctrl_pts(:, i), sample_pos, 'linear');
end
end

function raw_cn0_matrix = interpolate_raw_cn0_matrix_local(obs_data, valid_sats, t_grid)
num_samples = numel(t_grid);
num_sats = numel(valid_sats);
raw_cn0_matrix = nan(num_samples, num_sats);

for s = 1:num_sats
    sid = valid_sats{s};
    target_code = '';
    for k = 1:min(60, numel(obs_data))
        if ~isfield(obs_data(k), 'data') || isempty(obs_data(k).data) || ~isfield(obs_data(k).data, sid)
            continue;
        end
        sat_item = obs_data(k).data.(sid);
        if isfield(sat_item, 'snr') && ~isempty(sat_item.snr)
            fds = fieldnames(sat_item.snr);
            if ~isempty(fds)
                target_code = fds{1};
                break;
            end
        end
    end

    if isempty(target_code)
        continue;
    end

    s_times = [];
    s_vals = [];
    for k = 1:numel(obs_data)
        if ~isfield(obs_data(k), 'data') || isempty(obs_data(k).data) || ~isfield(obs_data(k).data, sid)
            continue;
        end
        sat_item = obs_data(k).data.(sid);
        if isfield(sat_item, 'snr') && isfield(sat_item.snr, target_code)
            val = sat_item.snr.(target_code);
            if isfinite(val) && val > 10
                s_times(end + 1, 1) = datenum(obs_data(k).time); %#ok<AGROW>
                s_vals(end + 1, 1) = val; %#ok<AGROW>
            end
        end
    end

    if numel(s_times) < 5
        continue;
    end

    [u_times, iu] = unique(s_times, 'stable');
    tq = datenum(t_grid);
    raw_cn0_matrix(:, s) = interp1(u_times, s_vals(iu), tq, 'linear', NaN);
end
end

function segment_idx = select_display_segment_local(step1)
segments = step1.segments;
if isempty(segments)
    segment_idx = 1;
    return;
end

if isfield(segments, 'peak_gvi')
    peak_vals = [segments.peak_gvi];
else
    peak_vals = [segments.end_idx] - [segments.start_idx] + 1;
end

[~, segment_idx] = max(peak_vals);
end

function segment_idx = select_display_segments_local(step1, max_segments)
segments = step1.segments;
if isempty(segments)
    segment_idx = 1;
    return;
end

if isfield(segments, 'peak_gvi')
    peak_vals = [segments.peak_gvi];
else
    peak_vals = [segments.end_idx] - [segments.start_idx] + 1;
end

[~, ord] = sort(peak_vals, 'descend');
segment_idx = sort(ord(1:min(max_segments, numel(ord))));
end

function [score, sat_order] = rank_display_satellites_local(clean_mat, vol_mat)
dev_score = sum(max(vol_mat, 0), 1, 'omitnan');
osc_score = nansum(abs(diff(clean_mat, 1, 1)), 1);

dev_n = normalize_nonnegative_local(dev_score(:));
osc_n = normalize_nonnegative_local(osc_score(:));
score = 0.55 * osc_n + 0.45 * dev_n;
score = score(:).';

[score_sorted, sat_order] = sort(score, 'descend');
sat_order = sat_order(score_sorted > 0);
end

function rel_seg = localize_segment_local(window_idx, seg)
rel_seg = struct('start_idx', max(1, seg.start_idx - window_idx(1) + 1), ...
    'end_idx', min(numel(window_idx), seg.end_idx - window_idx(1) + 1));
end

function rel_segments = localize_segments_local(window_idx, segments)
rel_segments = repmat(struct('start_idx', NaN, 'end_idx', NaN), numel(segments), 1);
for i = 1:numel(segments)
    rel_segments(i) = localize_segment_local(window_idx, segments(i));
end
end

function window_idx = build_display_window_local(t_grid, seg, window_pad_sec)
dt_sec = median(seconds(diff(t_grid)), 'omitnan');
if ~isfinite(dt_sec) || dt_sec <= 0
    dt_sec = 1 / 25;
end
pad_pts = max(1, round(window_pad_sec / dt_sec));
start_idx = min([seg.start_idx]);
end_idx = max([seg.end_idx]);
window_idx = max(1, start_idx - pad_pts):min(numel(t_grid), end_idx + pad_pts);
window_idx = window_idx(:);
end

function feature_demo = build_feature_demo_local(t_window, raw_wave, clean_wave, segment_rel, cfg)
feature_demo = struct();

t_s = seconds(t_window - t_window(1));
baseline = compute_single_baseline_local(clean_wave);
deviation = max(0, baseline - clean_wave);
oscillation = movmean_local(abs([0; diff(clean_wave)]), 5);
continuity = movmean_local(double(deviation > cfg.feature_demo_depth_floor), cfg.feature_demo_continuity_pts);

coarse_mask = false(numel(t_s), 1);
coarse_mask(max(1, segment_rel.start_idx):min(numel(t_s), segment_rel.end_idx)) = true;
deviation_norm = normalize_by_mask_local(deviation(:), coarse_mask);
oscillation_norm = normalize_by_mask_local(oscillation(:), coarse_mask);
continuity_norm = normalize_by_mask_local(continuity(:), coarse_mask);

segment_rel = refine_feature_segment_local( ...
    deviation_norm, oscillation_norm, continuity_norm, coarse_mask, segment_rel, cfg);

segment_mask = false(numel(t_s), 1);
segment_mask(segment_rel.start_idx:segment_rel.end_idx) = true;
display_amp = deviation_norm;
display_amp(~segment_mask) = 0;
display_osc = oscillation_norm;
display_osc(~segment_mask) = 0;
display_cont = continuity_norm;
display_cont(~segment_mask) = 0;

crop_idx = build_feature_crop_idx_local(t_s, segment_rel, cfg.feature_demo_display_pad_sec);
t_s = t_s(crop_idx);
raw_wave = raw_wave(crop_idx);
clean_wave = clean_wave(crop_idx);
deviation = deviation(crop_idx);
oscillation = oscillation(crop_idx);
continuity = continuity(crop_idx);
deviation_norm = deviation_norm(crop_idx);
oscillation_norm = oscillation_norm(crop_idx);
continuity_norm = continuity_norm(crop_idx);
display_amp = display_amp(crop_idx);
display_osc = display_osc(crop_idx);
display_cont = display_cont(crop_idx);
segment_rel = localize_segment_local(crop_idx, segment_rel);

seg_slice = segment_rel.start_idx:segment_rel.end_idx;

segment_amp_norm = clamp01_local(mean(deviation(seg_slice), 'omitnan') / max(eps, max(deviation, [], 'omitnan')));
segment_osc_norm = clamp01_local(max(oscillation(seg_slice), [], 'omitnan') / max(eps, max(oscillation, [], 'omitnan')));
segment_cont_norm = clamp01_local(mean(continuity(seg_slice), 'omitnan'));
segment_score_norm = clamp01_local(0.45 * segment_amp_norm + 0.35 * segment_osc_norm + 0.20 * segment_cont_norm);

seg_rows = struct( ...
    'segment_id', "e_{i,m}", ...
    'start_idx', segment_rel.start_idx, ...
    'end_idx', segment_rel.end_idx, ...
    'amplitude', segment_amp_norm, ...
    'oscillation', segment_osc_norm, ...
    'continuity', segment_cont_norm, ...
    'score', segment_score_norm, ...
    'level', "High");

feature_demo.segment_table = struct2table(seg_rows);
feature_demo.segment_amp_norm = segment_amp_norm;
feature_demo.segment_osc_norm = segment_osc_norm;
feature_demo.segment_cont_norm = segment_cont_norm;
feature_demo.segment_score_norm = segment_score_norm;

feature_demo.t_s = t_s;
feature_demo.raw_wave = raw_wave(:);
feature_demo.clean_wave = clean_wave(:);
feature_demo.baseline = baseline;
feature_demo.deviation = deviation(:);
feature_demo.oscillation = oscillation(:);
feature_demo.continuity = continuity(:);
feature_demo.deviation_norm = display_amp(:);
feature_demo.oscillation_norm = display_osc(:);
feature_demo.continuity_norm = display_cont(:);
feature_demo.segment_rel = segment_rel;
end

function geometry_demo = build_geometry_demo_local(obs_data, nav_data, step1, step1_shaped, window_idx, sat_order, traj_x, traj_y, traj_t, traj_conf, cfg)
geometry_demo = struct();

window_energy = sum(step1_shaped.volatility_matrix(window_idx, :), 1, 'omitnan');
active_sat_idx = sat_order(window_energy(sat_order) > 0);
active_sat_idx = active_sat_idx(1:min(cfg.geometry_active_satellites, numel(active_sat_idx)));
if isempty(active_sat_idx)
    active_sat_idx = sat_order(1:min(cfg.geometry_active_satellites, numel(sat_order)));
end

obs_time_num = datenum([obs_data.time]);
window_time_num = datenum(step1.t_grid(window_idx));
epoch_idx_window = round(interp1(obs_time_num, 1:numel(obs_time_num), window_time_num, 'nearest', 'extrap'));
epoch_idx_window = min(max(epoch_idx_window, 1), numel(obs_data));

event_frame_mask = max(step1_shaped.volatility_matrix(window_idx, active_sat_idx), [], 2) > 0;
candidate_frames = window_idx(event_frame_mask);
if numel(candidate_frames) < 4
    candidate_frames = window_idx;
end

sample_count = min(cfg.geometry_sample_frames, numel(candidate_frames));
sample_sel = round(linspace(1, numel(candidate_frames), sample_count));
sample_t_idx = candidate_frames(sample_sel);
sample_t_idx = unique(sample_t_idx, 'stable');
sample_epoch_idx = epoch_idx_window(ismember(window_idx, sample_t_idx));
sample_epoch_idx = unique(sample_epoch_idx, 'stable');

ref_epochs = sample_epoch_idx;
if numel(ref_epochs) > cfg.geometry_ref_epoch_count
    ref_sel = round(linspace(1, numel(ref_epochs), cfg.geometry_ref_epoch_count));
    ref_epochs = ref_epochs(ref_sel);
end
geo_ref = estimate_geometry_reference_local(obs_data, nav_data, ref_epochs);

mid_epoch = sample_epoch_idx(max(1, ceil(numel(sample_epoch_idx) / 2)));
mid_proj_all = compute_epoch_projection_local(obs_data, nav_data, mid_epoch, step1.valid_sats, geo_ref, cfg);

silent_candidates = find(all(isfinite(mid_proj_all), 2));
silent_candidates = setdiff(silent_candidates, active_sat_idx, 'stable');
if ~isempty(silent_candidates)
    [~, silent_ord] = sort(window_energy(silent_candidates), 'ascend');
    silent_candidates = silent_candidates(silent_ord);
end
silent_sat_idx = silent_candidates(1:min(cfg.geometry_silent_satellites, numel(silent_candidates)));

active_traces = nan(numel(sample_epoch_idx), numel(active_sat_idx), 2);
for i = 1:numel(sample_epoch_idx)
    proj_i = compute_epoch_projection_local(obs_data, nav_data, sample_epoch_idx(i), step1.valid_sats, geo_ref, cfg);
    if isempty(proj_i)
        continue;
    end
    active_traces(i, :, :) = proj_i(active_sat_idx, :);
end

traj_mask = (traj_t >= window_idx(1)) & (traj_t <= window_idx(end));
if nnz(traj_mask) < 8
    traj_mask = true(size(traj_t));
end

traj_window_idx = traj_t(traj_mask);
traj_window_x = traj_x(traj_mask);
traj_window_y = traj_y(traj_mask);
traj_window_conf = traj_conf(traj_mask);

active_proj_on_traj = nan(numel(traj_window_idx), numel(active_sat_idx), 2);
silent_proj_on_traj = nan(numel(traj_window_idx), numel(silent_sat_idx), 2);
for i = 1:numel(traj_window_idx)
    epoch_i = epoch_idx_window(find(window_idx == traj_window_idx(i), 1, 'first'));
    if isempty(epoch_i)
        [~, nearest_pos] = min(abs(window_idx - traj_window_idx(i)));
        epoch_i = epoch_idx_window(nearest_pos);
    end
    proj_i = compute_epoch_projection_local(obs_data, nav_data, epoch_i, step1.valid_sats, geo_ref, cfg);
    if isempty(proj_i)
        continue;
    end
    if ~isempty(active_sat_idx)
        active_proj_on_traj(i, :, :) = proj_i(active_sat_idx, :);
    end
    if ~isempty(silent_sat_idx)
        silent_proj_on_traj(i, :, :) = proj_i(silent_sat_idx, :);
    end
end

if isempty(silent_sat_idx)
    silent_points = zeros(0, 2);
else
    silent_points = mid_proj_all(silent_sat_idx, :);
end

event_segments = cell(numel(active_sat_idx), 1);
active_support_centers_s = nan(numel(active_sat_idx), 1);
for i = 1:numel(active_sat_idx)
    local_mask = step1_shaped.volatility_matrix(window_idx, active_sat_idx(i)) > 0;
    segs = mask_to_segments_local(local_mask);
    rows = repmat(struct('t_start', NaN, 't_end', NaN), numel(segs), 1);
    for j = 1:numel(segs)
        rows(j).t_start = seconds(step1.t_grid(window_idx(segs(j).start_idx)) - step1.t_grid(window_idx(1)));
        rows(j).t_end = seconds(step1.t_grid(window_idx(segs(j).end_idx)) - step1.t_grid(window_idx(1)));
    end
    event_segments{i} = rows;
    if ~isempty(rows)
        dur = arrayfun(@(r) r.t_end - r.t_start, rows);
        [~, best_idx] = max(dur);
        active_support_centers_s(i) = 0.5 * (rows(best_idx).t_start + rows(best_idx).t_end);
    end
end

sample_time_s = seconds(step1.t_grid(sample_t_idx) - step1.t_grid(window_idx(1)));
active_contact_points = nan(numel(active_sat_idx), 2);
for i = 1:numel(active_sat_idx)
    active_contact_points(i, :) = representative_contact_point_local( ...
        squeeze(active_traces(:, i, :)), sample_time_s, active_support_centers_s(i));
end

active_keep = all(isfinite(active_contact_points), 2);
active_keep = active_keep & (hypot(active_contact_points(:, 1), active_contact_points(:, 2)) <= ...
    (cfg.geometry_interaction_half_span + 0.18));
if nnz(active_keep) < min(2, numel(active_keep))
    active_keep = all(isfinite(active_contact_points), 2);
end

active_sat_idx = active_sat_idx(active_keep);
active_traces = active_traces(:, active_keep, :);
event_segments = event_segments(active_keep);
active_support_centers_s = active_support_centers_s(active_keep);
active_contact_points = active_contact_points(active_keep, :);
active_proj_on_traj = active_proj_on_traj(:, active_keep, :);

if isempty(silent_sat_idx)
    silent_points = zeros(0, 2);
else
    silent_points = mid_proj_all(silent_sat_idx, :);
end

if isempty(active_contact_points)
    active_centroid = [0, 0];
else
    active_centroid = mean(active_contact_points, 1, 'omitnan');
end

if ~isempty(silent_sat_idx)
    silent_keep = all(isfinite(silent_points), 2);
    silent_keep = silent_keep & (sqrt(sum((silent_points - active_centroid) .^ 2, 2)) <= cfg.geometry_silent_max_radius);
    if ~any(silent_keep)
        finite_idx = find(all(isfinite(silent_points), 2));
        if ~isempty(finite_idx)
            [~, nearest_ord] = sort(sqrt(sum((silent_points(finite_idx, :) - active_centroid) .^ 2, 2)), 'ascend');
            silent_keep(finite_idx(nearest_ord(1:min(2, numel(nearest_ord))))) = true;
        end
    end
    silent_sat_idx = silent_sat_idx(silent_keep);
    silent_points = silent_points(silent_keep, :);
    silent_proj_on_traj = silent_proj_on_traj(:, silent_keep, :);
end

[traj_window_x, traj_window_y, traj_window_time_s] = build_feasible_trajectory_local( ...
    active_contact_points, active_support_centers_s, cfg);
traj_pts = [traj_window_x(:), traj_window_y(:)];
traj_window_conf = ones(size(traj_window_x(:)));

active_match_time_s = active_support_centers_s(:);
active_match_dist_m = zeros(numel(active_support_centers_s), 1);

silent_match_time_s = nan(numel(silent_sat_idx), 1);
silent_match_dist_m = nan(numel(silent_sat_idx), 1);
for i = 1:numel(silent_sat_idx)
    if isempty(traj_pts)
        break;
    end
    d_i = sqrt(sum((traj_pts - silent_points(i, :)) .^ 2, 2));
    [silent_match_dist_m(i), idx_min] = min(d_i, [], 'omitnan');
    if ~isempty(idx_min) && isfinite(idx_min)
        silent_match_time_s(i) = traj_window_time_s(idx_min);
    end
end

geometry_demo.active_sat_ids = step1.valid_sats(active_sat_idx);
geometry_demo.silent_sat_ids = step1.valid_sats(silent_sat_idx);
geometry_demo.active_traces = active_traces;
geometry_demo.active_mid_points = active_contact_points;
geometry_demo.silent_points = silent_points;
geometry_demo.traj_x = traj_window_x(:);
geometry_demo.traj_y = traj_window_y(:);
geometry_demo.traj_conf = traj_window_conf(:);
geometry_demo.traj_time_s = traj_window_time_s(:);
geometry_demo.event_segments = event_segments;
geometry_demo.active_match_time_s = active_match_time_s;
geometry_demo.active_match_dist_m = active_match_dist_m;
geometry_demo.silent_match_time_s = silent_match_time_s;
geometry_demo.silent_match_dist_m = silent_match_dist_m;
geometry_demo.window_time_s = seconds(step1.t_grid(window_idx) - step1.t_grid(window_idx(1)));
geometry_demo.interaction_half_span = cfg.geometry_interaction_half_span;
geometry_demo.effective_radius_m = cfg.geometry_effective_radius_m;
geometry_demo.sample_epoch_count = numel(sample_epoch_idx);
geometry_demo.sample_obs_epochs = sample_epoch_idx(:);
geometry_demo.active_support_centers_s = active_support_centers_s(:);
end

function geo_ref = estimate_geometry_reference_local(obs_data, nav_data, epoch_list)
geo_ref = struct('ok', false, 'rec_pos_mean', [NaN, NaN, NaN], 'lat_deg', NaN, 'lon_deg', NaN, 'alt_m', NaN);

if isempty(epoch_list)
    return;
end

rec_buf = nan(numel(epoch_list), 3);
count = 0;
for i = 1:numel(epoch_list)
    try
        [rec_pos, ~, ~, lat_deg, lon_deg, alt_m] = calculate_receiver_position(obs_data, nav_data, epoch_list(i));
    catch
        continue;
    end
    if all(isfinite(rec_pos)) && all(isfinite([lat_deg, lon_deg, alt_m]))
        count = count + 1;
        rec_buf(count, :) = rec_pos(:)';
        geo_ref.lat_deg = lat_deg;
        geo_ref.lon_deg = lon_deg;
        geo_ref.alt_m = alt_m;
    end
end

if count < 1
    return;
end

geo_ref.ok = true;
geo_ref.rec_pos_mean = mean(rec_buf(1:count, :), 1, 'omitnan');
end

function proj_xy = compute_epoch_projection_local(obs_data, nav_data, epoch_idx, valid_sats, geo_ref, cfg)
proj_xy = nan(numel(valid_sats), 2);

try
    [rec_pos, ~, sat_states, lat_deg, lon_deg, alt_m] = calculate_receiver_position(obs_data, nav_data, epoch_idx);
catch
    return;
end

if ~all(isfinite(rec_pos))
    return;
end

if geo_ref.ok
    ref_pos = geo_ref.rec_pos_mean(:);
    ref_lat = geo_ref.lat_deg;
    ref_lon = geo_ref.lon_deg;
    ref_alt = geo_ref.alt_m;
else
    ref_pos = rec_pos(:);
    ref_lat = lat_deg;
    ref_lon = lon_deg;
    ref_alt = alt_m;
end

for s = 1:numel(valid_sats)
    sid = valid_sats{s};
    if ~isfield(sat_states, sid)
        continue;
    end

    sat_p = sat_states.(sid).position(:);
    [e, n, u] = ecef2enu( ...
        sat_p(1) - ref_pos(1), ...
        sat_p(2) - ref_pos(2), ...
        sat_p(3) - ref_pos(3), ...
        ref_lat, ref_lon, ref_alt);

    dist = norm([e, n, u]);
    if dist <= 0 || u <= 0
        continue;
    end

    elev_deg = asind(u / dist);
    if elev_deg < cfg.geometry_min_elevation_deg
        continue;
    end

    scale = cfg.geometry_gesture_height_m / u;
    px = scale * e;
    py = scale * n;

    if hypot(px, py) > cfg.geometry_max_projection_radius
        continue;
    end

    proj_xy(s, :) = [px, py];
end
end

function fig = plot_signal_preprocessing_local(demo, cfg)
sat_idx = 1;
sat_color = [0.05 0.47 0.80];
baseline_i = round(demo.top_sat_baselines(sat_idx));
y_lo = floor(min([demo.top_sat_raw(:, sat_idx); demo.top_sat_clean(:, sat_idx); baseline_i])) - 1;
y_hi = ceil(max([demo.top_sat_raw(:, sat_idx); demo.top_sat_clean(:, sat_idx); baseline_i])) + 1;
x = demo.sample_idx;
[x_raw, y_raw] = densify_integer_series_local(x, demo.top_sat_raw(:, sat_idx), 40);
[x_clean, y_clean] = densify_curve_local(x, demo.top_sat_clean(:, sat_idx), 28);

fig = figure('Color', 'w', 'Position', [140 120 1380 430]);
ax = axes(fig);
hold(ax, 'on');
stairs(ax, x_raw, y_raw, '-', 'Color', [0.72 0.72 0.72], 'LineWidth', 1.0);
plot(ax, x_clean, y_clean, '-', 'Color', sat_color, 'LineWidth', 2.1);
yline(ax, baseline_i, '--', 'Color', [0.45 0.45 0.45], 'LineWidth', 1.1);
shade_segments_local(ax, x, demo.display_segments_rel);
grid(ax, 'on');
xlabel(ax, 'Sample point');
ylabel(ax, sprintf('%s\nC/N0', demo.top_sat_ids{sat_idx}), 'Interpreter', 'none');
ylim(ax, [y_lo y_hi]);
yticks(ax, y_lo:y_hi);
set_sample_axis_ticks_local(ax, x, cfg.sample_rate_hz);
apply_axis_typography_local(ax);
end

function fig = plot_diffraction_feature_extraction_local(demo, cfg)
fd = demo.feature_demo;
x = fd.sample_idx;
wave_color = [0.05 0.47 0.80];
tick_font_size = 30;
label_font_size = 30;
annotation_font_size = 30;
[x_raw, y_raw] = densify_integer_series_local(x, fd.raw_wave, 260);
[x_clean, y_clean] = densify_curve_local(x, fd.clean_wave, 120);
fig = figure('Color', 'w', 'Position', [140 120 1740 460]);
ax1 = axes(fig);
hold(ax1, 'on');
h_raw = stairs(ax1, x_raw, y_raw, '-', 'Color', [0.78 0.78 0.78], 'LineWidth', 1.1);
h_clean = plot(ax1, x_clean, y_clean, '-', 'Color', wave_color, 'LineWidth', 2.1);
yline(ax1, fd.baseline, '--', 'Color', [0.25 0.25 0.25], 'LineWidth', 1.0, ...
    'HandleVisibility', 'off');
shade_segments_local(ax1, x, fd.event_segments);
for i = 1:numel(fd.event_segments)
    seg_i = fd.event_segments(i);
    x_i = fd.event_center_idx(i);
    y_i = min(fd.clean_wave(seg_i.start_idx:seg_i.end_idx)) - 0.12;
    text(ax1, x_i, y_i, fd.event_labels{i}, ...
        'Color', [0.18 0.47 0.80], ...
        'FontName', plot_font_name_local(), ...
        'FontSize', annotation_font_size, ...
        'FontWeight', 'bold', ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'top', ...
        'Interpreter', 'tex');
end
grid(ax1, 'on');
ylabel(ax1, 'C/N0 (dB-Hz)');
xlabel(ax1, 'Sample point');
ylim(ax1, [30 50]);
yticks(ax1, [30 40 45 50]);
set_sample_axis_ticks_local(ax1, x, cfg.sample_rate_hz);
lgd = legend(ax1, [h_raw, h_clean], ...
    {'Raw diffraction waveform', 'Smoothed diffraction envelope'}, ...
    'Location', 'northoutside', ...
    'Orientation', 'horizontal');
set(lgd, 'FontName', plot_font_name_local(), 'FontSize', tick_font_size);
apply_axis_typography_local(ax1, tick_font_size, label_font_size);
compact_axes_horizontal_local(ax1, [0.006 0.006]);
end

function fig = plot_diffraction_event_features_local(demo, cfg)
fd = demo.feature_demo;
x = 1:numel(fd.event_labels);
tick_font_size = 30;
label_font_size = 30;
fig = figure('Color', 'w', 'Position', [140 120 1740 460]);
ax = axes(fig);
hold(ax, 'on');
feat_bar = bar(ax, x, fd.feature_mat, 'grouped', 'BarWidth', 0.80);
feat_bar(1).FaceColor = [0.18 0.47 0.80];
feat_bar(2).FaceColor = [0.92 0.56 0.16];
feat_bar(3).FaceColor = [0.20 0.60 0.32];
feat_bar(1).EdgeColor = 'none';
feat_bar(2).EdgeColor = 'none';
feat_bar(3).EdgeColor = 'none';
feat_bar(1).FaceAlpha = 1;
feat_bar(2).FaceAlpha = 1;
feat_bar(3).FaceAlpha = 1;
ylim(ax, [0 1.05]);
grid(ax, 'on');
set(ax, 'Layer', 'bottom');
ylabel(ax, {'Normalized', 'feature value'});
xticks(ax, x);
xticklabels(ax, fd.event_labels);
xlabel(ax, 'Diffraction event');
xlim(ax, [0.45, numel(fd.event_labels) + 0.55]);
lgd = legend(ax, feat_bar, {'Amplitude', 'Oscillation', 'Continuity'}, ...
    'Location', 'northwest');
set(lgd, 'FontName', plot_font_name_local(), 'FontSize', tick_font_size);
apply_axis_typography_local(ax, tick_font_size, label_font_size);
compact_axes_to_figure_local(ax, [0.025 0.010 0.010 0.010]);
end

function fig = plot_diffraction_saliency_scores_local(demo, cfg)
fd = demo.feature_demo;
x = 1:numel(fd.event_labels);
tick_font_size = 30;
label_font_size = 30;
fig = figure('Color', 'w', 'Position', [140 120 1740 460]);
ax = axes(fig);
hold(ax, 'on');
bar(ax, x, fd.saliency_norm, 0.68, ...
    'FaceColor', [0.46 0.29 0.68], 'EdgeColor', 'none', 'FaceAlpha', 1);
xticks(ax, x);
xticklabels(ax, fd.event_labels);
grid(ax, 'on');
set(ax, 'Layer', 'bottom');
ylim(ax, [0 1.10]);
xlim(ax, [0.45, numel(fd.event_labels) + 0.55]);
xlabel(ax, 'Diffraction event');
ylabel(ax, {'Normalized', 'fused score'});
yline(ax, 0.70, '--', 'Color', [0.45 0.45 0.45], 'LineWidth', 1.0);
apply_axis_typography_local(ax, tick_font_size, label_font_size);
compact_axes_to_figure_local(ax, [0.025 0.010 0.010 0.010]);
end

function fig = plot_geometric_inversion_local(demo, cfg)
gd = demo.geometry_demo;
colors = lines(max(5, numel(gd.active_sat_ids) + numel(gd.silent_sat_ids)));
tick_font_size = 30;
label_font_size = 30;
annotation_font_size = 30;

fig = figure('Color', 'w', 'Position', [160 100 760 760]);
ax1 = axes(fig);
hold(ax1, 'on');
rectangle(ax1, 'Position', [-gd.interaction_half_span, -gd.interaction_half_span, ...
    2 * gd.interaction_half_span, 2 * gd.interaction_half_span], ...
    'EdgeColor', [0.25 0.25 0.25], 'LineStyle', '--', 'LineWidth', 1.2);

for i = 1:numel(gd.active_sat_ids)
    pts = squeeze(gd.active_traces(:, i, :));
    keep = all(isfinite(pts), 2);
    if any(keep)
        plot(ax1, pts(keep, 1), pts(keep, 2), ':', 'Color', colors(i, :), 'LineWidth', 1.2);
        lbl_idx = find(keep, min(cfg.geometry_contact_label_frames, nnz(keep)), 'last');
        lbl_pt = mean(pts(lbl_idx, :), 1, 'omitnan');
        if all(isfinite(lbl_pt))
            text(ax1, lbl_pt(1), lbl_pt(2), gd.active_sat_ids{i}, 'Color', colors(i, :), ...
                'FontName', plot_font_name_local(), ...
                'FontSize', annotation_font_size, ...
                'FontWeight', 'bold', 'Interpreter', 'none', 'HorizontalAlignment', 'center');
        end
    end
    mid_pt = gd.active_mid_points(i, :);
    if all(isfinite(mid_pt))
        theta = linspace(0, 2 * pi, 120);
        plot(ax1, mid_pt(1) + gd.effective_radius_m * cos(theta), ...
            mid_pt(2) + gd.effective_radius_m * sin(theta), '-', ...
            'Color', colors(i, :), 'LineWidth', 1.4);
        scatter(ax1, mid_pt(1), mid_pt(2), 82, colors(i, :), 'filled', ...
            'MarkerEdgeColor', [0.12 0.12 0.12]);
    end
end

if ~isempty(gd.silent_points)
    scatter(ax1, gd.silent_points(:, 1), gd.silent_points(:, 2), 42, ...
        repmat([0.60 0.60 0.60], size(gd.silent_points, 1), 1), 'x', 'LineWidth', 1.4);
    for i = 1:size(gd.silent_points, 1)
        if all(isfinite(gd.silent_points(i, :)))
            theta = linspace(0, 2 * pi, 120);
            plot(ax1, gd.silent_points(i, 1) + gd.effective_radius_m * cos(theta), ...
                gd.silent_points(i, 2) + gd.effective_radius_m * sin(theta), '--', ...
                'Color', [0.55 0.55 0.55], 'LineWidth', 1.0);
            text(ax1, gd.silent_points(i, 1), gd.silent_points(i, 2), ['  ', gd.silent_sat_ids{i}], ...
                'Color', [0.45 0.45 0.45], 'Interpreter', 'none', ...
                'FontName', plot_font_name_local(), ...
                'FontSize', annotation_font_size);
        end
    end
end

if ~isempty(gd.traj_x)
    plot(ax1, gd.traj_x, gd.traj_y, '-', 'Color', [0.08 0.08 0.08], 'LineWidth', 2.8);
    scatter(ax1, gd.traj_x, gd.traj_y, 12, [0.10 0.10 0.10], 'filled', 'MarkerFaceAlpha', 0.18);
    plot(ax1, gd.traj_x(1), gd.traj_y(1), 'o', 'MarkerFaceColor', [0.10 0.10 0.10], ...
        'MarkerEdgeColor', 'k', 'MarkerSize', 6);
    plot(ax1, gd.traj_x(end), gd.traj_y(end), 's', 'MarkerFaceColor', [0.10 0.10 0.10], ...
        'MarkerEdgeColor', 'k', 'MarkerSize', 6);
end

plot(ax1, 0, 0, 'kp', 'MarkerFaceColor', [0.95 0.78 0.18], 'MarkerSize', 11);
grid(ax1, 'on');
axis(ax1, 'equal');
xlim(ax1, [-0.55 0.55]);
ylim(ax1, [-0.55 0.55]);
xlabel(ax1, 'East on gesture plane (m)');
ylabel(ax1, 'North on gesture plane (m)');
apply_axis_typography_local(ax1, tick_font_size, label_font_size);
compact_axes_to_figure_local(ax1, [0.006 0.006 0.006 0.006]);
end

function fig = plot_geometric_constraint_rows_local(demo, cfg)
gd = demo.geometry_demo;
colors = lines(max(5, numel(gd.active_sat_ids) + numel(gd.silent_sat_ids)));
tick_font_size = 30;
label_font_size = 30;
segment_font_size = 30;

fig = figure('Color', 'w', 'Position', [160 100 760 760]);
ax2 = axes(fig);
hold(ax2, 'on');
if isfield(gd, 'event_segments_sample')
    x_rows = gd.event_segments_sample;
    x_active = gd.active_match_sample_idx;
    x_silent = gd.silent_match_sample_idx;
    x_window = gd.window_sample_idx;
else
    x_rows = gd.event_segments;
    x_active = gd.active_match_time_s;
    x_silent = gd.silent_match_time_s;
    x_window = gd.window_time_s;
end
xlim(ax2, [min(x_window) max(x_window)]);
active_order = 1:numel(gd.active_sat_ids);
finite_mask = isfinite(x_active(active_order));
active_order_finite = active_order(finite_mask);
if ~isempty(active_order_finite)
    [~, sort_idx] = sort(x_active(active_order_finite), 'ascend');
    active_order = [active_order_finite(sort_idx), active_order(~finite_mask)];
end

for row_idx = 1:numel(active_order)
    src_idx = active_order(row_idx);
    row_segments = x_rows{src_idx};
    if isstruct(row_segments)
        row_segments = row_segments(:).';
    end
    amp_i = 0.34 + 0.04 * mod(src_idx + 1, 3);
    if any(src_idx == 3)
        amp_i = amp_i + 0.06;
    end
    [x_dense_i, y_wave_i] = build_constraint_row_waveform_local( ...
        x_window, row_segments, x_active(src_idx), row_idx, amp_i, src_idx, false);
    plot(ax2, x_dense_i, row_idx * ones(size(x_dense_i)), '-', ...
        'Color', min(0.90, colors(src_idx, :) + 0.34), 'LineWidth', 0.8);
    annotate_constraint_action_segments_local( ...
        ax2, row_segments, x_active(src_idx), row_idx, colors(src_idx, :), ...
        gd.active_sat_ids{src_idx}, segment_font_size);
    plot(ax2, x_dense_i, y_wave_i, '-', 'Color', colors(src_idx, :), 'LineWidth', 2.0);
end
base_y = numel(active_order);
for i = 1:numel(gd.silent_sat_ids)
    y_i = base_y + i;
    [x_dense_i, y_wave_i] = build_constraint_row_waveform_local( ...
        x_window, struct([]), x_silent(i), y_i, 0.06, i, true);
    plot(ax2, x_dense_i, y_i * ones(size(x_dense_i)), '--', ...
        'Color', [0.82 0.82 0.82], 'LineWidth', 0.8);
    plot(ax2, x_dense_i, y_wave_i, '-', 'Color', [0.48 0.48 0.48], 'LineWidth', 1.2);
end
grid(ax2, 'on');
active_labels = gd.active_sat_ids(active_order);
all_labels = [active_labels(:); gd.silent_sat_ids(:)];
yticks(ax2, 1:numel(all_labels));
yticklabels(ax2, all_labels);
xlabel(ax2, 'Sample point');
ylabel(ax2, 'Constraint rows');
ylim(ax2, [0.4 max(1.8, numel(all_labels) + 0.8)]);
xlim(ax2, [min(x_window) max(x_window)]);
set_sample_axis_ticks_local(ax2, x_window, cfg.sample_rate_hz);
apply_axis_typography_local(ax2, tick_font_size, label_font_size);
compact_axes_to_figure_local(ax2, [0.006 0.006 0.006 0.006]);
end

function annotate_constraint_action_segments_local(ax, event_segments, match_pos, row_level, color_rgb, sat_label, font_size)
if isempty(event_segments) || ~isstruct(event_segments)
    return;
end

event_segments = event_segments(:).';
for j = 1:numel(event_segments)
    seg_j = event_segments(j);
    if isfield(seg_j, 'start_idx')
        seg_start = seg_j.start_idx;
        seg_end = seg_j.end_idx;
    else
        seg_start = seg_j.t_start;
        seg_end = seg_j.t_end;
    end
    if ~isfinite(seg_start) || ~isfinite(seg_end)
        continue;
    end
    if seg_end < seg_start
        tmp = seg_start;
        seg_start = seg_end;
        seg_end = tmp;
    end

    y_low = row_level - 0.44;
    y_high = row_level + 0.34;
    patch(ax, [seg_start seg_end seg_end seg_start], [y_low y_low y_high y_high], color_rgb, ...
        'FaceAlpha', 0.10, 'EdgeColor', color_rgb, 'LineStyle', ':', 'LineWidth', 1.0);

    if isfinite(match_pos)
        line(ax, [match_pos match_pos], [y_low y_high], ...
            'Color', color_rgb, 'LineStyle', '--', 'LineWidth', 1.1);
    end

    x_mid = 0.5 * (seg_start + seg_end);
    if isfinite(match_pos)
        x_mid = 0.45 * x_mid + 0.55 * match_pos;
    end
    x_limits = xlim(ax);
    x_pad = 0.24 * diff(x_limits);
    x_mid = min(max(x_mid, x_limits(1) + x_pad), x_limits(2) - x_pad);
    text(ax, x_mid, row_level + 0.37, sprintf('%s action', sat_label), ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
        'Color', color_rgb, 'FontName', plot_font_name_local(), ...
        'FontSize', font_size, 'FontWeight', 'bold', ...
        'BackgroundColor', 'w', 'Margin', 1.0, 'Interpreter', 'none');
end
end

function [x_dense, y_wave] = build_constraint_row_waveform_local(x_window, event_segments, match_pos, row_level, amp, style_id, is_silent)
x_window = x_window(:);
if isempty(x_window)
    x_dense = x_window;
    y_wave = x_window;
    return;
end

dense_factor = 14;
x_dense = linspace(x_window(1), x_window(end), dense_factor * max(numel(x_window) - 1, 1) + 1).';
y_wave = row_level * ones(size(x_dense));

ctrl_count = max(10, ceil(numel(x_dense) / 120));
ctrl_idx = unique(round(linspace(1, numel(x_dense), ctrl_count)).');
ctrl_wave = 0.55 * sin(0.18 * ctrl_idx + 0.70 * style_id) ...
    + 0.45 * cos(0.07 * ctrl_idx + 0.35 * style_id);
noise_amp = 0.010 + 0.003 * mod(style_id, 3);
if is_silent
    noise_amp = 0.018;
end
noise_dense = interp1(ctrl_idx, ctrl_wave, (1:numel(x_dense)).', 'pchip');
noise_dense = noise_amp * noise_dense / max(max(abs(noise_dense)), 1e-6);
y_wave = y_wave + noise_dense;

if isempty(event_segments) || ~isfinite(match_pos) || is_silent
    return;
end

if ~isstruct(event_segments)
    return;
end

for j = 1:numel(event_segments)
    seg_j = event_segments(j);
    if isfield(seg_j, 'start_idx')
        seg_start = seg_j.start_idx;
        seg_end = seg_j.end_idx;
    else
        seg_start = seg_j.t_start;
        seg_end = seg_j.t_end;
    end
    seg_len = max(8, seg_end - seg_start + 1);
    seg_center = 0.5 * (seg_start + seg_end);
    if isfinite(match_pos)
        seg_center = 0.60 * match_pos + 0.40 * seg_center;
    end
    profile_scale = 0.55 * seg_len;
    u = (x_dense - seg_center) / max(profile_scale, 1e-3);
    profile_j = diffraction_profile_local(u, style_id + j - 1);
    profile_j = profile_j / max(max(profile_j), 1e-6);
    taper_j = exp(-0.5 * ((x_dense - seg_center) / max(0.70 * seg_len, 1e-3)).^2);
    y_wave = y_wave - amp * profile_j .* taper_j;
end
end

function compact_axes_to_figure_local(ax, pad)
if nargin < 2 || isempty(pad)
    pad = [0.01 0.01 0.01 0.01];
end
pad = pad(:).';
if numel(pad) == 1
    pad = repmat(pad, 1, 4);
end
drawnow;
set(ax, 'Units', 'normalized');
tight_inset = get(ax, 'TightInset');
outer_box = [pad(1), pad(2), 1 - pad(1) - pad(3), 1 - pad(2) - pad(4)];
pos = [ ...
    outer_box(1) + tight_inset(1), ...
    outer_box(2) + tight_inset(2), ...
    outer_box(3) - tight_inset(1) - tight_inset(3), ...
    outer_box(4) - tight_inset(2) - tight_inset(4)];
pos(3) = max(pos(3), 0.05);
pos(4) = max(pos(4), 0.05);
set(ax, 'Position', pos, 'LooseInset', tight_inset);
end

function compact_axes_horizontal_local(ax, pad_lr)
if nargin < 2 || isempty(pad_lr)
    pad_lr = [0.01 0.01];
end
pad_lr = pad_lr(:).';
if numel(pad_lr) == 1
    pad_lr = repmat(pad_lr, 1, 2);
end
drawnow;
set(ax, 'Units', 'normalized');
tight_inset = get(ax, 'TightInset');
pos = get(ax, 'Position');
new_left = pad_lr(1) + tight_inset(1);
new_width = 1 - pad_lr(1) - pad_lr(2) - tight_inset(1) - tight_inset(3);
pos(1) = new_left;
pos(3) = max(new_width, 0.05);
set(ax, 'Position', pos, 'LooseInset', tight_inset);
end

function postprocess_pdf_groups_local(manifest_tbl)
script_path = fullfile(fileparts(mfilename('fullpath')), 'normalize_pdf_group_pages.py');
if exist(script_path, 'file') ~= 2
    warning('PDF post-processing script not found: %s', script_path);
    return;
end

group_ids = { ...
    ["B1_diffraction_waveforms", "B2_diffraction_feature_values", "B3_diffraction_saliency_scores"], ...
    ["C1_gesture_plane_geometric_inversion", "C2_geometric_constraint_rows"]};

for group_idx = 1:numel(group_ids)
    ids = group_ids{group_idx};
    pdf_paths = strings(1, numel(ids));
    for i = 1:numel(ids)
        match_idx = find(manifest_tbl.figure_id == ids(i), 1);
        if isempty(match_idx)
            warning('Missing figure id for PDF group normalization: %s', ids(i));
            pdf_paths = strings(0);
            break;
        end
        pdf_paths(i) = manifest_tbl.pdf_path(match_idx);
    end
    if isempty(pdf_paths)
        continue;
    end
    cmd_parts = ["python", quote_cmd_arg_local(script_path), "--pad", "8"];
    for i = 1:numel(pdf_paths)
        cmd_parts(end + 1) = quote_cmd_arg_local(pdf_paths(i)); %#ok<AGROW>
    end
    [status, cmd_out] = system(strjoin(cmd_parts, " "));
    if status ~= 0
        warning('PDF group normalization failed for group %d:\n%s', group_idx, cmd_out);
    end
end
end

function arg = quote_cmd_arg_local(value)
arg = """" + string(value) + """";
end

function row = save_figure_bundle_local(fig, figure_id, fmt_dirs, cfg)
png_path = fullfile(fmt_dirs.png, figure_id + ".png");
pdf_path = fullfile(fmt_dirs.pdf, figure_id + ".pdf");
fig_path = fullfile(fmt_dirs.fig, figure_id + ".fig");

if ~isgraphics(fig)
    error('Invalid figure handle for %s.', figure_id);
end

fixed_canvas_ids = [ ...
    "B1_diffraction_waveforms", ...
    "B2_diffraction_feature_values", ...
    "B3_diffraction_saliency_scores", ...
    "C1_gesture_plane_geometric_inversion", ...
    "C2_geometric_constraint_rows"];
use_fixed_canvas = any(strcmp(string(figure_id), fixed_canvas_ids));

if use_fixed_canvas
    configure_fixed_canvas_output_local(fig, cfg);
else
    set(fig, 'PaperPositionMode', 'auto');
end
suppress_axes_toolbar_local(fig);
force_figure_font_local(fig);
savefig(fig, char(fig_path));
drawnow;
if use_fixed_canvas
    try
        print(fig, char(png_path), '-dpng', ['-r', num2str(cfg.figure_resolution)]);
    catch
        saveas(fig, char(png_path));
    end
    try
        exportgraphics(fig, char(pdf_path), 'ContentType', 'vector');
    catch
        pdf_resolution = max(1200, cfg.figure_resolution);
        print(fig, char(pdf_path), '-dpdf', '-painters', ['-r', num2str(pdf_resolution)]);
    end
else
    try
        exportgraphics(fig, char(png_path), 'Resolution', cfg.figure_resolution);
    catch
        saveas(fig, char(png_path));
    end
    try
        exportgraphics(fig, char(pdf_path), 'ContentType', 'vector');
    catch
        set(0, 'CurrentFigure', fig);
        print(char(pdf_path), '-dpdf', '-painters', '-bestfit');
    end
end
if isgraphics(fig)
    close(fig);
end

row = struct( ...
    'figure_id', string(figure_id), ...
    'png_path', string(png_path), ...
    'pdf_path', string(pdf_path), ...
    'fig_path', string(fig_path));
end

function force_figure_font_local(fig)
font_name = plot_font_name_local();
font_objects = findall(fig, '-property', 'FontName');
for i = 1:numel(font_objects)
    try
        set(font_objects(i), 'FontName', font_name);
    catch
    end
end

try
    set(fig, ...
        'DefaultAxesFontName', font_name, ...
        'DefaultTextFontName', font_name, ...
        'DefaultLegendFontName', font_name, ...
        'DefaultColorbarFontName', font_name);
catch
end
end

function suppress_axes_toolbar_local(fig)
axes_handles = findall(fig, 'Type', 'axes');
for i = 1:numel(axes_handles)
    try
        axes_handles(i).Toolbar.Visible = 'off';
    catch
    end
    try
        disableDefaultInteractivity(axes_handles(i));
    catch
    end
end
end

function configure_fixed_canvas_output_local(fig, cfg)
set(fig, 'Units', 'pixels', 'InvertHardcopy', 'off');
fig_pos = get(fig, 'Position');
fig_w_px = max(1, fig_pos(3));
fig_h_px = max(1, fig_pos(4));
target_dpi = get(groot, 'ScreenPixelsPerInch');
if ~isfinite(target_dpi) || target_dpi <= 0
    target_dpi = 96;
end
paper_w_in = fig_w_px / target_dpi;
paper_h_in = fig_h_px / target_dpi;

set(fig, ...
    'PaperUnits', 'inches', ...
    'PaperPositionMode', 'manual', ...
    'PaperPosition', [0 0 paper_w_in paper_h_in], ...
    'PaperSize', [paper_w_in paper_h_in]);
try
    set(fig, 'PaperType', '<custom>');
catch
end
end

function baselines = compute_baselines_local(mat)
baselines = nan(1, size(mat, 2));
for i = 1:size(mat, 2)
    baselines(i) = compute_single_baseline_local(mat(:, i));
end
end

function baseline = compute_single_baseline_local(col)
valid = col(isfinite(col));
if isempty(valid)
    baseline = NaN;
else
    baseline = mode(round(valid));
end
end

function y = movmean_local(x, win)
if isempty(x)
    y = x;
    return;
end
win = max(1, round(win));
y = movmean(x, win, 'omitnan');
end

function x_n = normalize_nonnegative_local(x)
x = x(:);
x(~isfinite(x)) = 0;
x = max(x, 0);
mx = max(x);
if ~isfinite(mx) || mx <= 0
    x_n = zeros(size(x));
else
    x_n = x / mx;
end
end

function mask_out = fill_short_gaps_local(mask_in, max_gap)
mask_out = logical(mask_in(:));
if isempty(mask_out)
    return;
end

gap_mask = ~mask_out;
edges = diff([false; gap_mask; false]);
gap_starts = find(edges == 1);
gap_ends = find(edges == -1) - 1;
for i = 1:numel(gap_starts)
    if gap_starts(i) == 1 || gap_ends(i) == numel(mask_out)
        continue;
    end
    if (gap_ends(i) - gap_starts(i) + 1) <= max_gap
        mask_out(gap_starts(i):gap_ends(i)) = true;
    end
end
end

function x_n = normalize_by_mask_local(x, mask)
x = x(:);
mask = logical(mask(:));
x_n = zeros(size(x));
if isempty(x)
    return;
end
ref = x(mask & isfinite(x));
ref = ref(ref >= 0);
if isempty(ref)
    ref = x(isfinite(x));
    ref = ref(ref >= 0);
end
mx = max(ref);
if ~isfinite(mx) || mx <= 0
    return;
end
x_n = clamp01_local(max(x, 0) / mx);
end

function y = clamp01_local(x)
y = min(max(x, 0), 1);
end

function segment_rel = refine_feature_segment_local(deviation_norm, oscillation_norm, continuity_norm, coarse_mask, coarse_seg, cfg)
segment_rel = coarse_seg;
coarse_idx = find(coarse_mask);
if isempty(coarse_idx)
    return;
end

support_trace = 0.55 * deviation_norm(:) + 0.30 * oscillation_norm(:) + 0.15 * continuity_norm(:);
support_trace(~coarse_mask) = 0;
coarse_support = support_trace(coarse_mask);
peak_score = max(coarse_support, [], 'omitnan');
if ~isfinite(peak_score) || peak_score <= 0
    return;
end

[~, peak_pos_local] = max(coarse_support);
peak_idx = coarse_idx(peak_pos_local);
thr = max(0.42, 0.72 * peak_score);
candidate_mask = coarse_mask & (support_trace >= thr);
candidate_mask = fill_short_gaps_local(candidate_mask, cfg.feature_demo_gap_fill);
segs = mask_to_segments_local(candidate_mask);

if isempty(segs)
    half_span = floor(cfg.feature_demo_min_pts / 2);
    segment_rel.start_idx = max(coarse_idx(1), peak_idx - half_span);
    segment_rel.end_idx = min(coarse_idx(end), peak_idx + half_span);
    return;
end

pick_idx = [];
for i = 1:numel(segs)
    if segs(i).start_idx <= peak_idx && peak_idx <= segs(i).end_idx
        pick_idx = i;
        break;
    end
end
if isempty(pick_idx)
    seg_score = zeros(numel(segs), 1);
    for i = 1:numel(segs)
        idx_i = segs(i).start_idx:segs(i).end_idx;
        seg_score(i) = mean(support_trace(idx_i), 'omitnan');
    end
    [~, pick_idx] = max(seg_score);
end

segment_rel = segs(pick_idx);
seg_len = segment_rel.end_idx - segment_rel.start_idx + 1;
if seg_len < cfg.feature_demo_min_pts
    need = cfg.feature_demo_min_pts - seg_len;
    left_pad = floor(need / 2);
    right_pad = ceil(need / 2);
    segment_rel.start_idx = max(coarse_idx(1), segment_rel.start_idx - left_pad);
    segment_rel.end_idx = min(coarse_idx(end), segment_rel.end_idx + right_pad);
end
end

function crop_idx = build_feature_crop_idx_local(t_s, seg, pad_sec)
if isempty(t_s)
    crop_idx = zeros(0, 1);
    return;
end
if numel(t_s) == 1
    crop_idx = 1;
    return;
end
dt = median(diff(t_s), 'omitnan');
if ~isfinite(dt) || dt <= 0
    dt = 1;
end
pad_pts = max(3, round(pad_sec / dt));
crop_idx = (max(1, seg.start_idx - pad_pts):min(numel(t_s), seg.end_idx + pad_pts)).';
end

function pt = representative_contact_point_local(trace_xy, sample_time_s, ref_time_s)
pt = [NaN, NaN];
if isempty(trace_xy)
    return;
end
trace_xy = reshape(trace_xy, [], 2);
finite_mask = all(isfinite(trace_xy), 2);
if ~any(finite_mask)
    return;
end

if nargin >= 3 && isfinite(ref_time_s)
    valid_idx = find(finite_mask);
    [~, ord] = min(abs(sample_time_s(valid_idx) - ref_time_s));
    pt = trace_xy(valid_idx(ord), :);
else
    pt = mean(trace_xy(finite_mask, :), 1, 'omitnan');
end
end

function [traj_x, traj_y, traj_t_s] = build_feasible_trajectory_local(contact_points, contact_times_s, cfg)
traj_x = [];
traj_y = [];
traj_t_s = [];
if isempty(contact_points)
    return;
end

keep = all(isfinite(contact_points), 2) & isfinite(contact_times_s(:));
contact_points = contact_points(keep, :);
contact_times_s = contact_times_s(keep);
if isempty(contact_points)
    return;
end

[contact_times_s, ord] = sort(contact_times_s, 'ascend');
contact_points = contact_points(ord, :);

unique_keep = true(size(contact_times_s));
for i = 2:numel(contact_times_s)
    if norm(contact_points(i, :) - contact_points(i - 1, :)) < 0.015
        unique_keep(i) = false;
    end
end
contact_points = contact_points(unique_keep, :);
contact_times_s = contact_times_s(unique_keep);

if size(contact_points, 1) == 1
    traj_x = contact_points(1, 1);
    traj_y = contact_points(1, 2);
    traj_t_s = contact_times_s(1);
    return;
end

for i = 1:(size(contact_points, 1) - 1)
    dt = max(0.3, contact_times_s(i + 1) - contact_times_s(i));
    n_link = max(8, round(cfg.geometry_path_link_samples * dt / max(0.8, range(contact_times_s) + eps)));
    tau = linspace(0, 1, n_link).';
    x_i = (1 - tau) * contact_points(i, 1) + tau * contact_points(i + 1, 1);
    y_i = (1 - tau) * contact_points(i, 2) + tau * contact_points(i + 1, 2);
    t_i = (1 - tau) * contact_times_s(i) + tau * contact_times_s(i + 1);
    if i > 1
        x_i = x_i(2:end);
        y_i = y_i(2:end);
        t_i = t_i(2:end);
    end
    traj_x = [traj_x; x_i]; %#ok<AGROW>
    traj_y = [traj_y; y_i]; %#ok<AGROW>
    traj_t_s = [traj_t_s; t_i]; %#ok<AGROW>
end

traj_x = movmean_local(traj_x, 3);
traj_y = movmean_local(traj_y, 3);
end

function segs = mask_to_segments_local(mask)
mask = logical(mask(:));
edges = diff([false; mask; false]);
starts = find(edges == 1);
ends = find(edges == -1) - 1;
segs = repmat(struct('start_idx', NaN, 'end_idx', NaN), numel(starts), 1);
for i = 1:numel(starts)
    segs(i).start_idx = starts(i);
    segs(i).end_idx = ends(i);
end
end

function sample_rate_hz = estimate_sampling_rate_local(t_grid)
dt = median(seconds(diff(t_grid)), 'omitnan');
if ~isfinite(dt) || dt <= 0
    sample_rate_hz = NaN;
else
    sample_rate_hz = 1 / dt;
end
end

function shade_segments_local(ax, t_s, segments)
if isempty(segments)
    return;
end
yl = ylim(ax);
for i = 1:numel(segments)
    t0 = t_s(segments(i).start_idx);
    t1 = t_s(segments(i).end_idx);
    patch(ax, [t0 t1 t1 t0], [yl(1) yl(1) yl(2) yl(2)], [0.18 0.47 0.80], ...
        'FaceAlpha', 0.08, 'EdgeColor', 'none');
end
uistack(findobj(ax, 'Type', 'patch'), 'bottom');
end

function dst = merge_cfg_local(dst, src)
keys = fieldnames(src);
for i = 1:numel(keys)
    key = keys{i};
    if isstruct(src.(key))
        if ~isfield(dst, key) || ~isstruct(dst.(key))
            dst.(key) = src.(key);
        else
            dst.(key) = merge_cfg_local(dst.(key), src.(key));
        end
    else
        dst.(key) = src.(key);
    end
end
end
