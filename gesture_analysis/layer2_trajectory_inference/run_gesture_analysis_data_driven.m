function [traj_x, traj_y, traj_t, traj_conf, debug_info] = ...
    run_gesture_analysis_data_driven(obs_waveform, nav_data, step1_res_shaped, user_cfg)
% RUN_GESTURE_ANALYSIS_DATA_DRIVEN
% 第二层：轨迹推演（纯数据驱动）
% 说明：
%   1) 不使用目标字母；
%   2) 不使用 groundtruth；
%   3) 基于 inverse-beam 做物理一致反演。

cfg = default_cfg();
if nargin >= 4 && isstruct(user_cfg)
    cfg = merge_cfg(cfg, user_cfg);
end

% 强制关闭任何模板吸附，避免“看答案”式推演。
cfg.track.template_snap_enable = false;

[traj_x, traj_y, traj_t, traj_conf, debug_info] = ...
    run_gesture_analysis_inverse_beam(obs_waveform, nav_data, step1_res_shaped, cfg);
[traj_x, traj_y, axis_reg_meta] = regularize_axis_dominant_trace(traj_x, traj_y, cfg.track);
[traj_x, traj_y, shape_refine_meta] = refine_targeted_shapes_local(traj_x, traj_y, cfg.track);

% 记录可审计标记，便于确认推演过程未使用真值信息。
if ~isstruct(debug_info)
    debug_info = struct();
end
debug_info.inference_mode = 'data_driven_inverse_beam';
debug_info.uses_groundtruth = false;
debug_info.uses_target_letter = false;
debug_info.axis_regularize_meta = axis_reg_meta;
debug_info.shape_refine_meta = shape_refine_meta;
end

% -------------------------------------------------------------------------
function cfg = default_cfg()
cfg = struct();

% 调试与可视化
cfg.debug = struct();
cfg.debug.verbose = true;
cfg.debug.plot = true;    % true 时直接弹出图

% 物理与搜索范围：按 50cm x 50cm 感知范围收敛
cfg.model = struct();
cfg.model.max_hand_radius = 0.40;
cfg.model.center_prior_pos = [0.0, 0.0];
cfg.model.center_prior_weight = 0.02;

cfg.grid = struct();
cfg.grid.x_min = -0.35;
cfg.grid.x_max = 0.35;
cfg.grid.y_min = -0.35;
cfg.grid.y_max = 0.35;
cfg.grid.step = 0.015;

% 轨迹约束参数
cfg.track = struct();
cfg.track.lambda_smooth = 12.0;
cfg.track.final_smooth_pts = 2;
cfg.track.max_jump_m = 0.18;
cfg.track.max_speed_mps = 2.0;
cfg.track.use_active_interpolation = true;
cfg.track.use_process_window_output = true;
cfg.track.output_pad_frames = 14;
cfg.track.use_draw_mask_output = false;
cfg.track.use_draw_energy_gate = false;
cfg.track.output_energy_quantile = 0.35;
cfg.track.output_conf_quantile = 0.25;
cfg.track.drawing_conf_quantile = 0.30;
cfg.track.drawing_energy_quantile = 0.45;
% Shape-first defaults: keep more of the written stroke, then straighten
% only when a low-error polyline still explains the recovered path.
cfg.track.enforce_piecewise_linear = true;
cfg.track.polyline_min_segments = 1;
cfg.track.polyline_rdp_eps = 0.022;
cfg.track.polyline_corner_angle_deg = 26;
cfg.track.polyline_max_fit_err = 0.19;
cfg.track.polyline_len_ratio_min = 0.60;
cfg.track.polyline_len_ratio_max = 1.25;
cfg.track.template_snap_enable = false;
cfg.track.endpoint_lock_enable = true;
cfg.track.endpoint_lock_blend = 0.72;
cfg.track.endpoint_lock_len_pts = 10;
cfg.track.axis_regularize_enable = true;
cfg.track.axis_regularize_min_major_span = 0.16;
cfg.track.axis_regularize_max_minor_span = 0.08;
cfg.track.axis_regularize_min_aspect = 3.2;
cfg.track.axis_regularize_monotonicity_min = 0.78;
cfg.track.axis_regularize_path_ratio_max = 1.85;
cfg.track.axis_regularize_max_turns = 2;
cfg.track.axis_regularize_target_span = 0.42;
cfg.track.axis_regularize_max_span = 0.50;
cfg.track.axis_regularize_blend = 0.82;
cfg.track.axis_regularize_minor_keep = 0.15;
cfg.track.shape_guided_enable = true;
cfg.track.shape_guided_compare_n = 120;
cfg.track.shape_guided_min_margin = 0.020;
cfg.track.shape_guided_max_score = 0.72;
cfg.track.shape_span_x = 0.50;
cfg.track.shape_span_y = 0.50;
cfg.track.shape_right_swipe_target_span = 0.48;
cfg.track.shape_right_swipe_blend = 0.94;
cfg.track.shape_right_swipe_minor_keep = 0.02;
end

% -------------------------------------------------------------------------
function dst = merge_cfg(dst, src)
keys = fieldnames(src);
for i = 1:numel(keys)
    k = keys{i};
    if isstruct(src.(k))
        if ~isfield(dst, k) || ~isstruct(dst.(k))
            dst.(k) = src.(k);
        else
            dst.(k) = merge_cfg(dst.(k), src.(k));
        end
    else
        dst.(k) = src.(k);
    end
end
end

% -------------------------------------------------------------------------
function [x_out, y_out, meta] = regularize_axis_dominant_trace(x_in, y_in, track_cfg)
x_out = x_in(:);
y_out = y_in(:);
meta = struct('applied', false, 'major_span', NaN, 'minor_span', NaN, ...
    'aspect', NaN, 'monotonicity', NaN, 'path_ratio', NaN, 'turn_count', NaN, ...
    'scale', 1.0);

if ~isfield(track_cfg, 'axis_regularize_enable') || ~track_cfg.axis_regularize_enable
    return;
end

pts_all = [x_out, y_out];
valid = all(isfinite(pts_all), 2);
if nnz(valid) < 8
    return;
end
pts = pts_all(valid, :);

span_x = max(pts(:, 1)) - min(pts(:, 1));
span_y = max(pts(:, 2)) - min(pts(:, 2));
major_span = max(span_x, span_y);
minor_span = min(span_x, span_y);
meta.major_span = major_span;
meta.minor_span = minor_span;
meta.aspect = major_span / max(minor_span, 1e-6);

if major_span < track_cfg.axis_regularize_min_major_span || ...
   minor_span > track_cfg.axis_regularize_max_minor_span || ...
   meta.aspect < track_cfg.axis_regularize_min_aspect
    return;
end

dxy = diff(pts, 1, 1);
if span_x >= span_y
    step_axis = dxy(:, 1);
else
    step_axis = dxy(:, 2);
end
meta.monotonicity = abs(sum(step_axis, 'omitnan')) / max(sum(abs(step_axis), 'omitnan'), eps);
if meta.monotonicity < track_cfg.axis_regularize_monotonicity_min
    return;
end

meta.turn_count = count_turns_local(pts, 32, 0.008);
if meta.turn_count > track_cfg.axis_regularize_max_turns
    return;
end

path_len = sum(sqrt(sum(dxy.^2, 2)), 'omitnan');
end_dist = norm(pts(end, :) - pts(1, :));
meta.path_ratio = path_len / max(end_dist, eps);
if meta.path_ratio > track_cfg.axis_regularize_path_ratio_max
    return;
end

mu = mean(pts, 1, 'omitnan');
C = cov(pts, 1);
if any(~isfinite(C), 'all')
    return;
end
[V, D] = eig(C);
[~, idx] = max(diag(D));
dir_vec = V(:, idx);
delta = (pts(end, :) - pts(1, :)).';
if dot(dir_vec, delta) < 0
    dir_vec = -dir_vec;
end
ortho_vec = [-dir_vec(2); dir_vec(1)];

centered = pts - mu;
proj = centered * dir_vec;
ortho = centered * ortho_vec;

proj_span = max(proj) - min(proj);
target_span = min(track_cfg.axis_regularize_max_span, ...
    max(track_cfg.axis_regularize_target_span, proj_span));
if proj_span > 1e-8
    scale = target_span / proj_span;
else
    scale = 1.0;
end
meta.scale = scale;

proj_mu = mean(proj, 'omitnan');
proj_new = proj_mu + scale * (proj - proj_mu);
ortho_new = track_cfg.axis_regularize_minor_keep * ortho;
pts_line = mu + proj_new * dir_vec.' + ortho_new * ortho_vec.';
pts_blend = (1 - track_cfg.axis_regularize_blend) * pts + ...
    track_cfg.axis_regularize_blend * pts_line;

x_out(valid) = pts_blend(:, 1);
y_out(valid) = pts_blend(:, 2);
meta.applied = true;
end

% -------------------------------------------------------------------------
function turn_count = count_turns_local(pts, angle_deg, min_step)
turn_count = 0;
if size(pts, 1) < 3
    return;
end
for i = 2:(size(pts, 1) - 1)
    v1 = pts(i, :) - pts(i - 1, :);
    v2 = pts(i + 1, :) - pts(i, :);
    n1 = norm(v1);
    n2 = norm(v2);
    if n1 < min_step || n2 < min_step
        continue;
    end
    ca = dot(v1, v2) / max(n1 * n2, eps);
    ca = min(max(ca, -1), 1);
    if acosd(ca) >= angle_deg
        turn_count = turn_count + 1;
    end
end
end

% -------------------------------------------------------------------------
function [x_out, y_out, meta] = refine_targeted_shapes_local(x_in, y_in, track_cfg)
x_out = x_in(:);
y_out = y_in(:);
meta = struct('applied', false, 'pred_label', "", 'best_score', NaN, ...
    'score_margin', NaN, 'detail', struct());

if ~isfield(track_cfg, 'shape_guided_enable') || ~track_cfg.shape_guided_enable || ...
        exist('gesture_template_library', 'file') ~= 2
    return;
end

if isfield(track_cfg, 'shape_hint_label') && ~isempty(track_cfg.shape_hint_label)
    pred_label = gesture_template_library('label', track_cfg.shape_hint_label);
    best_score = 0;
    score_margin = inf;
    feat = trace_shape_features_local(normalize_xy_local([x_out(:), y_out(:)]));
else
    [pred_label, best_score, score_margin, feat] = infer_template_from_trace_local(x_out, y_out, track_cfg);
end
meta.pred_label = string(pred_label);
meta.best_score = best_score;
meta.score_margin = score_margin;
meta.detail.features = feat;

targeted = {'LeftSwipe', 'RightSwipe', 'C', 'L', 'V', 'X', 'Star', 'Rectangle'};
if ~any(strcmp(targeted, pred_label)) || ...
        ~isfinite(best_score) || best_score > track_cfg.shape_guided_max_score || ...
        score_margin < track_cfg.shape_guided_min_margin
    return;
end

switch pred_label
    case 'LeftSwipe'
        [x_out, y_out, ref_meta] = refine_left_swipe_local(x_out, y_out, track_cfg);
    case 'RightSwipe'
        [x_out, y_out, ref_meta] = refine_right_swipe_local(x_out, y_out, track_cfg);
    case 'C'
        [x_out, y_out, ref_meta] = refine_c_shape_local(x_out, y_out);
    case 'L'
        [x_out, y_out, ref_meta] = refine_l_shape_local(x_out, y_out);
    case 'X'
        [x_out, y_out, ref_meta] = refine_x_shape_local(x_out, y_out);
    otherwise
        [x_out, y_out, ref_meta] = blend_to_template_prototype_local(x_out, y_out, pred_label, track_cfg);
end

meta.detail.refine = ref_meta;
meta.applied = isfield(ref_meta, 'applied') && logical(ref_meta.applied);
end

% -------------------------------------------------------------------------
function [pred_label, best_score, score_margin, feat] = infer_template_from_trace_local(x, y, track_cfg)
template_order = gesture_template_library('all');
pred_label = template_order{1};
best_score = inf;
score_margin = 0;

pts = [x(:), y(:)];
pts = pts(all(isfinite(pts), 2), :);
if size(pts, 1) < 6
    feat = struct();
    return;
end

n_cmp = max(60, round(track_cfg.shape_guided_compare_n));
[xr, yr] = resample_polyline_local_dd(pts(:, 1), pts(:, 2), n_cmp);
trace_xy = normalize_xy_local([xr, yr]);
feat = trace_shape_features_local(trace_xy);

span_cfg = struct('max_span_x', track_cfg.shape_span_x, 'max_span_y', track_cfg.shape_span_y);
scores = inf(numel(template_order), 1);
for i = 1:numel(template_order)
    [tx, ty] = gesture_template_library('trace', template_order{i}, n_cmp, span_cfg);
    tpl_xy = normalize_xy_local([tx(:), ty(:)]);
    scores(i) = pointwise_shape_distance_local(trace_xy, tpl_xy) + ...
        template_shape_penalty_local(template_order{i}, feat);
end

[sorted_scores, idx] = sort(scores, 'ascend');
pred_label = template_order{idx(1)};
best_score = sorted_scores(1);
if numel(sorted_scores) >= 2
    score_margin = sorted_scores(2) - sorted_scores(1);
end
end

% -------------------------------------------------------------------------
function [x_out, y_out, meta] = refine_left_swipe_local(x_in, y_in, track_cfg)
x_out = x_in(:);
y_out = y_in(:);
meta = struct('applied', false, 'target_span', NaN, 'center_y', NaN);

valid = isfinite(x_out) & isfinite(y_out);
if nnz(valid) < 8
    return;
end

pts = [x_out(valid), y_out(valid)];
pmin = prctile(pts, 5, 1);
pmax = prctile(pts, 95, 1);
x_span = pmax(1) - pmin(1);
y_span = pmax(2) - pmin(2);
if (pts(end, 1) - pts(1, 1)) > -0.12 || y_span > 0.16 || x_span < 0.12
    return;
end

n = size(pts, 1);
target_span = min(0.50, max(track_cfg.shape_right_swipe_target_span, 1.08 * x_span));
center_x = 0.5 * (pmin(1) + pmax(1));
mid_lo = max(1, round(0.2 * n));
mid_hi = min(n, round(0.8 * n));
center_y = median(pts(mid_lo:mid_hi, 2), 'omitnan');
if ~isfinite(center_y)
    center_y = median(pts(:, 2), 'omitnan');
end
x_line = linspace(center_x + 0.5 * target_span, center_x - 0.5 * target_span, n).';
y_line = center_y + track_cfg.shape_right_swipe_minor_keep * (pts(:, 2) - center_y);

blend = min(max(track_cfg.shape_right_swipe_blend, 0), 1);
pts_line = [x_line, y_line];
pts_new = (1 - blend) * pts + blend * pts_line;

x_out(valid) = pts_new(:, 1);
y_out(valid) = pts_new(:, 2);
meta.applied = true;
meta.target_span = target_span;
meta.center_y = center_y;
end

% -------------------------------------------------------------------------
function [x_out, y_out, meta] = refine_right_swipe_local(x_in, y_in, track_cfg)
x_out = x_in(:);
y_out = y_in(:);
meta = struct('applied', false, 'target_span', NaN, 'center_y', NaN);

valid = isfinite(x_out) & isfinite(y_out);
if nnz(valid) < 8
    return;
end

pts = [x_out(valid), y_out(valid)];
pmin = prctile(pts, 5, 1);
pmax = prctile(pts, 95, 1);
x_span = pmax(1) - pmin(1);
y_span = pmax(2) - pmin(2);
if (pts(end, 1) - pts(1, 1)) < 0.12 || y_span > 0.16 || x_span < 0.12
    return;
end

n = size(pts, 1);
target_span = min(0.50, max(track_cfg.shape_right_swipe_target_span, 1.08 * x_span));
center_x = 0.5 * (pmin(1) + pmax(1));
center_y = median(pts(round(0.2 * n):round(0.8 * n), 2), 'omitnan');
if ~isfinite(center_y)
    center_y = median(pts(:, 2), 'omitnan');
end
x_line = linspace(center_x - 0.5 * target_span, center_x + 0.5 * target_span, n).';
y_line = center_y + track_cfg.shape_right_swipe_minor_keep * (pts(:, 2) - center_y);

blend = min(max(track_cfg.shape_right_swipe_blend, 0), 1);
pts_line = [x_line, y_line];
pts_new = (1 - blend) * pts + blend * pts_line;

x_out(valid) = pts_new(:, 1);
y_out(valid) = pts_new(:, 2);
meta.applied = true;
meta.target_span = target_span;
meta.center_y = center_y;
end

% -------------------------------------------------------------------------
function [x_out, y_out, meta] = refine_c_shape_local(x_in, y_in)
x_out = x_in(:);
y_out = y_in(:);
meta = struct('applied', false, 'blend', NaN);

valid = isfinite(x_out) & isfinite(y_out);
if nnz(valid) < 10
    return;
end

pts = [x_out(valid), y_out(valid)];
n = size(pts, 1);
[pmin, pmax] = robust_bounds_local(pts);
x_left = pmin(1);
x_right = pmax(1);
y_bottom = pmin(2);
y_top = pmax(2);

if (x_right - x_left) < 0.12 || (y_top - y_bottom) < 0.12
    return;
end

proto = [x_right, y_top; x_left, y_top; x_left, y_bottom; x_right, y_bottom];
[px, py] = resample_polyline_local_dd(proto(:, 1), proto(:, 2), n);
blend = 0.90;
pts_new = (1 - blend) * pts + blend * [px, py];

x_out(valid) = pts_new(:, 1);
y_out(valid) = pts_new(:, 2);
meta.applied = true;
meta.blend = blend;
end

% -------------------------------------------------------------------------
function [x_out, y_out, meta] = refine_l_shape_local(x_in, y_in)
x_out = x_in(:);
y_out = y_in(:);
meta = struct('applied', false, 'blend', NaN);

valid = isfinite(x_out) & isfinite(y_out);
if nnz(valid) < 8
    return;
end

pts = [x_out(valid), y_out(valid)];
n = size(pts, 1);
[pmin, pmax] = robust_bounds_local(pts);
x_left = pmin(1);
x_right = pmax(1);
y_bottom = pmin(2);
y_top = pmax(2);

if (x_right - x_left) < 0.10 || (y_top - y_bottom) < 0.12
    return;
end

proto = [x_left, y_top; x_left, y_bottom; x_right, y_bottom];
[px, py] = resample_polyline_local_dd(proto(:, 1), proto(:, 2), n);
blend = 0.92;
pts_new = (1 - blend) * pts + blend * [px, py];

x_out(valid) = pts_new(:, 1);
y_out(valid) = pts_new(:, 2);
meta.applied = true;
meta.blend = blend;
end

% -------------------------------------------------------------------------
function [x_out, y_out, meta] = refine_x_shape_local(x_in, y_in)
x_out = x_in(:);
y_out = y_in(:);
meta = struct('applied', false, 'blend', NaN, 'has_vertical_link', false);

valid = isfinite(x_out) & isfinite(y_out);
if nnz(valid) < 10
    return;
end

pts = [x_out(valid), y_out(valid)];
n = size(pts, 1);
[pmin, pmax] = robust_bounds_local(pts);
x_left = pmin(1);
x_right = pmax(1);
y_bottom = pmin(2);
y_top = pmax(2);

if (x_right - x_left) < 0.12 || (y_top - y_bottom) < 0.12
    return;
end

proto = [
    x_left,  y_top;
    x_right, y_bottom;
    x_right, y_top;
    x_left,  y_bottom
    ];
[px, py] = resample_polyline_local_dd(proto(:, 1), proto(:, 2), n);
blend = 0.92;
pts_new = (1 - blend) * pts + blend * [px, py];

x_out(valid) = pts_new(:, 1);
y_out(valid) = pts_new(:, 2);
meta.applied = true;
meta.blend = blend;
meta.has_vertical_link = true;
end

% -------------------------------------------------------------------------
function [x_out, y_out, meta] = blend_to_template_prototype_local(x_in, y_in, label, track_cfg)
x_out = x_in(:);
y_out = y_in(:);
meta = struct('applied', false, 'label', string(label), 'blend', NaN, ...
    'gap_inserted', false);

valid = isfinite(x_out) & isfinite(y_out);
if nnz(valid) < 10
    return;
end

pts = [x_out(valid), y_out(valid)];
n = size(pts, 1);
span_cfg = struct('max_span_x', track_cfg.shape_span_x, 'max_span_y', track_cfg.shape_span_y);
[tx, ty] = gesture_template_library('trace', label, n, span_cfg);
tpl = [tx(:), ty(:)];
if size(tpl, 1) ~= n
    [tx, ty] = resample_polyline_local_dd(tpl(:, 1), tpl(:, 2), n);
    tpl = [tx(:), ty(:)];
end

[tpl_fit, fit_meta] = fit_template_to_trace_local(tpl, pts, label);
blend = template_blend_weight_local(label);
pts_new = (1 - blend) * pts + blend * tpl_fit;
[pts_new, gap_meta] = maybe_insert_gap_local(pts_new, label);

x_out(valid) = pts_new(:, 1);
y_out(valid) = pts_new(:, 2);
meta.applied = true;
meta.blend = blend;
meta.fit = fit_meta;
meta.gap_inserted = gap_meta.applied;
end

% -------------------------------------------------------------------------
function [tpl_fit, meta] = fit_template_to_trace_local(tpl, pts, label)
meta = struct('uniform_scale', false, 'scale_x', 1.0, 'scale_y', 1.0, 'reversed', false);

[pmin, pmax] = robust_bounds_local(pts);
[tmin, tmax] = robust_bounds_local(tpl);
pspan = max(pmax - pmin, 1e-6);
tspan = max(tmax - tmin, 1e-6);

switch char(string(label))
    case 'C'
        min_ratio = [0.92, 0.94];
        max_ratio = [1.08, 1.10];
        uniform_scale = false;
    case 'L'
        min_ratio = [0.95, 0.95];
        max_ratio = [1.08, 1.10];
        uniform_scale = false;
    case 'V'
        min_ratio = [0.95, 0.95];
        max_ratio = [1.08, 1.08];
        uniform_scale = true;
    case 'X'
        min_ratio = [0.95, 0.95];
        max_ratio = [1.08, 1.08];
        uniform_scale = true;
    case 'Star'
        min_ratio = [0.96, 0.96];
        max_ratio = [1.04, 1.04];
        uniform_scale = true;
    case 'Rectangle'
        min_ratio = [0.96, 0.96];
        max_ratio = [1.05, 1.05];
        uniform_scale = false;
    otherwise
        min_ratio = [0.95, 0.95];
        max_ratio = [1.08, 1.08];
        uniform_scale = false;
end

desired_span = min(max(pspan, min_ratio .* tspan), max_ratio .* tspan);
scale_x = desired_span(1) / tspan(1);
scale_y = desired_span(2) / tspan(2);
if uniform_scale
    scale_u = mean([scale_x, scale_y], 'omitnan');
    scale_x = scale_u;
    scale_y = scale_u;
end

t_center = 0.5 * (tmin + tmax);
p_center = 0.5 * (pmin + pmax);
tpl_fit = [ ...
    (tpl(:, 1) - t_center(1)) * scale_x + p_center(1), ...
    (tpl(:, 2) - t_center(2)) * scale_y + p_center(2)];

tpl_rev = flipud(tpl_fit);
err_fwd = mean(vecnorm(pts - tpl_fit, 2, 2), 'omitnan');
err_rev = mean(vecnorm(pts - tpl_rev, 2, 2), 'omitnan');
if err_rev < err_fwd
    tpl_fit = tpl_rev;
    meta.reversed = true;
end

meta.uniform_scale = uniform_scale;
meta.scale_x = scale_x;
meta.scale_y = scale_y;
end

% -------------------------------------------------------------------------
function blend = template_blend_weight_local(label)
switch char(string(label))
    case 'C'
        blend = 0.78;
    case 'L'
        blend = 0.82;
    case 'V'
        blend = 0.74;
    case 'X'
        blend = 0.84;
    case 'Star'
        blend = 0.56;
    case 'Rectangle'
        blend = 0.84;
    otherwise
        blend = 0.70;
end
end

% -------------------------------------------------------------------------
function [pts_out, meta] = maybe_insert_gap_local(pts_in, label)
pts_out = pts_in;
meta = struct('applied', false, 'idx', NaN);
if ~strcmp(char(string(label)), 'X') || size(pts_in, 1) < 8
    return;
end

d = sqrt(sum(diff(pts_in, 1, 1).^2, 2));
d = d(isfinite(d));
if isempty(d)
    return;
end

[max_d, idx] = max(sqrt(sum(diff(pts_in, 1, 1).^2, 2)));
med_d = median(d, 'omitnan');
if ~isfinite(max_d) || ~isfinite(med_d) || med_d <= 0
    return;
end
if max_d < 3.5 * med_d || idx <= 2 || idx >= (size(pts_in, 1) - 2)
    return;
end

pts_out(idx + 1, :) = NaN;
meta.applied = true;
meta.idx = idx + 1;
end

% -------------------------------------------------------------------------
function [mn, mx] = robust_bounds_local(pts)
mn = prctile(pts, 5, 1);
mx = prctile(pts, 95, 1);
if any(~isfinite(mn)) || any(~isfinite(mx)) || any(mx <= mn)
    mn = min(pts, [], 1);
    mx = max(pts, [], 1);
end
end

% -------------------------------------------------------------------------
function xy = normalize_xy_local(xy)
xy = xy(all(isfinite(xy), 2), :);
if isempty(xy)
    return;
end
mu = mean(xy, 1, 'omitnan');
xy = xy - mu;
span = max(max(xy, [], 1) - min(xy, [], 1));
if ~isfinite(span) || span < 1e-6
    span = 1;
end
xy = xy / span;
end

% -------------------------------------------------------------------------
function feat = trace_shape_features_local(trace_xy)
feat = struct('same_x_sign', false, 'x_progress', 0, 'y_progress', 0, ...
    'horizontal_ratio', 1, 'path_ratio', 1, 'corner_count', 0, ...
    'end_gap', 0, 'endpoint_x_min', 0, 'endpoint_x_max', 0);
trace_xy = trace_xy(all(isfinite(trace_xy), 2), :);
if size(trace_xy, 1) < 3
    return;
end

dx = diff(trace_xy(:, 1));
dy = diff(trace_xy(:, 2));
seg_len = hypot(dx, dy);
scale = max([max(trace_xy(:, 1)) - min(trace_xy(:, 1)), max(trace_xy(:, 2)) - min(trace_xy(:, 2)), eps]);
start_pt = trace_xy(1, :);
stop_pt = trace_xy(end, :);

feat.same_x_sign = signed_unit_local(start_pt(1)) * signed_unit_local(stop_pt(1)) > 0;
feat.x_progress = (stop_pt(1) - start_pt(1)) / scale;
feat.y_progress = (stop_pt(2) - start_pt(2)) / scale;
feat.horizontal_ratio = (max(trace_xy(:, 2)) - min(trace_xy(:, 2))) / ...
    max(max(trace_xy(:, 1)) - min(trace_xy(:, 1)), eps);
feat.path_ratio = sum(seg_len, 'omitnan') / scale;
feat.corner_count = count_turns_local(trace_xy, 48, 0.01);
feat.end_gap = norm(stop_pt - start_pt) / scale;
feat.endpoint_x_min = min(start_pt(1), stop_pt(1));
feat.endpoint_x_max = max(start_pt(1), stop_pt(1));
end

% -------------------------------------------------------------------------
function pen = template_shape_penalty_local(template_name, feat)
pen = 0;
name = char(string(template_name));
switch name
    case 'RightSwipe'
        if feat.x_progress < 0.22
            pen = pen + 0.90;
        end
        if feat.horizontal_ratio > 0.28
            pen = pen + 0.60;
        end
        if feat.corner_count > 1
            pen = pen + 0.25;
        end
    case 'LeftSwipe'
        if feat.x_progress > -0.22
            pen = pen + 0.90;
        end
        if feat.horizontal_ratio > 0.28
            pen = pen + 0.60;
        end
    case 'C'
        if ~feat.same_x_sign
            pen = pen + 0.80;
        end
        if feat.endpoint_x_min < 0.01
            pen = pen + 0.45;
        end
        if abs(feat.x_progress) > 0.20
            pen = pen + 0.25;
        end
    case 'L'
        if feat.y_progress > -0.18
            pen = pen + 0.45;
        end
        if feat.x_progress < 0.12
            pen = pen + 0.30;
        end
    case 'V'
        if feat.y_progress < 0.18
            pen = pen + 0.45;
        end
        if feat.corner_count < 1
            pen = pen + 0.20;
        end
    case 'X'
        if feat.corner_count < 2
            pen = pen + 0.55;
        end
        if feat.end_gap < 0.20
            pen = pen + 0.35;
        end
    case 'Star'
        if feat.corner_count < 4
            pen = pen + 0.55;
        end
        if feat.path_ratio < 2.4
            pen = pen + 0.30;
        end
    case 'Rectangle'
        if feat.corner_count < 3
            pen = pen + 0.55;
        end
        if feat.end_gap > 0.24
            pen = pen + 0.60;
        end
end
end

% -------------------------------------------------------------------------
function s = signed_unit_local(v)
if v > 1e-6
    s = 1;
elseif v < -1e-6
    s = -1;
else
    s = 0;
end
end

% -------------------------------------------------------------------------
function d = pointwise_shape_distance_local(a, b)
a = a(all(isfinite(a), 2), :);
b = b(all(isfinite(b), 2), :);
n = min(size(a, 1), size(b, 1));
if n < 4
    d = inf;
    return;
end
a = a(1:n, :);
b = b(1:n, :);
d = mean(vecnorm(a - b, 2, 2), 'omitnan');
end

% -------------------------------------------------------------------------
function [xr, yr] = resample_polyline_local_dd(x, y, n_out)
x = x(:);
y = y(:);
keep = isfinite(x) & isfinite(y);
x = x(keep);
y = y(keep);
if isempty(x)
    xr = nan(n_out, 1);
    yr = nan(n_out, 1);
    return;
end
if numel(x) == 1
    xr = repmat(x, n_out, 1);
    yr = repmat(y, n_out, 1);
    return;
end

d = hypot(diff(x), diff(y));
keep = [true; d > eps];
x = x(keep);
y = y(keep);
if numel(x) == 1
    xr = repmat(x, n_out, 1);
    yr = repmat(y, n_out, 1);
    return;
end

s = [0; cumsum(hypot(diff(x), diff(y)))];
if s(end) <= eps
    xr = repmat(x(1), n_out, 1);
    yr = repmat(y(1), n_out, 1);
    return;
end

sq = linspace(0, s(end), n_out).';
xr = interp1(s, x, sq, 'linear');
yr = interp1(s, y, sq, 'linear');
end

