function [traj_x, traj_y, traj_t, traj_conf, debug_info] = run_gesture_analysis_template_inverse(obs_in, nav_data, step1_res, user_cfg)
% RUN_GESTURE_ANALYSIS_TEMPLATE_INVERSE
% Global template-based inversion tailored for synthetic ideal data.
%
% This method uses the same GNSS occlusion geometry as simulation and
% solves a global objective over:
%   1) letter template identity,
%   2) temporal alignment (shift + scale),
%   3) spatial transform (anisotropic scale + translation).
%
% It is intentionally optimized for synthetic/ideal data first.

addpath(genpath('calculate_clock_bias_and_positon'));
addpath(genpath('nav_parse'));

cfg = default_cfg();
if nargin >= 4 && isstruct(user_cfg)
    cfg = merge_cfg(cfg, user_cfg);
end

required_fields = {'volatility_matrix', 't_grid', 'valid_sats', 'segments'};
for i = 1:numel(required_fields)
    if ~isfield(step1_res, required_fields{i})
        error('step1_res missing required field: %s', required_fields{i});
    end
end
if isempty(obs_in)
    [traj_x, traj_y, traj_t, traj_conf, debug_info] = empty_outputs();
    return;
end

vol_mat = step1_res.volatility_matrix;
t_grid = step1_res.t_grid(:);
valid_sats = step1_res.valid_sats(:);
segments = step1_res.segments;

[num_samples, num_sats] = size(vol_mat);
if num_samples ~= numel(t_grid)
    error('Size mismatch: volatility_matrix rows and t_grid length are different.');
end
if num_sats ~= numel(valid_sats)
    error('Size mismatch: volatility_matrix columns and valid_sats length are different.');
end

if cfg.debug.verbose
    fprintf('--> [TemplateInverse] Start global inversion...\n');
    fprintf('    Samples=%d, Sats=%d\n', num_samples, num_sats);
end

% Build active mask from Step1 segments.
active_mask = false(num_samples, 1);
if isempty(segments)
    active_mask(:) = true;
else
    pad = max(0, round(cfg.track.segment_pad_frames));
    for i = 1:numel(segments)
        s = max(1, segments(i).start_idx - pad);
        e = min(num_samples, segments(i).end_idx + pad);
        if e >= s
            active_mask(s:e) = true;
        end
    end
end

frame_energy = sum(max(vol_mat, 0), 2, 'omitnan');
process_mask = active_mask & (frame_energy > cfg.track.min_frame_energy);
if nnz(process_mask) < cfg.track.min_process_frames
    process_mask = active_mask & (frame_energy > 0);
end
if nnz(process_mask) < cfg.track.min_process_frames
    process_mask = active_mask;
end

proc_indices = find(process_mask);
if isempty(proc_indices)
    [traj_x, traj_y, traj_t, traj_conf, debug_info] = empty_outputs();
    return;
end

if cfg.track.frame_stride > 1
    proc_indices = proc_indices(1:cfg.track.frame_stride:end);
end
if numel(proc_indices) > cfg.track.max_frames
    proc_indices = proc_indices(1:cfg.track.max_frames);
end
num_proc = numel(proc_indices);

if cfg.debug.verbose
    fprintf('    ActiveFrames=%d, ProcessFrames=%d\n', nnz(active_mask), num_proc);
end

% Epoch mapping.
obs_times = [obs_in.time];
obs_posix = posixtime(obs_times(:));
proc_posix = posixtime(t_grid(proc_indices));
epoch_idx_proc = round(interp1(obs_posix, 1:numel(obs_posix), proc_posix, 'nearest', 'extrap'));
epoch_idx_proc = min(max(epoch_idx_proc, 1), numel(obs_times));

% Build projection geometry.
unique_epochs = unique(epoch_idx_proc);
geo_cache = repmat(struct('ok', false, 'proj_xy', [], 'elev_deg', []), numel(obs_times), 1);
geo_ref = build_geometry_reference(obs_in, nav_data, unique_epochs, cfg);
for i = 1:numel(unique_epochs)
    ep = unique_epochs(i);
    [geo_cache(ep).ok, geo_cache(ep).proj_xy, geo_cache(ep).elev_deg] = ...
        build_epoch_geometry(obs_in, nav_data, ep, valid_sats, cfg, geo_ref);
end

[q_obs, sat_w, geom_x, geom_y, geom_ok] = build_frame_observation_mats( ...
    vol_mat, proc_indices, epoch_idx_proc, geo_cache, cfg);

valid_frames = find(sum(geom_ok, 2) >= cfg.obs.min_sats);
if numel(valid_frames) < cfg.track.min_process_frames
    if cfg.debug.verbose
        fprintf('    Not enough valid frames for inversion.\n');
    end
    [traj_x, traj_y, traj_t, traj_conf, debug_info] = empty_outputs();
    return;
end

letters = cfg.template.letters;
best = struct('score', inf, 'letter', '', 'params', [], 'time_score', inf);
candidate_log = repmat(struct('letter', '', 'best_score', inf, 'best_params', []), numel(letters), 1);

energy_proc = frame_energy(proc_indices);
energy_proc = normalize_vec(energy_proc);
energy_proc = smoothdata(energy_proc, 'movmean', cfg.search.energy_smooth_win);

for li = 1:numel(letters)
    letter = upper(strtrim(char(letters{li})));
    [time_cands, time_cand_score] = select_time_candidates(letter, energy_proc, cfg);
    local_best_score = inf;
    local_best_params = [];

    for ti = 1:size(time_cands, 1)
        shift0 = time_cands(ti, 1);
        tscale0 = time_cands(ti, 2);
        [bx, by, pen] = synth_template(letter, num_proc, shift0, tscale0);

        for sx = cfg.search.sx_grid
            for sy = cfg.search.sy_grid
                x0 = sx * bx;
                y0 = sy * by;
                for tx = cfg.search.tx_grid
                    for ty = cfg.search.ty_grid
                        xh = x0 + tx;
                        yh = y0 + ty;
                        [~, score] = forward_loss(xh, yh, pen, q_obs, sat_w, geom_x, geom_y, geom_ok, cfg);
                        if score < local_best_score
                            local_best_score = score;
                            local_best_params = [sx, sy, tx, ty, shift0, tscale0];
                        end
                    end
                end
            end
        end
    end

    if isempty(local_best_params)
        continue;
    end

    [local_best_params, local_best_score] = refine_params_coordinate( ...
        local_best_params, letter, q_obs, sat_w, geom_x, geom_y, geom_ok, cfg);

    candidate_log(li).letter = letter;
    candidate_log(li).best_score = local_best_score;
    candidate_log(li).best_params = local_best_params;

    if local_best_score < best.score
        best.score = local_best_score;
        best.letter = letter;
        best.params = local_best_params;
        best.time_score = time_cand_score;
    end
end

if isempty(best.params)
    [traj_x, traj_y, traj_t, traj_conf, debug_info] = empty_outputs();
    return;
end

[bx, by, pen] = synth_template(best.letter, num_proc, best.params(5), best.params(6));
traj_proc_x = best.params(1) * bx + best.params(3);
traj_proc_y = best.params(2) * by + best.params(4);

[frame_loss, ~] = forward_loss(traj_proc_x, traj_proc_y, pen, q_obs, sat_w, geom_x, geom_y, geom_ok, cfg);
loss_scale = max(cfg.conf.loss_scale, 1e-6);
conf_proc = exp(-frame_loss / loss_scale);
conf_proc = min(max(conf_proc, 0), 1);

traj_full_x = nan(num_samples, 1);
traj_full_y = nan(num_samples, 1);
conf_full = nan(num_samples, 1);
traj_full_x(proc_indices) = traj_proc_x;
traj_full_y(proc_indices) = traj_proc_y;
conf_full(proc_indices) = conf_proc;

output_mask = active_mask;
out_idx = find(output_mask);

fit_idx = find(isfinite(traj_full_x) & isfinite(traj_full_y));
if numel(fit_idx) >= 2
    mid = out_idx(out_idx >= fit_idx(1) & out_idx <= fit_idx(end));
    head = out_idx(out_idx < fit_idx(1));
    tail = out_idx(out_idx > fit_idx(end));

    if ~isempty(mid)
        traj_full_x(mid) = interp1(fit_idx, traj_full_x(fit_idx), mid, 'linear');
        traj_full_y(mid) = interp1(fit_idx, traj_full_y(fit_idx), mid, 'linear');
    end
    if ~isempty(head)
        traj_full_x(head) = traj_full_x(fit_idx(1));
        traj_full_y(head) = traj_full_y(fit_idx(1));
    end
    if ~isempty(tail)
        traj_full_x(tail) = traj_full_x(fit_idx(end));
        traj_full_y(tail) = traj_full_y(fit_idx(end));
    end

    cf_idx = find(isfinite(conf_full));
    if numel(cf_idx) >= 2
        conf_full(mid) = interp1(cf_idx, conf_full(cf_idx), mid, 'linear');
    elseif numel(cf_idx) == 1
        conf_full(mid) = conf_full(cf_idx);
    end
    if ~isempty(head) && isfinite(conf_full(fit_idx(1)))
        conf_full(head) = conf_full(fit_idx(1));
    end
    if ~isempty(tail) && isfinite(conf_full(fit_idx(end)))
        conf_full(tail) = conf_full(fit_idx(end));
    end
elseif numel(fit_idx) == 1
    traj_full_x(out_idx) = traj_full_x(fit_idx);
    traj_full_y(out_idx) = traj_full_y(fit_idx);
    conf_full(out_idx) = conf_full(fit_idx);
end

if cfg.track.final_smooth_pts > 1 && ~isempty(out_idx)
    traj_full_x(out_idx) = smoothdata(traj_full_x(out_idx), 'movmean', cfg.track.final_smooth_pts);
    traj_full_y(out_idx) = smoothdata(traj_full_y(out_idx), 'movmean', cfg.track.final_smooth_pts);
end

traj_t = out_idx(isfinite(traj_full_x(out_idx)) & isfinite(traj_full_y(out_idx)));
traj_x = traj_full_x(traj_t);
traj_y = traj_full_y(traj_t);
traj_conf = conf_full(traj_t);
traj_conf = min(max(traj_conf, 0), 1);

debug_info = struct();
debug_info.cfg = cfg;
debug_info.best_letter = best.letter;
debug_info.best_score = best.score;
debug_info.best_params = best.params;
debug_info.candidates = candidate_log;
debug_info.process_mask = process_mask;
debug_info.active_mask = active_mask;
debug_info.proc_indices = proc_indices;
debug_info.epoch_idx_proc = epoch_idx_proc;
debug_info.frame_energy_proc = energy_proc;
debug_info.traj_full_x = traj_full_x;
debug_info.traj_full_y = traj_full_y;
debug_info.conf_full = conf_full;
debug_info.output_mask = output_mask;
debug_info.output_idx = traj_t;

if cfg.debug.plot
    plot_debug_template(t_grid, active_mask, traj_full_x, traj_full_y, conf_full, best.letter, cfg);
end

if cfg.debug.verbose
    fprintf('    BestLetter=%s, Score=%.5f, Output=%d\n', best.letter, best.score, numel(traj_t));
end
end

% -------------------------------------------------------------------------
function [q_obs, sat_w, geom_x, geom_y, geom_ok] = build_frame_observation_mats( ...
    vol_mat, proc_indices, epoch_idx_proc, geo_cache, cfg)

num_proc = numel(proc_indices);
num_sats = size(vol_mat, 2);

q_obs = nan(num_proc, num_sats);
geom_x = nan(num_proc, num_sats);
geom_y = nan(num_proc, num_sats);
geom_el = nan(num_proc, num_sats);
geom_ok = false(num_proc, num_sats);

for n = 1:num_proc
    t_idx = proc_indices(n);
    ep = epoch_idx_proc(n);

    q_obs(n, :) = vol_mat(t_idx, :);
    if ep < 1 || ep > numel(geo_cache) || ~geo_cache(ep).ok
        continue;
    end

    gx = geo_cache(ep).proj_xy(:, 1)';
    gy = geo_cache(ep).proj_xy(:, 2)';
    ge = geo_cache(ep).elev_deg(:)';

    geom_x(n, :) = gx;
    geom_y(n, :) = gy;
    geom_el(n, :) = ge;
    geom_ok(n, :) = isfinite(gx) & isfinite(gy) & isfinite(ge) & (ge >= cfg.obs.min_elevation_deg);
end

% Sat-wise normalization into [0,1].
for s = 1:num_sats
    v = q_obs(:, s);
    good = isfinite(v);
    if nnz(good) < 4
        q_obs(:, s) = 0;
        continue;
    end

    ql = quick_quantile(v(good), cfg.obs.norm_q_low);
    qh = quick_quantile(v(good), cfg.obs.norm_q_high);
    span = max(qh - ql, 1e-8);
    q_obs(:, s) = (v - ql) / span;
end
q_obs(~isfinite(q_obs)) = 0;
q_obs = min(max(q_obs, 0), 1);

% Frame/sat weights.
sat_w = zeros(num_proc, num_sats);
for n = 1:num_proc
    idx = find(geom_ok(n, :));
    if isempty(idx)
        continue;
    end

    el = geom_el(n, idx);
    we = (sind(max(el, 0))).^cfg.obs.elevation_power;
    wq = 0.2 + 0.8 * q_obs(n, idx);
    w = we .* wq;
    w = w / (sum(w) + eps);
    sat_w(n, idx) = w;
end
end

% -------------------------------------------------------------------------
function [time_cands, best_time_score] = select_time_candidates(letter, energy_proc, cfg)
sh = cfg.search.time_shift_grid(:);
sc = cfg.search.time_scale_grid(:);

num_proc = numel(energy_proc);
score_list = nan(numel(sh) * numel(sc), 3);
kk = 0;

for i = 1:numel(sh)
    for j = 1:numel(sc)
        kk = kk + 1;
        [~, ~, pen] = synth_template(letter, num_proc, sh(i), sc(j));
        pen = double(pen(:));
        pen = smoothdata(pen, 'movmean', cfg.search.pen_smooth_win);
        pen = normalize_vec(pen);

        corr_val = corr(energy_proc(:), pen(:), 'Rows', 'pairwise');
        if ~isfinite(corr_val)
            corr_val = 0;
        end
        mse_val = mean((energy_proc(:) - pen(:)).^2, 'omitnan');
        score = -corr_val + cfg.search.time_mse_weight * mse_val;
        score_list(kk, :) = [score, sh(i), sc(j)];
    end
end

[~, ord] = sort(score_list(:, 1), 'ascend');
take = min(cfg.search.top_time_candidates, numel(ord));
time_cands = score_list(ord(1:take), 2:3);
best_time_score = score_list(ord(1), 1);
end

% -------------------------------------------------------------------------
function [params, best_score] = refine_params_coordinate(params0, letter, q_obs, sat_w, geom_x, geom_y, geom_ok, cfg)
params = params0(:)';
[bx, by, pen] = synth_template(letter, size(q_obs, 1), params(5), params(6));
[~, best_score] = forward_loss(params(1) * bx + params(3), params(2) * by + params(4), ...
    pen, q_obs, sat_w, geom_x, geom_y, geom_ok, cfg);

step_sets = cfg.search.refine_steps;
for si = 1:numel(step_sets)
    step = step_sets{si};
    improved = true;
    while improved
        improved = false;
        for d = 1:6
            cand1 = params;
            cand2 = params;
            cand1(d) = cand1(d) - step(d);
            cand2(d) = cand2(d) + step(d);
            cand1 = clamp_params(cand1, cfg);
            cand2 = clamp_params(cand2, cfg);

            [s1] = eval_param_score(cand1, letter, q_obs, sat_w, geom_x, geom_y, geom_ok, cfg);
            [s2] = eval_param_score(cand2, letter, q_obs, sat_w, geom_x, geom_y, geom_ok, cfg);

            if s1 < best_score || s2 < best_score
                if s1 <= s2
                    params = cand1;
                    best_score = s1;
                else
                    params = cand2;
                    best_score = s2;
                end
                improved = true;
            end
        end
    end
end
end

% -------------------------------------------------------------------------
function score = eval_param_score(params, letter, q_obs, sat_w, geom_x, geom_y, geom_ok, cfg)
[bx, by, pen] = synth_template(letter, size(q_obs, 1), params(5), params(6));
 [~, score] = forward_loss(params(1) * bx + params(3), params(2) * by + params(4), ...
    pen, q_obs, sat_w, geom_x, geom_y, geom_ok, cfg);
end

% -------------------------------------------------------------------------
function params = clamp_params(params, cfg)
params(1) = min(max(params(1), cfg.search.sx_range(1)), cfg.search.sx_range(2));
params(2) = min(max(params(2), cfg.search.sy_range(1)), cfg.search.sy_range(2));
params(3) = min(max(params(3), cfg.search.tx_range(1)), cfg.search.tx_range(2));
params(4) = min(max(params(4), cfg.search.ty_range(1)), cfg.search.ty_range(2));
params(5) = min(max(params(5), cfg.search.time_shift_range(1)), cfg.search.time_shift_range(2));
params(6) = min(max(params(6), cfg.search.time_scale_range(1)), cfg.search.time_scale_range(2));
end

% -------------------------------------------------------------------------
function [frame_loss, total_loss] = forward_loss(xh, yh, pen, q_obs, sat_w, geom_x, geom_y, geom_ok, cfg)
num_proc = numel(xh);
frame_loss = ones(num_proc, 1) * cfg.loss.missing_frame_penalty;
B = cfg.model.body_pos;
sigma2 = max((cfg.model.arm_width * cfg.model.sigma_ratio)^2, 1e-8);

for n = 1:num_proc
    idx = find(geom_ok(n, :));
    if numel(idx) < cfg.obs.min_sats
        continue;
    end

    if ~pen(n)
        p_occ = zeros(numel(idx), 1);
    else
        Ax = xh(n);
        Ay = yh(n);
        ABx = B(1) - Ax;
        ABy = B(2) - Ay;
        len2 = ABx * ABx + ABy * ABy + 1e-9;

        Px = geom_x(n, idx)';
        Py = geom_y(n, idx)';
        APx = Px - Ax;
        APy = Py - Ay;
        t = (APx * ABx + APy * ABy) / len2;
        t = min(max(t, 0), 1);
        Cx = Ax + t * ABx;
        Cy = Ay + t * ABy;
        d2 = (Px - Cx).^2 + (Py - Cy).^2;
        p_occ = exp(-0.5 * d2 / sigma2);
    end

    q = q_obs(n, idx)';
    w = sat_w(n, idx)';
    if ~any(isfinite(w)) || sum(w) <= 0
        w = ones(size(q)) / max(numel(q), 1);
    end
    w = w / (sum(w) + eps);

    err = q - p_occ;
    frame_loss(n) = sum(w .* (err.^2));
end

if cfg.loss.use_energy_weight
    fw = normalize_vec(sum(q_obs, 2));
    fw = cfg.loss.energy_floor + (1 - cfg.loss.energy_floor) * fw;
    total_loss = sum(fw .* frame_loss, 'omitnan') / max(sum(fw, 'omitnan'), eps);
else
    total_loss = mean(frame_loss, 'omitnan');
end
end

% -------------------------------------------------------------------------
function [x, y, pen] = synth_template(letter, n, shift, tscale)
p = linspace(0, 1, n)';
ph_raw = (p - shift) ./ max(tscale, 1e-6);
ph = min(max(ph_raw, 0), 1);
[x, y, pen] = sample_template_phase(letter, ph);
outside = (ph_raw < 0) | (ph_raw > 1);
pen(outside) = false;
end

% -------------------------------------------------------------------------
function [x, y, pen] = sample_template_phase(letter, phase)
phase = min(max(phase(:), 0), 1);
stages = letter_stages(letter);
n_stage = size(stages, 1);

dur = zeros(n_stage, 1);
cum = zeros(n_stage + 1, 1);
for i = 1:n_stage
    dur(i) = stages{i, 3};
    cum(i + 1) = cum(i) + dur(i);
end
tot = cum(end);
t = phase * tot;

x = nan(size(phase));
y = nan(size(phase));
pen = false(size(phase));

for k = 1:numel(phase)
    tk = t(k);
    sid = find(tk <= cum(2:end), 1, 'first');
    if isempty(sid)
        sid = n_stage;
    end

    t0 = cum(sid);
    t1 = cum(sid + 1);
    a = (tk - t0) / max(t1 - t0, eps);
    a = min(max(a, 0), 1);

    p1 = stages{sid, 1};
    p2 = stages{sid, 2};
    pp = p1 + a * (p2 - p1);
    x(k) = pp(1);
    y(k) = pp(2);
    pen(k) = logical(stages{sid, 4});
end
end

% -------------------------------------------------------------------------
function stages = letter_stages(letter)
switch upper(string(letter))
    case "A"
        P1 = [-0.40, -0.50]; P2 = [0.00, 0.50]; P3 = [0.40, -0.50];
        P4 = [-0.20, -0.10]; P5 = [0.20, -0.10];
        stages = {
            P1, P2, 1.5, true;
            P2, P2, 3.0, true;
            P2, P3, 1.5, true;
            P3, P4, 3.0, false;
            P4, P5, 1.5, true
        };
    case "B"
        P1 = [-0.40, -1.00]; P2 = [-0.40, 1.00]; P3 = [1.50, 0.40];
        P4 = [-0.40, 0.00]; P5 = [1.50, 0.00]; P6 = [-0.40, -1.00];
        stages = {
            P1, P2, 1.5, true;  P2, P2, 3.0, true;
            P2, P3, 1.5, true;  P3, P3, 3.0, true;
            P3, P4, 3.0, true;  P4, P4, 3.0, true;
            P4, P5, 1.5, true;  P5, P5, 3.0, true;
            P5, P6, 1.5, true
        };
    case "M"
        P1 = [-0.40, -0.40]; P2 = [-0.40, 0.40]; P3 = [0.00, 0.00];
        P4 = [0.40, 0.40]; P5 = [0.40, -0.40];
        stages = {
            P1, P2, 1.5, true;  P2, P2, 3.0, true;
            P2, P3, 1.5, true;  P3, P3, 3.0, true;
            P3, P4, 1.5, true;  P4, P4, 3.0, true;
            P4, P5, 1.5, true
        };
    case "STAR"
        P1 = [-0.30, -0.45]; P2 = [0.00, 0.55]; P3 = [0.30, -0.45];
        P4 = [-0.48, 0.15]; P5 = [0.48, 0.15];
        stages = {
            P1, P2, 1.5, true;  P2, P2, 3.0, true;
            P2, P3, 1.5, true;  P3, P3, 3.0, true;
            P3, P4, 1.5, true;  P4, P4, 3.0, true;
            P4, P5, 1.5, true;  P5, P5, 3.0, true;
            P5, P1, 1.5, true
        };
    case "L"
        P1 = [-0.3, 0.4]; P2 = [-0.3, -0.4]; P3 = [0.3, -0.4];
        stages = {
            P1, P2, 1.5, true;
            P2, P2, 3.0, true;
            P2, P3, 1.5, true
        };
    case "X"
        P1 = [-0.3, 0.4]; P2 = [0.3, -0.4]; P3 = [0.3, 0.4]; P4 = [-0.3, -0.4];
        stages = {
            P1, P2, 1.5, true;
            P2, P3, 3.0, false;
            P3, P4, 1.5, true
        };
    case "Z"
        P1 = [-0.3, 0.4]; P2 = [0.3, 0.4]; P3 = [-0.3, -0.4]; P4 = [0.7, -0.4];
        stages = {
            P1, P2, 1.5, true;
            P2, P2, 3.0, true;
            P2, P3, 1.5, true;
            P3, P3, 3.0, true;
            P3, P4, 1.5, true
        };
    case "N"
        P1 = [-0.3, -0.4]; P2 = [-0.3, 0.4]; P3 = [0.3, -0.4]; P4 = [0.3, 0.4];
        stages = {
            P1, P2, 1.5, true;
            P2, P2, 3.0, true;
            P2, P3, 1.5, true;
            P3, P3, 3.0, true;
            P3, P4, 1.5, true
        };
    otherwise
        error('Unsupported template letter: %s', char(letter));
end
end

% -------------------------------------------------------------------------
function geo_ref = build_geometry_reference(obs_in, nav_data, epoch_list, cfg)
geo_ref = struct( ...
    'ok', false, ...
    'rec_pos_mean', [NaN, NaN, NaN], ...
    'lat0', NaN, ...
    'lon0', NaN, ...
    'alt0', NaN, ...
    'used_epochs', 0);

if isempty(epoch_list)
    return;
end

epochs = unique(epoch_list(:));
if isempty(epochs)
    return;
end

target_n = max(1, cfg.model.ref_epoch_samples);
if numel(epochs) > target_n
    sel = round(linspace(1, numel(epochs), target_n));
    epochs = epochs(sel);
end

rec_buf = nan(numel(epochs), 3);
n_ok = 0;
for i = 1:numel(epochs)
    ep = epochs(i);
    try
        [rec_pos, ~, ~] = calculate_receiver_position(obs_in, nav_data, ep);
    catch
        continue;
    end
    if all(isfinite(rec_pos))
        n_ok = n_ok + 1;
        rec_buf(n_ok, :) = rec_pos(:)';
    end
end

if n_ok < cfg.model.min_ref_epochs
    return;
end

rec_mean = mean(rec_buf(1:n_ok, :), 1);
if any(~isfinite(rec_mean))
    return;
end

[lat0, lon0, alt0] = ecef2geodetic(rec_mean(1), rec_mean(2), rec_mean(3));
if any(~isfinite([lat0, lon0, alt0]))
    return;
end

geo_ref.ok = true;
geo_ref.rec_pos_mean = rec_mean;
geo_ref.lat0 = lat0;
geo_ref.lon0 = lon0;
geo_ref.alt0 = alt0;
geo_ref.used_epochs = n_ok;
end

% -------------------------------------------------------------------------
function [ok, proj_xy, elev_deg] = build_epoch_geometry(obs_in, nav_data, epoch_idx, valid_sats, cfg, geo_ref)
ok = false;
num_sats = numel(valid_sats);
proj_xy = nan(num_sats, 2);
elev_deg = nan(num_sats, 1);

try
    [rec_pos, ~, sat_states, lat0, lon0, alt0] = calculate_receiver_position(obs_in, nav_data, epoch_idx);
catch
    return;
end

if any(~isfinite(rec_pos))
    return;
end

if geo_ref.ok
    ref_pos = geo_ref.rec_pos_mean;
    ref_lat = geo_ref.lat0;
    ref_lon = geo_ref.lon0;
    ref_alt = geo_ref.alt0;
else
    ref_pos = rec_pos;
    ref_lat = lat0;
    ref_lon = lon0;
    ref_alt = alt0;
end

for s = 1:num_sats
    sid = valid_sats{s};
    if ~isfield(sat_states, sid)
        continue;
    end

    sat_p = sat_states.(sid).position;
    [e, n, u] = ecef2enu( ...
        sat_p(1) - ref_pos(1), ...
        sat_p(2) - ref_pos(2), ...
        sat_p(3) - ref_pos(3), ...
        ref_lat, ref_lon, ref_alt);

    dist = norm([e, n, u]);
    if dist <= 0
        continue;
    end

    el = asind(u / dist);
    elev_deg(s) = el;
    if u <= 0 || el < cfg.obs.min_elevation_deg
        continue;
    end

    sc = cfg.model.gesture_height / u;
    px = sc * e;
    py = sc * n;
    if hypot(px, py) > cfg.obs.max_projection_radius
        continue;
    end

    proj_xy(s, :) = [px, py];
end

ok = true;
end

% -------------------------------------------------------------------------
function v = normalize_vec(v)
v = v(:);
good = isfinite(v);
if ~any(good)
    v(:) = 0;
    return;
end
lo = min(v(good));
hi = max(v(good));
if hi > lo
    v = (v - lo) ./ (hi - lo);
else
    v(:) = 0;
end
v(~isfinite(v)) = 0;
v = min(max(v, 0), 1);
end

% -------------------------------------------------------------------------
function q = quick_quantile(x, p)
x = x(isfinite(x));
if isempty(x)
    q = 0;
    return;
end
x = sort(x);
n = numel(x);
pos = 1 + (n - 1) * p;
lo = floor(pos);
hi = ceil(pos);
if lo == hi
    q = x(lo);
else
    a = pos - lo;
    q = (1 - a) * x(lo) + a * x(hi);
end
end

% -------------------------------------------------------------------------
function cfg = default_cfg()
cfg = struct();

cfg.model = struct();
cfg.model.gesture_height = 0.30;
cfg.model.arm_width = 0.15;
cfg.model.sigma_ratio = 0.55;
cfg.model.body_pos = [0.0, -1.0];
cfg.model.ref_epoch_samples = 24;
cfg.model.min_ref_epochs = 4;

cfg.obs = struct();
cfg.obs.min_elevation_deg = 5;
cfg.obs.max_projection_radius = 2.5;
cfg.obs.min_sats = 4;
cfg.obs.elevation_power = 2;
cfg.obs.norm_q_low = 0.15;
cfg.obs.norm_q_high = 0.85;

cfg.track = struct();
cfg.track.min_frame_energy = 0.15;
cfg.track.min_process_frames = 24;
cfg.track.segment_pad_frames = 20;
cfg.track.frame_stride = 1;
cfg.track.max_frames = 3000;
cfg.track.final_smooth_pts = 5;

cfg.template = struct();
cfg.template.letters = {'A', 'B', 'M', 'Star', 'L', 'X', 'Z', 'N'};

cfg.search = struct();
cfg.search.energy_smooth_win = 11;
cfg.search.pen_smooth_win = 9;
cfg.search.time_mse_weight = 0.35;
cfg.search.top_time_candidates = 5;

cfg.search.time_shift_grid = -0.50:0.10:0.50;
cfg.search.time_scale_grid = 0.25:0.10:1.50;
cfg.search.sx_grid = [0.86, 1.00, 1.14];
cfg.search.sy_grid = [0.86, 1.00, 1.14];
cfg.search.tx_grid = [-0.12, 0.00, 0.12];
cfg.search.ty_grid = [-0.14, 0.00, 0.14];

cfg.search.sx_range = [0.70, 1.35];
cfg.search.sy_range = [0.70, 1.35];
cfg.search.tx_range = [-0.40, 0.40];
cfg.search.ty_range = [-0.40, 0.40];
cfg.search.time_shift_range = [-0.80, 0.80];
cfg.search.time_scale_range = [0.15, 2.20];
cfg.search.refine_steps = {
    [0.05, 0.05, 0.05, 0.05, 0.06, 0.10], ...
    [0.02, 0.02, 0.02, 0.02, 0.02, 0.04]
};

cfg.loss = struct();
cfg.loss.missing_frame_penalty = 0.8;
cfg.loss.use_energy_weight = true;
cfg.loss.energy_floor = 0.05;

cfg.conf = struct();
cfg.conf.loss_scale = 0.18;

cfg.debug = struct();
cfg.debug.verbose = true;
cfg.debug.plot = false;
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
function plot_debug_template(t_grid, active_mask, x_full, y_full, conf_full, best_letter, cfg)
idx = find(active_mask);
if isempty(idx)
    return;
end

figure('Name', sprintf('TemplateInverse Debug [%s]', best_letter), ...
    'Position', [100, 80, 1120, 520], 'Color', 'w');

subplot(1, 2, 1);
hold on; grid on; axis equal;
plot(0, 0, '^k', 'MarkerFaceColor', 'k', 'DisplayName', 'Receiver');
plot(cfg.model.body_pos(1), cfg.model.body_pos(2), 'bs', 'MarkerFaceColor', 'b', 'DisplayName', 'Body Ref');
plot(x_full(idx), y_full(idx), 'r-', 'LineWidth', 2, 'DisplayName', 'TemplateInverse');
plot(x_full(idx(1)), y_full(idx(1)), 'go', 'MarkerFaceColor', 'g', 'DisplayName', 'Start');
plot(x_full(idx(end)), y_full(idx(end)), 'mo', 'MarkerFaceColor', 'm', 'DisplayName', 'End');
xlabel('East (m)'); ylabel('North (m)');
title(sprintf('Trajectory (%s)', best_letter));
legend('Location', 'best');

subplot(2, 2, 2);
plot(t_grid(idx), conf_full(idx), 'b-', 'LineWidth', 1.5);
grid on; ylim([0, 1.05]);
title('Confidence');
datetick('x', 'HH:MM:ss', 'keepticks', 'keeplimits');

subplot(2, 2, 4);
plot(t_grid(idx), x_full(idx), 'k-', 'LineWidth', 1.4); hold on;
plot(t_grid(idx), y_full(idx), 'r-', 'LineWidth', 1.2);
grid on;
legend({'X', 'Y'}, 'Location', 'best');
title('Coordinates');
datetick('x', 'HH:MM:ss', 'keepticks', 'keeplimits');
end

% -------------------------------------------------------------------------
function [traj_x, traj_y, traj_t, traj_conf, debug_info] = empty_outputs()
traj_x = [];
traj_y = [];
traj_t = [];
traj_conf = [];
debug_info = struct();
end
