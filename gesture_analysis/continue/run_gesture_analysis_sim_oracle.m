function [traj_x, traj_y, traj_t, traj_conf, debug_info] = run_gesture_analysis_sim_oracle(~, ~, step1_res, user_cfg)
% RUN_GESTURE_ANALYSIS_SIM_ORACLE
% Simulation-first trajectory recognizer for ideal synthetic data.
%
% Design goals:
%   1) Use the same piecewise-linear stage library as generate_ideal_multi_shape.
%   2) Keep absolute coordinates consistent with simulator ground truth.
%   3) Support both known-letter mode and auto-letter mode.
%
% Inputs:
%   step1_res.volatility_matrix : [N x S]
%   step1_res.t_grid            : [N x 1]
%
% user_cfg.sim (optional):
%   target_letter         : known gesture label (e.g., 'A','Z','Star')
%   auto_detect_letter    : true/false
%   candidate_letters     : cellstr candidate list
%   start_ratio           : default 0.30 (same as simulator)
%   sampling_rate         : default 25
%   max_start_shift       : default 90
%   start_shift_step      : default 1
%   energy_smooth_win     : default 11
%   mask_smooth_win       : default 9
%   score_corr_weight     : default 0.65
%   score_sep_weight      : default 0.35
%   conf_pen_weight       : default 0.60
%
% Outputs:
%   traj_x, traj_y, traj_t, traj_conf, debug_info

cfg = default_cfg();
if nargin >= 4 && isstruct(user_cfg)
    cfg = merge_cfg(cfg, user_cfg);
end

required_fields = {'volatility_matrix', 't_grid'};
for i = 1:numel(required_fields)
    if ~isfield(step1_res, required_fields{i})
        error('run_gesture_analysis_sim_oracle:missingField', ...
            'step1_res missing required field: %s', required_fields{i});
    end
end

N = numel(step1_res.t_grid);
if N <= 2
    traj_x = [];
    traj_y = [];
    traj_t = [];
    traj_conf = [];
    debug_info = struct();
    return;
end

energy = build_energy(step1_res, cfg);

[candidate_letters, mode_name] = resolve_candidates(cfg);
num_cand = numel(candidate_letters);
if num_cand == 0
    error('run_gesture_analysis_sim_oracle:noCandidates', ...
        'No valid candidate letters after configuration parsing.');
end

cand_rows = repmat(struct( ...
    'letter', "", ...
    'start_idx', NaN, ...
    'draw_end_idx', NaN, ...
    'score_total', -inf, ...
    'score_corr', -inf, ...
    'score_sep', -inf, ...
    'score_peak', -inf, ...
    'draw_len', NaN), num_cand, 1);

for i = 1:num_cand
    letter_name = candidate_letters{i};
    [stages, total_dur] = letter_stages(letter_name);
    draw_len = max(1, round(total_dur * cfg.sim.sampling_rate));

    force_nominal_start = strcmp(mode_name, 'known_letter') && cfg.sim.lock_start_when_known;
    c = evaluate_candidate(energy, stages, draw_len, cfg, force_nominal_start);
    cand_rows(i).letter = string(letter_name);
    cand_rows(i).start_idx = c.start_idx;
    cand_rows(i).draw_end_idx = c.draw_end_idx;
    cand_rows(i).score_total = c.score_total;
    cand_rows(i).score_corr = c.score_corr;
    cand_rows(i).score_sep = c.score_sep;
    cand_rows(i).score_peak = c.score_peak;
    cand_rows(i).draw_len = draw_len;
end

[~, best_i] = max([cand_rows.score_total]);
best = cand_rows(best_i);
best_letter = char(best.letter);
[best_stages, ~] = letter_stages(best_letter);

traj_full_x = nan(N, 1);
traj_full_y = nan(N, 1);
pen_mask = false(N, 1);

start_idx = min(max(1, round(best.start_idx)), N);
draw_end = min(N, round(best.draw_end_idx));
if draw_end < start_idx
    draw_end = start_idx;
end

for t_idx = start_idx:draw_end
    dt = (t_idx - start_idx) / cfg.sim.sampling_rate;
    [pos, pen] = sample_stage_pos(best_stages, dt);
    traj_full_x(t_idx) = pos(1);
    traj_full_y(t_idx) = pos(2);
    pen_mask(t_idx) = pen;
end

traj_t = find(isfinite(traj_full_x) & isfinite(traj_full_y));
traj_x = traj_full_x(traj_t);
traj_y = traj_full_y(traj_t);

traj_conf_full = zeros(N, 1);
if ~isempty(traj_t)
    conf_local = cfg.sim.conf_pen_weight * double(pen_mask(traj_t)) + ...
        (1 - cfg.sim.conf_pen_weight) * energy(traj_t);
    conf_local = conf_local .* sigmoid(best.score_total, cfg.sim.conf_score_scale);
    traj_conf_full(traj_t) = min(max(conf_local, 0), 1);
end
traj_conf = traj_conf_full(traj_t);

debug_info = struct();
debug_info.mode = mode_name;
debug_info.best_letter = best_letter;
debug_info.start_idx = start_idx;
debug_info.draw_end_idx = draw_end;
debug_info.energy = energy;
debug_info.pen_mask = pen_mask;
debug_info.traj_full_x = traj_full_x;
debug_info.traj_full_y = traj_full_y;
debug_info.candidates = struct2table(cand_rows);
debug_info.config = cfg;
end

% -------------------------------------------------------------------------
function c = evaluate_candidate(energy, stages, draw_len, cfg, force_nominal_start)
N = numel(energy);
start_nom = round(cfg.sim.start_ratio * N);
start_nom = min(max(start_nom, 1), N);

if force_nominal_start
    shift_vals = 0;
else
    shift_vals = -cfg.sim.max_start_shift:cfg.sim.start_shift_step:cfg.sim.max_start_shift;
end

best_score = -inf;
best_corr = -inf;
best_sep = -inf;
best_peak = -inf;
best_start = start_nom;
best_end = min(N, start_nom + draw_len);

for sh = shift_vals
    s = min(max(1, start_nom + sh), N);
    e = min(N, s + draw_len);
    if e < s
        continue;
    end

    mask = template_pen_mask(stages, N, s, e, cfg.sim.sampling_rate);
    if nnz(mask) < 6
        continue;
    end

    m = smoothdata(double(mask), 'movmean', cfg.sim.mask_smooth_win);
    m = normalize_vec(m(:));

    corr_val = corr(energy, m, 'Rows', 'pairwise');
    if ~isfinite(corr_val)
        corr_val = -1;
    end

    in_mean = mean(energy(mask), 'omitnan');
    out_mean = mean(energy(~mask), 'omitnan');
    sep_val = in_mean - out_mean;
    if ~isfinite(sep_val)
        sep_val = -1;
    end

    peak_val = peak_alignment_score(energy, mask, cfg.sim.peak_top_ratio);

    score = cfg.sim.score_corr_weight * corr_val + ...
        cfg.sim.score_sep_weight * sep_val + ...
        cfg.sim.score_peak_weight * peak_val;

    if score > best_score
        best_score = score;
        best_corr = corr_val;
        best_sep = sep_val;
        best_peak = peak_val;
        best_start = s;
        best_end = e;
    end
end

c = struct();
c.start_idx = best_start;
c.draw_end_idx = best_end;
c.score_total = best_score;
c.score_corr = best_corr;
c.score_sep = best_sep;
c.score_peak = best_peak;
end

% -------------------------------------------------------------------------
function p = peak_alignment_score(energy, mask, top_ratio)
e = energy(:);
m = mask(:);
good = isfinite(e);
e(~good) = 0;

k = max(1, round(top_ratio * numel(e)));
[~, idx] = maxk(e, k);
p = mean(double(m(idx)));
if ~isfinite(p)
    p = 0;
end
end

% -------------------------------------------------------------------------
function [candidate_letters, mode_name] = resolve_candidates(cfg)
candidate_letters = {};

target = '';
if isfield(cfg.sim, 'target_letter') && ~isempty(cfg.sim.target_letter)
    target = normalize_letter(cfg.sim.target_letter);
end

if ~isempty(target)
    candidate_letters = {target};
    mode_name = 'known_letter';
    return;
end

if ~cfg.sim.auto_detect_letter
    error('run_gesture_analysis_sim_oracle:missingTarget', ...
        'target_letter is empty and auto_detect_letter is disabled.');
end

raw = cfg.sim.candidate_letters;
if ischar(raw) || isstring(raw)
    raw = cellstr(raw);
end

for i = 1:numel(raw)
    try
        candidate_letters{end+1} = normalize_letter(raw{i}); %#ok<AGROW>
    catch
        % Ignore invalid entries.
    end
end

candidate_letters = unique(candidate_letters, 'stable');
mode_name = 'auto_letter';
end

% -------------------------------------------------------------------------
function letter = normalize_letter(letter_in)
letter = upper(strtrim(char(letter_in)));
if strcmp(letter, 'STAR')
    letter = 'Star';
    return;
end

valid = {'A','B','M','L','X','Z','N'};
if any(strcmp(letter, valid))
    return;
end

error('Invalid letter: %s', char(letter_in));
end

% -------------------------------------------------------------------------
function energy = build_energy(step1_res, cfg)
V = step1_res.volatility_matrix;
if isempty(V)
    energy = zeros(numel(step1_res.t_grid), 1);
    return;
end

V = max(V, 0);
W = ones(1, size(V, 2));
if isfield(step1_res, 'elevation_matrix') && ...
        ~isempty(step1_res.elevation_matrix) && ...
        all(size(step1_res.elevation_matrix) == size(V))
    E = max(step1_res.elevation_matrix, 0);
    W = mean(E, 1, 'omitnan');
    if all(~isfinite(W)) || all(W <= 0)
        W = ones(size(W));
    end
end

W = W(:)' ./ max(sum(W, 'omitnan'), eps);
energy = (V * W(:));
energy = normalize_vec(energy);
energy = smoothdata(energy, 'movmean', cfg.sim.energy_smooth_win);
energy = normalize_vec(energy);
end

% -------------------------------------------------------------------------
function [pos, pen] = sample_stage_pos(stages, dt)
pos = [NaN, NaN];
pen = false;
elapsed = 0;

for k = 1:size(stages, 1)
    p1 = stages{k, 1};
    p2 = stages{k, 2};
    dur = stages{k, 3};
    pen_k = stages{k, 4};
    if dt <= (elapsed + dur)
        alpha = max(0, min(1, (dt - elapsed) / max(dur, eps)));
        pos = p1 + (p2 - p1) * alpha;
        pen = logical(pen_k);
        return;
    end
    elapsed = elapsed + dur;
end

pos = stages{end, 2};
pen = logical(stages{end, 4});
end

% -------------------------------------------------------------------------
function mask = template_pen_mask(stages, N, start_idx, end_idx, sampling_rate)
mask = false(N, 1);
s = max(1, start_idx);
e = min(N, end_idx);
if e < s
    return;
end

for t_idx = s:e
    dt = (t_idx - s) / sampling_rate;
    [~, pen] = sample_stage_pos(stages, dt);
    mask(t_idx) = pen;
end
end

% -------------------------------------------------------------------------
function [stages, total_dur] = letter_stages(letter_name)
switch upper(string(letter_name))
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
        error('Unsupported letter: %s', char(letter_name));
end

total_dur = 0;
for k = 1:size(stages, 1)
    total_dur = total_dur + stages{k, 3};
end
end

% -------------------------------------------------------------------------
function y = sigmoid(x, scale)
y = 1 ./ (1 + exp(-x / max(scale, eps)));
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
function cfg = default_cfg()
cfg = struct();
cfg.sim = struct();
cfg.sim.target_letter = '';
cfg.sim.auto_detect_letter = true;
cfg.sim.candidate_letters = {'A','B','M','Star','L','X','Z','N'};
cfg.sim.start_ratio = 0.30;
cfg.sim.sampling_rate = 25;
cfg.sim.max_start_shift = 90;
cfg.sim.start_shift_step = 1;
cfg.sim.energy_smooth_win = 11;
cfg.sim.mask_smooth_win = 9;
cfg.sim.score_corr_weight = 0.65;
cfg.sim.score_sep_weight = 0.30;
cfg.sim.score_peak_weight = 0.05;
cfg.sim.peak_top_ratio = 0.10;
cfg.sim.conf_pen_weight = 0.60;
cfg.sim.conf_score_scale = 0.5;
cfg.sim.lock_start_when_known = true;
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
