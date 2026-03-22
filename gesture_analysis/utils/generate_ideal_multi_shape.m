% =========================================================================
% generate_ideal_multi_shape.m
% 功能：
%   在真实 obs/nav 基础上注入理想手势遮挡信号，并绘制 GroundTruth 轨迹。
%
% 关键点：
%   1) 使用分段直线（含停顿段）描述字母轨迹；
%   2) 支持“感知范围约束”，默认左右/上下总跨度均为 50cm；
%   3) 保留可视化，便于与推演结果直接对比。
%
% 调用：
%   obs_sim = generate_ideal_multi_shape(obs_data, nav_data, target_letter);
%   obs_sim = generate_ideal_multi_shape(obs_data, nav_data, target_letter, sim_cfg);
% =========================================================================
function obs_data = generate_ideal_multi_shape(obs_data, nav_data, target_letter, sim_cfg)

addpath(genpath('calculate_clock_bias_and_positon'));
addpath(genpath('nav_parse'));

if nargin < 4 || isempty(sim_cfg)
    sim_cfg = struct();
end

% -------------------- 仿真参数 --------------------
SIM = struct();
SIM.baseline_db = get_cfg(sim_cfg, 'baseline_db', 45);      % 基线信号强度
SIM.drop_depth_db = get_cfg(sim_cfg, 'drop_depth_db', 15);  % 遮挡下跌深度
SIM.noise_sigma = get_cfg(sim_cfg, 'noise_sigma', 0.02);    % 微噪声

SIM.gesture_height = get_cfg(sim_cfg, 'gesture_height', 0.30); % 手势平面高度
SIM.arm_width = get_cfg(sim_cfg, 'arm_width', 0.15);           % 手臂遮挡宽度
SIM.body_pos = get_cfg(sim_cfg, 'body_pos', [0.0, -1.0]);      % 身体参考点

SIM.start_ratio = get_cfg(sim_cfg, 'start_ratio', 0.30);
SIM.sampling_rate = get_cfg(sim_cfg, 'sampling_rate', 25);
SIM.plot = get_cfg(sim_cfg, 'plot', true);

% 感知范围约束（总跨度）
SIM.max_span_x = get_cfg(sim_cfg, 'max_span_x', 0.50); % 左右 50cm
SIM.max_span_y = get_cfg(sim_cfg, 'max_span_y', 0.50); % 上下 50cm

fprintf('--> 启动仿真 V7.2: 目标 [%s]，范围约束 %.2fm x %.2fm\n', ...
    char(target_letter), SIM.max_span_x, SIM.max_span_y);

% -------------------- 字母模板 --------------------
if exist('gesture_template_library', 'file') == 2
    span_cfg = struct('max_span_x', SIM.max_span_x, 'max_span_y', SIM.max_span_y);
    [stages, span_meta] = gesture_template_library('stages', target_letter, span_cfg);
else
    stages = letter_stages(target_letter);
    [stages, span_meta] = normalize_stages_to_span(stages, SIM.max_span_x, SIM.max_span_y);
end
fprintf('    轨迹缩放系数: %.3f, 约束后跨度: [%.3fm, %.3fm]\n', ...
    span_meta.scale, span_meta.span_x, span_meta.span_y);

total_sim_duration = 0;
for k = 1:size(stages, 1)
    total_sim_duration = total_sim_duration + stages{k, 3};
end

% -------------------- 数据准备 --------------------
ideal_obs = obs_data;
num_samples = numel(ideal_obs);
start_idx = round(num_samples * SIM.start_ratio);
end_idx = min(num_samples, start_idx + round(total_sim_duration * SIM.sampling_rate));

all_sat_ids = {};
for i = 1:min(100, num_samples)
    if ~isempty(ideal_obs(i).data)
        all_sat_ids = [all_sat_ids, fieldnames(ideal_obs(i).data)']; %#ok<AGROW>
    end
end
unique_sat_ids = unique(all_sat_ids);
valid_sats = {};
for i = 1:numel(unique_sat_ids)
    sid = unique_sat_ids{i};
    if ismember(sid(1), ['G', 'C', 'R', 'E', 'J'])
        valid_sats{end+1} = sid; %#ok<AGROW>
    end
end

% 估计接收机平均坐标（用于卫星投影）
rec_pos_acc = [0, 0, 0];
count = 0;
for t = start_idx:10:end_idx
    try
        [rp, ~, ~] = calculate_receiver_position(obs_data, nav_data, t);
        rec_pos_acc = rec_pos_acc + rp;
        count = count + 1;
    catch
    end
end
if count == 0
    warning('接收机位置估计失败，返回原始数据。');
    obs_data = obs_data;
    return;
end

rec_pos_mean = rec_pos_acc / count;
[lat0, lon0, alt0] = ecef2geodetic(rec_pos_mean(1), rec_pos_mean(2), rec_pos_mean(3));

% -------------------- 信号注入 --------------------
gt_trace_x = nan(num_samples, 1);
gt_trace_y = nan(num_samples, 1);
gt_pen_down = false(num_samples, 1);

for t_idx = 1:num_samples
    is_pen_down = false;
    current_hand_pos = [NaN, NaN];

    if t_idx >= start_idx && t_idx <= end_idx
        dt = (t_idx - start_idx) / SIM.sampling_rate;
        [current_hand_pos, is_pen_down] = sample_stage_pos(stages, dt);
        gt_trace_x(t_idx) = current_hand_pos(1);
        gt_trace_y(t_idx) = current_hand_pos(2);
        gt_pen_down(t_idx) = is_pen_down;
    end

    if isempty(ideal_obs(t_idx).data)
        continue;
    end

    try
        [~, ~, sat_states] = calculate_receiver_position(obs_data, nav_data, t_idx);
    catch
        continue;
    end

    for s = 1:numel(valid_sats)
        sid = valid_sats{s};
        if ~isfield(ideal_obs(t_idx).data, sid)
            continue;
        end

        sim_val = SIM.baseline_db + randn() * SIM.noise_sigma;

        if is_pen_down && all(isfinite(current_hand_pos)) && isfield(sat_states, sid)
            sat_p = sat_states.(sid).position;
            [e, n, u] = ecef2enu( ...
                sat_p(1) - rec_pos_mean(1), ...
                sat_p(2) - rec_pos_mean(2), ...
                sat_p(3) - rec_pos_mean(3), ...
                lat0, lon0, alt0);

            if u > 0
                scale = SIM.gesture_height / u;
                P = [scale * e, scale * n];
                A = current_hand_pos;
                B = SIM.body_pos;

                vec_AB = B - A;
                vec_AP = P - A;
                len_sq = sum(vec_AB.^2);
                if len_sq > 0
                    t_proj = max(0, min(1, dot(vec_AP, vec_AB) / len_sq));
                    dist_to_arm = norm(P - (A + t_proj * vec_AB));
                    if dist_to_arm < SIM.arm_width
                        sim_val = SIM.baseline_db - SIM.drop_depth_db + randn() * SIM.noise_sigma;
                    end
                end
            end
        end

        snr_struct = ideal_obs(t_idx).data.(sid).snr;
        fields = fieldnames(snr_struct);
        for f = 1:numel(fields)
            ideal_obs(t_idx).data.(sid).snr.(fields{f}) = sim_val;
        end
    end
end

obs_data = ideal_obs;

% -------------------- GroundTruth 绘图 --------------------
if SIM.plot
    fprintf('--> 生成 GroundTruth 轨迹图...\n');
    figure('Name', sprintf('Ideal GroundTruth [%s]', char(target_letter)), ...
        'Position', [100, 100, 620, 620], 'Color', 'w');
    ax = axes;
    hold(ax, 'on');
    grid(ax, 'on');
    axis(ax, 'equal');
    xlabel('East (m)');
    ylabel('North (m)');

    plot(ax, 0, 0, '^', 'MarkerSize', 12, 'MarkerFaceColor', 'k', ...
        'MarkerEdgeColor', 'k', 'DisplayName', 'Receiver');
    plot(ax, SIM.body_pos(1), SIM.body_pos(2), 'bs', 'MarkerSize', 10, ...
        'MarkerFaceColor', 'b', 'DisplayName', 'Body Ref');

    plot_x_up = gt_trace_x;
    plot_y_up = gt_trace_y;
    plot_x_up(gt_pen_down) = NaN;
    plot_y_up(gt_pen_down) = NaN;

    plot_x_dn = gt_trace_x;
    plot_y_dn = gt_trace_y;
    plot_x_dn(~gt_pen_down) = NaN;
    plot_y_dn(~gt_pen_down) = NaN;

    plot(ax, plot_x_up, plot_y_up, '--', 'Color', [0.65 0.65 0.65], ...
        'LineWidth', 1.1, 'DisplayName', 'GT Pen Up');
    plot(ax, plot_x_dn, plot_y_dn, '-', 'Color', [0.10 0.45 0.95], ...
        'LineWidth', 2.8, 'DisplayName', 'GT Pen Down');

    f_idx = find(isfinite(gt_trace_x) & isfinite(gt_trace_y), 1, 'first');
    l_idx = find(isfinite(gt_trace_x) & isfinite(gt_trace_y), 1, 'last');
    if ~isempty(f_idx)
        plot(ax, gt_trace_x(f_idx), gt_trace_y(f_idx), 'go', ...
            'MarkerFaceColor', 'g', 'MarkerSize', 9, 'DisplayName', 'Start');
    end
    if ~isempty(l_idx)
        plot(ax, gt_trace_x(l_idx), gt_trace_y(l_idx), 'rs', ...
            'MarkerFaceColor', 'r', 'MarkerSize', 9, 'DisplayName', 'End');
    end

    title({sprintf('Ideal GroundTruth (Letter %s)', char(target_letter)), ...
        sprintf('Range Constraint: %.0fcm x %.0fcm', 100 * SIM.max_span_x, 100 * SIM.max_span_y)});
    legend('Location', 'best');

    valid_x = gt_trace_x(isfinite(gt_trace_x));
    valid_y = gt_trace_y(isfinite(gt_trace_y));
    if ~isempty(valid_x)
        m_range = max([max(abs(valid_x)), max(abs(valid_y)), 0.30]);
    else
        m_range = 0.30;
    end
    xlim(ax, [-m_range * 1.25, m_range * 1.25]);
    ylim(ax, [-m_range * 1.25, m_range * 1.25] + SIM.body_pos(2) * 0.1);
end

fprintf('✓ 仿真注入完成。\n');
end

% =========================================================================
% 局部函数
% =========================================================================
function v = get_cfg(s, k, v0)
if isstruct(s) && isfield(s, k)
    v = s.(k);
else
    v = v0;
end
end

% -------------------------------------------------------------------------
function [stages, meta] = normalize_stages_to_span(stages_in, max_span_x, max_span_y)
stages = stages_in;
meta = struct('scale', 1.0, 'span_x', NaN, 'span_y', NaN);
if isempty(stages_in)
    return;
end

pts = zeros(0, 2);
for i = 1:size(stages_in, 1)
    pts(end+1, :) = stages_in{i, 1}; %#ok<AGROW>
    pts(end+1, :) = stages_in{i, 2}; %#ok<AGROW>
end

min_xy = min(pts, [], 1);
max_xy = max(pts, [], 1);
span = max_xy - min_xy;
center = 0.5 * (min_xy + max_xy);

sx = max_span_x / max(span(1), eps);
sy = max_span_y / max(span(2), eps);
s = min([sx, sy, 1.0]); % 只缩小，不放大

for i = 1:size(stages_in, 1)
    p1 = stages_in{i, 1};
    p2 = stages_in{i, 2};
    stages{i, 1} = (p1 - center) * s;
    stages{i, 2} = (p2 - center) * s;
end

pts2 = zeros(0, 2);
for i = 1:size(stages, 1)
    pts2(end+1, :) = stages{i, 1}; %#ok<AGROW>
    pts2(end+1, :) = stages{i, 2}; %#ok<AGROW>
end
min2 = min(pts2, [], 1);
max2 = max(pts2, [], 1);
span2 = max2 - min2;

meta.scale = s;
meta.span_x = span2(1);
meta.span_y = span2(2);
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
        a = max(0, min(1, (dt - elapsed) / max(dur, eps)));
        pos = p1 + (p2 - p1) * a;
        pen = logical(pen_k);
        return;
    end
    elapsed = elapsed + dur;
end

pos = stages{end, 2};
pen = logical(stages{end, 4});
end

% -------------------------------------------------------------------------
function stages = letter_stages(letter_in)
letter = upper(strtrim(char(letter_in)));
if strcmp(letter, 'STAR')
    letter = 'STAR';
end

switch letter
    case 'A'
        P1 = [-0.40, -0.50]; P2 = [0.00, 0.50]; P3 = [0.40, -0.50];
        P4 = [-0.20, -0.10]; P5 = [0.20, -0.10];
        stages = {
            P1, P2, 1.5, true;
            P2, P2, 3.0, true;
            P2, P3, 1.5, true;
            P3, P4, 3.0, false;
            P4, P5, 1.5, true
        };

    case 'B'
        P1 = [-0.40, -1.00]; P2 = [-0.40, 1.00]; P3 = [1.50, 0.40];
        P4 = [-0.40, 0.00]; P5 = [1.50, 0.00]; P6 = [-0.40, -1.00];
        stages = {
            P1, P2, 1.5, true; P2, P2, 3.0, true;
            P2, P3, 1.5, true; P3, P3, 3.0, true;
            P3, P4, 3.0, true; P4, P4, 3.0, true;
            P4, P5, 1.5, true; P5, P5, 3.0, true;
            P5, P6, 1.5, true
        };

    case 'M'
        P1 = [-0.40, -0.40]; P2 = [-0.40, 0.40]; P3 = [0.00, 0.00];
        P4 = [0.40, 0.40]; P5 = [0.40, -0.40];
        stages = {
            P1, P2, 1.5, true; P2, P2, 3.0, true;
            P2, P3, 1.5, true; P3, P3, 3.0, true;
            P3, P4, 1.5, true; P4, P4, 3.0, true;
            P4, P5, 1.5, true
        };

    case 'STAR'
        P1 = [-0.30, -0.45]; P2 = [0.00, 0.55]; P3 = [0.30, -0.45];
        P4 = [-0.48, 0.15]; P5 = [0.48, 0.15];
        stages = {
            P1, P2, 1.5, true; P2, P2, 3.0, true;
            P2, P3, 1.5, true; P3, P3, 3.0, true;
            P3, P4, 1.5, true; P4, P4, 3.0, true;
            P4, P5, 1.5, true; P5, P5, 3.0, true;
            P5, P1, 1.5, true
        };

    case 'L'
        P1 = [-0.3, 0.4]; P2 = [-0.3, -0.4]; P3 = [0.3, -0.4];
        stages = {
            P1, P2, 1.5, true;
            P2, P2, 3.0, true;
            P2, P3, 1.5, true
        };

    case 'X'
        P1 = [-0.3, 0.4]; P2 = [0.3, -0.4]; P3 = [0.3, 0.4]; P4 = [-0.3, -0.4];
        stages = {
            P1, P2, 1.5, true;
            P2, P3, 3.0, false;
            P3, P4, 1.5, true
        };

    case 'Z'
        P1 = [-0.3, 0.4]; P2 = [0.3, 0.4]; P3 = [-0.3, -0.4]; P4 = [0.7, -0.4];
        stages = {
            P1, P2, 1.5, true;
            P2, P2, 3.0, true;
            P2, P3, 1.5, true;
            P3, P3, 3.0, true;
            P3, P4, 1.5, true
        };

    case 'N'
        P1 = [-0.3, -0.4]; P2 = [-0.3, 0.4]; P3 = [0.3, -0.4]; P4 = [0.3, 0.4];
        stages = {
            P1, P2, 1.5, true;
            P2, P2, 3.0, true;
            P2, P3, 1.5, true;
            P3, P3, 3.0, true;
            P3, P4, 1.5, true
        };

    otherwise
        error('未定义的字母: %s', char(letter_in));
end
end

