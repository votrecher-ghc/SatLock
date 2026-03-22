function plot_layer2_compare(traj_x, traj_y, traj_t, traj_conf, t_grid, target_letter, span_cfg)
% PLOT_LAYER2_COMPARE
% 第二层结果可视化：只对比 GroundTruth 与新算法估计轨迹。
%
% 注意：
%   target_letter 仅用于“仿真验证可视化”，不参与轨迹推演。

if nargin < 6
    target_letter = '';
end
if nargin < 7 || isempty(span_cfg)
    span_cfg = struct();
end
if ~isfield(span_cfg, 'max_span_x')
    span_cfg.max_span_x = 0.50;
end
if ~isfield(span_cfg, 'max_span_y')
    span_cfg.max_span_y = 0.50;
end

N = numel(t_grid);
[est_full_x, est_full_y, ~] = to_full_series(traj_x, traj_y, traj_t, traj_conf, N);

has_gt = ~isempty(target_letter);
if has_gt
    [gt_x, gt_y, gt_pen] = build_ground_truth(target_letter, N, span_cfg);
end

figure('Name', 'GroundTruth vs 新算法轨迹', 'Color', 'w', 'Position', [120, 80, 900, 760]);
hold on;
grid on;
axis equal;
plot(0, 0, '^k', 'MarkerFaceColor', 'k', 'DisplayName', '接收机');

if has_gt
    idx_up = ~gt_pen & isfinite(gt_x) & isfinite(gt_y);
    idx_dn = gt_pen & isfinite(gt_x) & isfinite(gt_y);
    plot(gt_x(idx_up), gt_y(idx_up), '--', 'Color', [0.70 0.70 0.70], ...
        'LineWidth', 1.0, 'DisplayName', 'GroundTruth 抬笔');
    plot(gt_x(idx_dn), gt_y(idx_dn), '-', 'Color', [0.10 0.45 0.95], ...
        'LineWidth', 3.0, 'DisplayName', 'GroundTruth 落笔');
end

plot(est_full_x, est_full_y, '-', 'Color', [0.90 0.25 0.20], ...
    'LineWidth', 2.2, 'DisplayName', '新算法估计轨迹');

idx_est = find(isfinite(est_full_x) & isfinite(est_full_y));
if ~isempty(idx_est)
    plot(est_full_x(idx_est(1)), est_full_y(idx_est(1)), 'go', ...
        'MarkerFaceColor', 'g', 'DisplayName', '估计起点');
    plot(est_full_x(idx_est(end)), est_full_y(idx_est(end)), 'mo', ...
        'MarkerFaceColor', 'm', 'DisplayName', '估计终点');
end

xlabel('East (m)');
ylabel('North (m)');
title(sprintf('轨迹对比（感知范围 %.0fcm x %.0fcm）', ...
    100 * span_cfg.max_span_x, 100 * span_cfg.max_span_y));
legend('Location', 'best');
end

% -------------------------------------------------------------------------
function [fx, fy, fc] = to_full_series(x, y, t_idx, c, N)
fx = nan(N, 1);
fy = nan(N, 1);
fc = nan(N, 1);

if isempty(x) || isempty(y) || isempty(t_idx)
    return;
end

x = x(:);
y = y(:);
t_idx = round(t_idx(:));
if isempty(c)
    c = nan(size(t_idx));
else
    c = c(:);
end

n = min([numel(x), numel(y), numel(t_idx), numel(c)]);
x = x(1:n);
y = y(1:n);
t_idx = t_idx(1:n);
c = c(1:n);

keep = isfinite(x) & isfinite(y) & isfinite(t_idx) & (t_idx >= 1) & (t_idx <= N);
x = x(keep);
y = y(keep);
t_idx = t_idx(keep);
c = c(keep);

[u_idx, ia] = unique(t_idx, 'stable');
fx(u_idx) = x(ia);
fy(u_idx) = y(ia);
fc(u_idx) = c(ia);
end

% -------------------------------------------------------------------------
function [gt_x, gt_y, gt_pen_down] = build_ground_truth(letter, num_samples, span_cfg)
if exist('gesture_template_library', 'file') == 2
    [gt_x, gt_y, gt_pen_down] = gesture_template_library('groundtruth', letter, num_samples, span_cfg);
    return;
end

sampling_rate = 25;
stages = letter_stages(letter);
stages = normalize_stages_to_span(stages, span_cfg.max_span_x, span_cfg.max_span_y);

gt_x = nan(num_samples, 1);
gt_y = nan(num_samples, 1);
gt_pen_down = false(num_samples, 1);

total_dur = 0;
for k = 1:size(stages, 1)
    total_dur = total_dur + stages{k, 3};
end
start_idx = round(num_samples * 0.3);
end_idx = min(num_samples, start_idx + round(total_dur * sampling_rate));

for t_idx = start_idx:end_idx
    dt = (t_idx - start_idx) / sampling_rate;
    elapsed = 0;
    for k = 1:size(stages, 1)
        p1 = stages{k, 1};
        p2 = stages{k, 2};
        dur = stages{k, 3};
        pen = stages{k, 4};
        if dt <= (elapsed + dur)
            alpha = max(0, min(1, (dt - elapsed) / max(dur, eps)));
            pos = p1 + (p2 - p1) * alpha;
            gt_x(t_idx) = pos(1);
            gt_y(t_idx) = pos(2);
            gt_pen_down(t_idx) = logical(pen);
            break;
        end
        elapsed = elapsed + dur;
    end
end
end

% -------------------------------------------------------------------------
function stages = letter_stages(letter)
switch upper(string(letter))
    case "A"
        P1 = [-0.40, -0.50]; P2 = [0.00, 0.50]; P3 = [0.40, -0.50];
        P4 = [-0.20, -0.10]; P5 = [0.20, -0.10];
        stages = {P1,P2,1.5,true; P2,P2,3.0,true; P2,P3,1.5,true; P3,P4,3.0,false; P4,P5,1.5,true};
    case "B"
        P1 = [-0.40,-1.00]; P2 = [-0.40,1.00]; P3 = [1.50,0.40];
        P4 = [-0.40,0.00]; P5 = [1.50,0.00]; P6 = [-0.40,-1.00];
        stages = {P1,P2,1.5,true; P2,P2,3.0,true; P2,P3,1.5,true; P3,P3,3.0,true; ...
                  P3,P4,3.0,true; P4,P4,3.0,true; P4,P5,1.5,true; P5,P5,3.0,true; P5,P6,1.5,true};
    case "M"
        P1=[-0.40,-0.40]; P2=[-0.40,0.40]; P3=[0.00,0.00]; P4=[0.40,0.40]; P5=[0.40,-0.40];
        stages = {P1,P2,1.5,true; P2,P2,3.0,true; P2,P3,1.5,true; P3,P3,3.0,true; ...
                  P3,P4,1.5,true; P4,P4,3.0,true; P4,P5,1.5,true};
    case "STAR"
        P1=[-0.30,-0.45]; P2=[0.00,0.55]; P3=[0.30,-0.45]; P4=[-0.48,0.15]; P5=[0.48,0.15];
        stages = {P1,P2,1.5,true; P2,P2,3.0,true; P2,P3,1.5,true; P3,P3,3.0,true; ...
                  P3,P4,1.5,true; P4,P4,3.0,true; P4,P5,1.5,true; P5,P5,3.0,true; P5,P1,1.5,true};
    case "L"
        P1=[-0.3,0.4]; P2=[-0.3,-0.4]; P3=[0.3,-0.4];
        stages = {P1,P2,1.5,true; P2,P2,3.0,true; P2,P3,1.5,true};
    case "X"
        P1=[-0.3,0.4]; P2=[0.3,-0.4]; P3=[0.3,0.4]; P4=[-0.3,-0.4];
        stages = {P1,P2,1.5,true; P2,P3,3.0,false; P3,P4,1.5,true};
    case "Z"
        P1=[-0.3,0.4]; P2=[0.3,0.4]; P3=[-0.3,-0.4]; P4=[0.7,-0.4];
        stages = {P1,P2,1.5,true; P2,P2,3.0,true; P2,P3,1.5,true; P3,P3,3.0,true; P3,P4,1.5,true};
    case "N"
        P1=[-0.3,-0.4]; P2=[-0.3,0.4]; P3=[0.3,-0.4]; P4=[0.3,0.4];
        stages = {P1,P2,1.5,true; P2,P2,3.0,true; P2,P3,1.5,true; P3,P3,3.0,true; P3,P4,1.5,true};
    otherwise
        error('不支持的字母: %s', char(letter));
end
end

% -------------------------------------------------------------------------
function stages_out = normalize_stages_to_span(stages_in, max_span_x, max_span_y)
stages_out = stages_in;
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
s = min([sx, sy, 1.0]); % 只缩小不放大

for i = 1:size(stages_in, 1)
    p1 = stages_in{i, 1};
    p2 = stages_in{i, 2};
    stages_out{i, 1} = (p1 - center) * s;
    stages_out{i, 2} = (p2 - center) * s;
end
end
