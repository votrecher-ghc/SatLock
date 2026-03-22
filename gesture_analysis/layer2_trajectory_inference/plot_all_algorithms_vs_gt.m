function plot_all_algorithms_vs_gt(alg_results, t_grid, target_letter, span_cfg, out_dir, show_plot)
% PLOT_ALL_ALGORITHMS_VS_GT
% 功能：
%   将“每种算法”的轨迹分别与 GroundTruth 做对比绘图。
%
% 输入：
%   alg_results : 结构体数组，每项包含
%                 .name (算法名)
%                 .x / .y / .t (轨迹)
%                 .conf (可选)
%   t_grid      : 统一时间轴
%   target_letter : GroundTruth 字母
%   span_cfg    : 范围约束配置（max_span_x/max_span_y）
%   out_dir     : 输出目录（png 保存位置）
%   show_plot   : true=弹窗显示，false=仅保存

if nargin < 6
    show_plot = true;
end
if nargin < 5 || isempty(out_dir)
    out_dir = fullfile('gesture_analysis', 'results');
end
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end
if nargin < 4 || isempty(span_cfg)
    span_cfg = struct();
end
if ~isfield(span_cfg, 'max_span_x')
    span_cfg.max_span_x = 0.50;
end
if ~isfield(span_cfg, 'max_span_y')
    span_cfg.max_span_y = 0.50;
end

N = numel(t_grid);
[gt_x, gt_y, gt_pen] = build_ground_truth(target_letter, N, span_cfg);

for i = 1:numel(alg_results)
    item = alg_results(i);
    if ~isfield(item, 'name') || isempty(item.name)
        item.name = sprintf('alg_%d', i);
    end

    if ~isfield(item, 'conf')
        item.conf = [];
    end

    [est_x, est_y] = extract_plot_series(item, t_grid);

    fig_vis = 'off';
    if show_plot
        fig_vis = 'on';
    end

    f = figure('Name', sprintf('GT vs %s', char(item.name)), ...
        'Color', 'w', 'Position', [120, 80, 900, 760], 'Visible', fig_vis);

    hold on;
    grid on;
    axis equal;

    plot(0, 0, '^k', 'MarkerFaceColor', 'k', 'DisplayName', '接收机');

    idx_up = ~gt_pen & isfinite(gt_x) & isfinite(gt_y);
    idx_dn = gt_pen & isfinite(gt_x) & isfinite(gt_y);
    plot(gt_x(idx_up), gt_y(idx_up), '--', 'Color', [0.70 0.70 0.70], ...
        'LineWidth', 1.0, 'DisplayName', 'GroundTruth 抬笔');
    plot(gt_x(idx_dn), gt_y(idx_dn), '-', 'Color', [0.10 0.45 0.95], ...
        'LineWidth', 3.0, 'DisplayName', 'GroundTruth 落笔');

    plot(est_x, est_y, '-', 'Color', [0.90 0.25 0.20], ...
        'LineWidth', 2.1, 'DisplayName', '算法估计轨迹');

    idx_est = find(isfinite(est_x) & isfinite(est_y));
    if ~isempty(idx_est)
        plot(est_x(idx_est(1)), est_y(idx_est(1)), 'go', ...
            'MarkerFaceColor', 'g', 'DisplayName', '估计起点');
        plot(est_x(idx_est(end)), est_y(idx_est(end)), 'mo', ...
            'MarkerFaceColor', 'm', 'DisplayName', '估计终点');
    end

    xlabel('East (m)');
    ylabel('North (m)');
    title(sprintf('%s vs GroundTruth (%.0fcm x %.0fcm)', ...
        char(item.name), 100 * span_cfg.max_span_x, 100 * span_cfg.max_span_y));
    legend('Location', 'best');

    fname = sprintf('%s_vs_groundtruth.png', sanitize_name(item.name));
    fpath = fullfile(out_dir, fname);
    exportgraphics(f, fpath, 'Resolution', 150);

    if ~show_plot
        close(f);
    end
end
end

% -------------------------------------------------------------------------
function [xp, yp] = extract_plot_series(item, t_grid)
xp = [];
yp = [];

if ~isfield(item, 'x') || ~isfield(item, 'y') || isempty(item.x) || isempty(item.y)
    return;
end

x = item.x(:);
y = item.y(:);
n = min(numel(x), numel(y));
x = x(1:n);
y = y(1:n);

% 优先按 t 对齐排序，避免时序乱序；不再插入完整时间轴 NaN。
if isfield(item, 't') && ~isempty(item.t)
    t = item.t(:);
    t = t(1:min(numel(t), n));
    idx = map_time_to_index(t, t_grid, n);
    idx = idx(1:min(numel(idx), n));
    m = min([numel(x), numel(y), numel(idx)]);
    x = x(1:m);
    y = y(1:m);
    idx = idx(1:m);
    keep = isfinite(x) & isfinite(y) & isfinite(idx);
    x = x(keep);
    y = y(keep);
    idx = idx(keep);
    [~, ord] = sort(idx, 'ascend');
    xp = x(ord);
    yp = y(ord);
else
    keep = isfinite(x) & isfinite(y);
    xp = x(keep);
    yp = y(keep);
end
end

% -------------------------------------------------------------------------
function idx = map_time_to_index(t, t_grid, n_ref)
N = numel(t_grid);

if isnumeric(t)
    t_num = t(:);
    if all(isfinite(t_num)) && all(t_num >= 1) && all(t_num <= N)
        idx = round(t_num);
        return;
    end
    if numel(t_num) == N
        idx = (1:N)';
        return;
    end
    idx = round(linspace(1, N, max(1, n_ref)))';
    return;
end

if isdatetime(t)
    tg = t_grid(:);
    tt = t(:);
    tg_sec = posixtime(tg);
    tt_sec = posixtime(tt);
    idx = round(interp1(tg_sec, 1:N, tt_sec, 'nearest', 'extrap'));
    idx = min(max(idx, 1), N);
    return;
end

idx = round(linspace(1, N, max(1, n_ref)))';
end

% -------------------------------------------------------------------------
function s = sanitize_name(name_in)
s = char(name_in);
s = strrep(s, ' ', '_');
s = regexprep(s, '[^a-zA-Z0-9_\-]', '_');
if isempty(s)
    s = 'algorithm';
end
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
s = min([sx, sy, 1.0]);

for i = 1:size(stages_in, 1)
    p1 = stages_in{i, 1};
    p2 = stages_in{i, 2};
    stages_out{i, 1} = (p1 - center) * s;
    stages_out{i, 2} = (p2 - center) * s;
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
