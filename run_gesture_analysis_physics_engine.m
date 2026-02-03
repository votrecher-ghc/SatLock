% =========================================================================
% run_gesture_analysis_physics_engine.m
% 功能: 物理引擎逆向匹配与几何约束融合算法 (Physics-Based Solver)
% 描述:
%   这是针对“手臂干扰”和“轨迹准确性”的根本性解决方案。
%   它不依赖数据的统计特征，而是直接逆向求解物理遮挡方程。
%
%   [核心流程]:
%   Step 1: 物理逆向求解 (Inverse Physics)
%       - 网格搜索: 遍历空间所有可能位置 P。
%       - 阴影合成: 对每个 P，构建虚拟手臂 (P->Body)，计算理论遮挡状态。
%       - 似然匹配: 寻找与实际 SNR 跌落图谱最吻合的 P，作为 Raw Point。
%       -> 彻底消除手臂对重心的拉扯，恢复绝对位置。
%
%   Step 2: 几何笔画精修 (Stroke Refinement)
%       - 速度/角度检测: 识别动作折点。
%       - PCA 线性化: 对每段笔画进行主成分投影，强制拉直。
%       -> 消除网格搜索带来的量化抖动，恢复完美的几何形状。
%
% [调用]:
%   [final_traj, debug_info] = run_gesture_analysis_physics_engine(obs_clean, nav_data, step1_res);
% =========================================================================

function [final_traj, debug_info] = run_gesture_analysis_physics_engine(obs_clean, nav_data, step1_res)

% --- 1. 数据与参数加载 ---
vol_mat = step1_res.volatility_matrix;
t_grid  = step1_res.t_grid;
valid_sats = step1_res.valid_sats;
segments = step1_res.segments;

if isempty(segments)
    fprintf('⚠️ 无有效分段，算法终止。\n');
    final_traj=[]; debug_info=[]; return;
end

addpath(genpath('calculate_clock_bias_and_positon'));
addpath(genpath('nav_parse'));

fprintf('--> [Ultimate Step 2] 启动物理引擎逆向求解 (Physics Engine Solver)...\n');

% === [关键参数] 必须与 generate_ideal_multi_shape 保持物理一致 ===
% 参考 中的定义
MODEL.body_pos       = [0.0, -1.0];   % 身体位置
MODEL.gesture_height = 0.20;          % 平面高度
MODEL.arm_width      = 0.40;          % 手臂宽度
% ==============================================================

% 网格搜索参数 (精度越高越准，但越慢)
GRID.range_x = [-0.6, 0.6];
GRID.range_y = [-0.6, 0.6];
GRID.step    = 0.02; % 2cm 网格

% 信号判决门限
PARA.hit_threshold   = 1.5; % dB

%% ================= [Part 1] 物理逆向求解 (Inverse Physics) =================
fprintf('    1. 执行全空间逆向阴影匹配 (解决手臂干扰)...\n');

% 1.1 构建网格
gx = GRID.range_x(1) : GRID.step : GRID.range_x(2);
gy = GRID.range_y(1) : GRID.step : GRID.range_y(2);
[XX, YY] = meshgrid(gx, gy);
grid_pts = [XX(:), YY(:)]; % M x 2
num_grid = size(grid_pts, 1);

% 1.2 预计算接收机位置 (取中间时刻)
mid_idx = round(length(obs_clean)/2);
[rec_xyz, ~, ~, lat0, lon0, alt0] = calculate_receiver_position(obs_clean, nav_data, mid_idx);

% 1.3 核心循环
num_samples = length(t_grid);
raw_x = NaN(num_samples, 1);
raw_y = NaN(num_samples, 1);
raw_conf = zeros(num_samples, 1); % 置信度

% 构建 Active Mask
is_active_mask = false(num_samples, 1);
for k = 1:length(segments)
    is_active_mask(segments(k).start_idx : segments(k).end_idx) = true;
end

% 缓存卫星投影 (每秒更新一次)
last_update = -999;
curr_sat_projs = []; % N x 2
curr_s_indices = [];

h_bar = waitbar(0, '物理引擎计算中...');

for t = 1 : num_samples
    if mod(t, 100) == 0, waitbar(t/num_samples, h_bar); end
    if ~is_active_mask(t), continue; end
    
    % --- A. 更新卫星投影几何 (Geometry) ---
    curr_time_val = posixtime(t_grid(t));
    if curr_time_val - last_update > 0.5
        [~, epoch_idx] = min(abs([obs_clean.time] - t_grid(t)));
        try [~, ~, sat_states] = calculate_receiver_position(obs_clean, nav_data, epoch_idx); catch, continue; end
        
        curr_sat_projs = [];
        curr_s_indices = [];
        
        for s = 1:length(valid_sats)
            sid = valid_sats{s};
            if isfield(sat_states, sid)
                sat_p = sat_states.(sid).position;
                [e, n, u] = ecef2enu(sat_p(1)-rec_xyz(1), sat_p(2)-rec_xyz(2), sat_p(3)-rec_xyz(3), lat0, lon0, alt0);
                
                if u > 0
                    % 投影公式: P_proj = (H / u) * [e, n]
                    scale = MODEL.gesture_height / u;
                    curr_sat_projs = [curr_sat_projs; scale*e, scale*n];
                    curr_s_indices = [curr_s_indices; s];
                end
            end
        end
        last_update = curr_time_val;
    end
    
    if isempty(curr_sat_projs), continue; end
    
    % --- B. 获取观测向量 (Observation) ---
    % obs_vec: 当前时刻各卫星的能量跌落值
    obs_vec = vol_mat(t, curr_s_indices)'; % N x 1
    
    % 二值化观测状态: 1=遮挡, -1=未遮挡, 0=不确定
    obs_state = zeros(size(obs_vec));
    obs_state(obs_vec > PARA.hit_threshold) = 1;
    obs_state(obs_vec < 0.5) = -1;
    
    % 如果全都没遮挡，跳过
    if ~any(obs_state == 1), continue; end
    
    % --- C. 矩阵化网格匹配 (Vectorized Matching) ---
    % 这是一个 N x M 的大矩阵运算
    % N = 卫星数, M = 网格点数
    
    % Grid (M x 2), Body (1 x 2)
    % 虚拟手臂线段: A=Grid, B=Body. 向量 AB = B - A
    A = grid_pts;
    B = MODEL.body_pos;
    AB_x = B(1) - A(:,1); % M x 1
    AB_y = B(2) - A(:,2);
    len_sq = AB_x.^2 + AB_y.^2;
    
    scores = zeros(num_grid, 1);
    
    for k = 1:length(obs_state)
        val = obs_state(k);
        if val == 0, continue; end % 不贡献分数
        
        P = curr_sat_projs(k, :); % 当前卫星投影 1 x 2
        
        % 计算 P 到所有 M 个线段的距离
        % AP = P - A
        AP_x = P(1) - A(:,1);
        AP_y = P(2) - A(:,2);
        
        % t = dot(AP, AB) / len_sq
        dot_v = AP_x .* AB_x + AP_y .* AB_y;
        t_proj = max(0, min(1, dot_v ./ (len_sq + 1e-6)));
        
        % 最近点 C
        C_x = A(:,1) + t_proj .* AB_x;
        C_y = A(:,2) + t_proj .* AB_y;
        
        dist_sq = (P(1) - C_x).^2 + (P(2) - C_y).^2;
        
        % 理论是否遮挡 (M x 1 logical)
        is_blocked_theory = dist_sq < (MODEL.arm_width^2);
        
        % 匹配得分: 
        % 实际=1, 理论=1 -> +1
        % 实际=-1, 理论=0 -> +1
        % 否则 -> -1
        match = (val .* (2*double(is_blocked_theory) - 1));
        scores = scores + match;
    end
    
    % --- D. 取最优解 ---
    [max_score, best_idx] = max(scores);
    
    % 只有当得分足够高（说明找到了合理解释）才采纳
    if max_score > 0
        raw_x(t) = grid_pts(best_idx, 1);
        raw_y(t) = grid_pts(best_idx, 2);
        raw_conf(t) = max_score;
    end
end
close(h_bar);

% 简单填补 NaN
raw_x = fillmissing(raw_x, 'linear');
raw_y = fillmissing(raw_y, 'linear');

%% ================= [Part 2] 几何笔画精修 (Stroke Refinement) =================
fprintf('    2. 执行几何笔画精修 (消除网格抖动 & 拉直)...\n');

% 1. 平滑数据准备切分
sx = smoothdata(raw_x, 'gaussian', 10);
sy = smoothdata(raw_y, 'gaussian', 10);

% 2. 基于角速度的折点检测
vx = gradient(sx); vy = gradient(sy);
vel = sqrt(vx.^2 + vy.^2);
ang = atan2(vy, vx);
ang_unwrap = unwrap(ang);
ang_vel = [0; abs(diff(ang_unwrap))];
% 仅在运动时检测角度变化
ang_vel(vel < max(vel)*0.1) = 0; 

% 寻找峰值 (折点)
[~, locs] = findpeaks(ang_vel, 'MinPeakHeight', 0.5, 'MinPeakDistance', 10);
cut_idx = unique([1; locs; length(t_grid)]);

% 3. PCA 分段拟合
final_traj.x = nan(size(raw_x));
final_traj.y = nan(size(raw_y));
final_traj.t = t_grid;
final_traj.strokes = [];

fprintf('       检测到 %d 个笔画，正在线性化...\n', length(cut_idx)-1);

for k = 1:length(cut_idx)-1
    idx_s = cut_idx(k);
    idx_e = cut_idx(k+1);
    rng = idx_s:idx_e;
    
    if ~any(is_active_mask(rng)), continue; end
    if length(rng) < 5, continue; end % 太短忽略
    
    % 提取该段
    seg_x = raw_x(rng);
    seg_y = raw_y(rng);
    
    % PCA 直线拟合
    pts = [seg_x, seg_y];
    mu = mean(pts, 1);
    [coeff, ~, ~] = svd(cov(pts - mu));
    dir_vec = coeff(:, 1); % 主方向
    
    % 投影
    pts_proj = (pts - mu) * dir_vec * dir_vec' + mu;
    
    % 存回
    final_traj.x(rng) = pts_proj(:, 1);
    final_traj.y(rng) = pts_proj(:, 2);
    
    % 记录笔画信息
    stroke.start = pts_proj(1,:);
    stroke.end   = pts_proj(end,:);
    stroke.vec   = dir_vec;
    final_traj.strokes = [final_traj.strokes; stroke];
end

%% ================= [Part 3] 绘图验证 =================
figure('Name', 'Physics Engine + Stroke Refinement', 'Color', 'w');
subplot(1,2,1);
hold on; grid on; axis equal;
plot(raw_x, raw_y, '.', 'Color', [0.7 0.7 0.7]);
title('Step 1: 物理逆向求解 (Raw)');
plot(MODEL.body_pos(1), MODEL.body_pos(2), 'bs', 'MarkerFaceColor', 'b');
xlabel('East'); ylabel('North');

subplot(1,2,2);
hold on; grid on; axis equal;
% 画 Body
plot(MODEL.body_pos(1), MODEL.body_pos(2), 'bs', 'MarkerFaceColor', 'b', 'DisplayName', 'Body');
% 画接收机
plot(0, 0, '^k', 'MarkerFaceColor', 'k', 'DisplayName', 'Receiver');
% 画最终轨迹
plot(final_traj.x, final_traj.y, 'r.-', 'LineWidth', 2, 'DisplayName', 'Final Trajectory');
% 画起点
valid_idx = find(~isnan(final_traj.x));
if ~isempty(valid_idx)
    plot(final_traj.x(valid_idx(1)), final_traj.y(valid_idx(1)), 'go', 'MarkerFaceColor', 'g', 'DisplayName', 'Start');
end

title('Step 2: 几何约束精修 (Final)');
legend;
xlabel('East'); ylabel('North');

fprintf('✅ 物理引擎分析完成。\n');

end