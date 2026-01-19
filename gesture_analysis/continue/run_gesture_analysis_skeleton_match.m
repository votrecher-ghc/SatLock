% =========================================================================
% run_gesture_analysis_skeleton_match.m (新函数)
% 功能: 鲁棒手势感知 Step 2 - 骨架匹配/正向建模版 (Skeleton Match)
% 描述:
%   该函数实现了一种基于"正向建模与假设验证" (Forward Modeling & Hypothesis Verification)
%   的轨迹追踪算法。与传统的从数据推位置(逆向)不同，本算法利用先验知识（人体骨架结构）
%   来"猜"指尖位置，并根据理论遮挡与实际观测的吻合度进行打分。
%
%   核心逻辑:
%   1. 建立虚拟肩膀 (Virtual Shoulder) 作为手臂锚点。
%   2. 在手势平面生成密集的候选指尖网格 (Candidate Grid)。
%   3. 对每个候选点构建"虚拟手臂"，计算其理论上会遮挡哪些卫星。
%   4. 将理论遮挡与实际方波观测进行比对：
%      - 命中 (Hit): 理论遮挡且实际有波动 -> 加分
%      - 冲突 (Conflict): 理论遮挡但实际无波动 -> 强力扣分 (排除法)
%   5. 选取每一时刻得分最高的候选点作为最优估计。
%
% [调用格式]:
%   [traj_x, traj_y, traj_t, traj_max_scores] = run_gesture_analysis_skeleton_match(obs_waveform, nav_data, step1_res_shaped);
%
% [输入参数]:
%   1. obs_waveform (struct): 
%      经过 waveform_reshaping 处理后的方波观测数据 (0 或 10)。
%   2. nav_data (struct): 
%      导航星历数据。
%   3. step1_res_shaped (struct): 
%      Step 1 结果包 (含整形后的 volatility_matrix)。
%
% [返回值说明]:
%   1. traj_x / traj_y (double列向量): 最优匹配的指尖轨迹坐标 (米)。
%   2. traj_t (double列向量): 轨迹点时间索引。
%   3. traj_max_scores (double列向量): 轨迹点的匹配得分 (置信度)。
% =========================================================================

function [traj_x, traj_y, traj_t, traj_max_scores] = run_gesture_analysis_skeleton_match(obs_waveform, nav_data, step1_res_shaped)

% --- 1. 数据解包 (Unpack Data) ---
segments = step1_res_shaped.segments;
volatility_matrix = step1_res_shaped.volatility_matrix;
t_grid = step1_res_shaped.t_grid;
valid_sats = step1_res_shaped.valid_sats;
% PARA = step1_res_shaped.PARA;

if isempty(segments)
    fprintf('⚠️ 警告: step1_res 中无有效分段，Step 2 (Skeleton Match) 终止。\n');
    traj_x=[]; traj_y=[]; traj_t=[]; traj_max_scores=[];
    return;
end

addpath(genpath('sky_plot')); 
addpath(genpath('calculate_clock_bias_and_positon'));
addpath(genpath('nav_parse'));

fprintf('--> [Step 2] 启动骨架匹配追踪 (Skeleton Match / Forward Modeling)...\n');

%% ================= [Part 1] 参数设置 =================

% 1. 骨架模型参数 (Skeleton Model)
%    假设坐标系以接收机为原点 (0,0,0)。
%    肩膀位于身后、下方。
MODEL.shoulder_pos     = [0.0, -0.50, -0.30]; % [x, y, z] 米
MODEL.arm_radius       = 0.10;               % [物理] 手臂/遮挡圆柱体的有效半径 (米)

% 2. 网格搜索参数 (Grid Search)
%    手势发生的活动平面区域
GRID.x_range           = [-0.5, 0.5];        % X轴搜索范围 (米)
GRID.y_range           = [-0.2, 0.6];        % Y轴搜索范围 (米)
GRID.z_height          = 0.20;               % Z轴固定高度 (米)
GRID.step_size         = 0.05;               % 网格步长 (米)，越小越精细但越慢

% 3. 打分权重 (Scoring Weights) [核心]
%    Score = (Hit * W_HIT) + (Conflict * W_CONFLICT)
SCORE.w_hit            = 1.0;   % [奖励] 理论与实际均遮挡 (吻合)
SCORE.w_conflict       = -3.0;  % [惩罚] 理论遮挡但实际未遮挡 (严重不符)
SCORE.w_miss           = 0.0;   % [中性] 理论未遮挡 (不提供信息)

% 4. 仰角安全过滤
ALG.zenith_safe_deg    = 15;    % [掩膜] 最低安全仰角 (度)

% 5. 轨迹平滑
TRAJ.traj_smooth_m     = 3;     % [平滑] 窗口大小

%% ================= [Part 2] 预计算网格 =================

fprintf('    生成候选网格...\n');
[gx, gy] = meshgrid(GRID.x_range(1):GRID.step_size:GRID.x_range(2), ...
                    GRID.y_range(1):GRID.step_size:GRID.y_range(2));
candidate_pts = [gx(:), gy(:), ones(numel(gx), 1) * GRID.z_height];
num_candidates = size(candidate_pts, 1);
fprintf('    -> 候选点数量: %d\n', num_candidates);

%% ================= [Part 3] 核心计算流程 =================

% 1. 构建激活掩膜
num_samples = length(t_grid);
is_active_mask = false(num_samples, 1);
for k = 1:length(segments)
    is_active_mask(segments(k).start_idx : segments(k).end_idx) = true;
end

% 2. 预计算接收机位置缓存
rec_pos_cache = NaN(num_samples, 3);
for t = 1:num_samples
    if is_active_mask(t)
        [~, epoch_idx] = min(abs([obs_waveform.time] - t_grid(t)));
        try [rp, ~, ~] = calculate_receiver_position(obs_waveform, nav_data, epoch_idx); rec_pos_cache(t,:) = rp; catch, end
    end
end

% 3. 主循环: 逐时刻匹配
track_results = struct('t_idx', {}, 'x', {}, 'y', {}, 'score', {});
track_cnt = 0;
is_tracking_started = false; 

fprintf('--> 执行骨架匹配 (Grid=%dx%d, ConflictWeight=%.1f)...\n', size(gx,1), size(gx,2), SCORE.w_conflict);

% 为了提高速度，这里不进行聚类跳步，或者步长设为1
% 如果太慢可以设 step = 2 或 3
time_step = 1; 

for t = 1 : time_step : num_samples
    
    if ~is_active_mask(t), continue; end
    
    % --- A. 准备卫星向量 ---
    rec_pos = rec_pos_cache(t, :); 
    if any(isnan(rec_pos)), continue; end
    [lat0, lon0, alt0] = ecef2geodetic(rec_pos(1), rec_pos(2), rec_pos(3));
    
    current_vols = volatility_matrix(t, :);
    
    [~, epoch_idx] = min(abs([obs_waveform.time] - t_grid(t)));
    try [~, ~, sat_states] = calculate_receiver_position(obs_waveform, nav_data, epoch_idx); catch, continue; end
    
    % 收集当前可见卫星的视线向量 (Unit Vectors) 和 观测状态 (Status)
    sat_vecs = [];      % [N x 3] 视线单位向量 (ENU)
    sat_status = [];    % [N x 1] 1=被遮挡(Vol=10), 0=未遮挡(Vol=0)
    
    for s = 1:length(valid_sats)
        sid = valid_sats{s};
        if ~isfield(sat_states, sid), continue; end
        
        sat_p = sat_states.(sid).position;
        [e, n, u] = ecef2enu(sat_p(1)-rec_pos(1), sat_p(2)-rec_pos(2), sat_p(3)-rec_pos(3), lat0, lon0, alt0);
        dist = norm([e, n, u]);
        vec_u = [e, n, u] / dist; 
        
        % 仰角过滤
        el_deg = asind(vec_u(3));
        if el_deg < ALG.zenith_safe_deg, continue; end
        
        % 状态判定 (0 或 10)
        is_blocked = (current_vols(s) > 5.0); % 阈值取中间值 5
        
        sat_vecs = [sat_vecs; vec_u]; 
        sat_status = [sat_status; is_blocked];
    end
    
    if isempty(sat_vecs), continue; end
    
    % --- B. 骨架匹配 (Batch Scoring) ---
    % 我们需要计算所有 Candidate 的得分，取最大值
    
    best_score = -inf;
    best_idx = -1;
    
    shoulder = MODEL.shoulder_pos; % [1x3]
    
    % 遍历每一个候选指尖 (Candidate Fingertip)
    for c = 1:num_candidates
        fingertip = candidate_pts(c, :); % [1x3]
        
        % 虚拟手臂向量: Shoulder -> Fingertip
        arm_vec = fingertip - shoulder; 
        arm_len_sq = dot(arm_vec, arm_vec);
        
        if arm_len_sq < 1e-6, continue; end % 防止除零
        
        current_score = 0;
        
        % 遍历每颗卫星，计算几何关系
        % 简化模型: 计算线段(Shoulder-Fingertip)到射线(Origin-SatVec)的距离
        % 由于射线过原点(0,0,0)，问题转化为:
        % 线段上的某一点 P(s) = Shoulder + s * ArmVec (0<=s<=1)
        % 是否非常靠近 射线 Line(k) = k * SatVec
        % 
        % 距离公式: d(s) = || P(s) - (P(s) dot SatVec) * SatVec ||
        % 我们需要找 s in [0,1] 使得 d(s) 最小
        
        for k = 1:size(sat_vecs, 1)
            u = sat_vecs(k, :); % Unit vector of LoS
            
            % P(s) = S + s * V
            % Project P(s) onto line u: Proj(s) = (S + sV) dot u * u
            % DistSq(s) = || (S + sV) - ((S + sV) dot u) * u ||^2
            % 这是一个关于 s 的二次函数。
            % 令 vector A = S - (S dot u)u
            % 令 vector B = V - (V dot u)u
            % DistSq(s) = || A + sB ||^2, 需在 s=[0,1] 内最小化
            
            dot_su = dot(shoulder, u);
            dot_vu = dot(arm_vec, u);
            
            A = shoulder - dot_su * u;
            B = arm_vec  - dot_vu * u;
            
            dot_BB = dot(B, B);
            dot_AB = dot(A, B);
            
            % 极值点 s_star (无约束)
            if dot_BB < 1e-6
                s_star = 0; % B is zero vector, distance is constant
            else
                s_star = -dot_AB / dot_BB;
            end
            
            % 截断到 [0, 1] 线段范围内
            s_clamped = max(0, min(1, s_star));
            
            % 计算最小距离
            vec_dist = A + s_clamped * B;
            min_dist = norm(vec_dist);
            
            % --- 判定与打分 ---
            is_intersect = (min_dist < MODEL.arm_radius);
            is_obs_blocked = sat_status(k);
            
            if is_intersect && is_obs_blocked
                % [Hit] 理论挡，实际挡 -> 加分
                current_score = current_score + SCORE.w_hit;
            elseif is_intersect && ~is_obs_blocked
                % [Conflict] 理论挡，实际没挡 -> 扣分
                current_score = current_score + SCORE.w_conflict;
            elseif ~is_intersect && is_obs_blocked
                % [Miss] 理论没挡，实际挡了 (可能是噪音，或者被身体其他部位挡了)
                % 暂时忽略，或者给个微小惩罚
                current_score = current_score + SCORE.w_miss;
            end
        end
        
        if current_score > best_score
            best_score = current_score;
            best_idx = c;
        end
    end
    
    % --- C. 记录结果 ---
    if best_idx ~= -1 && best_score > 0 % 至少要有正向匹配
        final_pt = candidate_pts(best_idx, :);
        
        % 动作触发锁
        if ~is_tracking_started
            % 只有当高分点离身体有一定距离时才开始
            if norm(final_pt(1:2)) > 0.05
                is_tracking_started = true;
            else
                continue;
            end
        end
        
        track_cnt = track_cnt + 1;
        track_results(track_cnt).t_idx = t;
        track_results(track_cnt).x = final_pt(1);
        track_results(track_cnt).y = final_pt(2);
        track_results(track_cnt).score = best_score;
    end
end

% 4. 轨迹提取与输出
if track_cnt > 0
    traj_x = [track_results.x]'; 
    traj_y = [track_results.y]'; 
    traj_t = [track_results.t_idx]'; 
    traj_max_scores = [track_results.score]';
    
    % 平滑处理
    if TRAJ.traj_smooth_m > 1
        traj_x = smoothdata(traj_x, 'movmean', TRAJ.traj_smooth_m);
        traj_y = smoothdata(traj_y, 'movmean', TRAJ.traj_smooth_m);
    end
else
    traj_x = []; traj_y = []; traj_t = []; traj_max_scores = [];
end


%% ================= [Part 4] 绘图流程 =================
fprintf('\n--> 开始生成图表...\n');

if ~isempty(traj_x)
    figure('Name', 'Skeleton Match Trajectory (Forward Model)', 'Position', [100, 200, 600, 600], 'Color', 'w');
    ax = axes; hold(ax, 'on'); grid(ax, 'on'); axis(ax, 'equal');
    xlabel('East (m)'); ylabel('North (m)');
    
    % 画接收机
    plot(ax, 0, 0, '^', 'MarkerSize', 12, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'DisplayName', 'Receiver');
    viscircles(ax, [0,0], 0.05, 'Color', [0.8 0.8 0.8], 'LineStyle', '--'); 
    
    % 画肩膀位置示意
    plot(ax, MODEL.shoulder_pos(1), MODEL.shoulder_pos(2), 'ms', 'MarkerFaceColor', 'm', 'DisplayName', 'Virtual Shoulder');
    
    % 轨迹线
    plot(ax, traj_x, traj_y, 'b.-', 'LineWidth', 1.5, 'MarkerSize', 8, 'DisplayName', 'Matched Fingertip');
    
    % 颜色映射得分
    % scatter(ax, traj_x, traj_y, 30, traj_max_scores, 'filled');
    % c = colorbar; c.Label.String = 'Match Score'; colormap(ax, 'parula');
    
    % 起终点
    plot(ax, traj_x(1), traj_y(1), 'go', 'MarkerFaceColor', 'g', 'MarkerSize', 8, 'DisplayName', 'Start');
    plot(ax, traj_x(end), traj_y(end), 'rs', 'MarkerFaceColor', 'r', 'MarkerSize', 8, 'DisplayName', 'End');
    
    title({'Skeleton Match Track', ...
           sprintf('Conflict W: %.1f | Arm Radius: %.2fm', SCORE.w_conflict, MODEL.arm_radius)});
    legend('Location', 'best');
    
    max_range = max([max(abs(traj_x)), max(abs(traj_y)), 0.5]);
    xlim([-max_range*1.2, max_range*1.2]); ylim([-max_range*1.2, max_range*1.2]);
else
    fprintf('⚠️ 本次未生成有效轨迹。\n');
end

fprintf('✅ 骨架匹配分析完成。\n');
end