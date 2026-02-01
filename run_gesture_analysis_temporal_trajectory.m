% =========================================================================
% run_gesture_analysis_linear_boundary_scan.m
% 功能: 基于“方向扫描”与“边界约束”的线性轨迹恢复
%
% 核心思想 (User Inspired):
%   1. [直线假设]: 认定短时手势为直线运动，不进行曲线/样条拟合。
%   2. [时序方向]: 利用 Hit 点的时间顺序计算主方向向量 (Primary Direction)。
%   3. [边界扫描]: 寻找一条直线，使其：
%      - 最大程度覆盖 Hit 点 (Attraction)
%      - 严格避开 Miss 点 (Boundary Repulsion)
%      - 类似于在 Miss 点构成的“墙壁”中寻找一条走廊。
%
% 输入:
%   - obs_data: 观测数据 (align_duration 处理后的 0/10 数据)
%   - nav_data: 导航电文
%
% 输出:
%   - analysis_results: 包含直线起终点、方向角等信息
% =========================================================================

function analysis_results = run_gesture_analysis_linear_boundary_scan(obs_data, nav_data)

    addpath(genpath('calculate_clock_bias_and_positon'));
    addpath(genpath('sky_plot'));

    fprintf('--> [Linear Scan] 启动线性边界扫描算法 (Linear Boundary Scan)...\n');

    %% 1. 参数配置
    PARA.gesture_height    = 0.30;   % 平面高度
    PARA.miss_radius       = 0.05;   % 障碍物半径
    PARA.min_elevation     = 15;     
    PARA.time_gap_split    = 1.0;    % 动作切分时间阈值 (秒)
    
    % 扫描参数
    PARA.scan_angle_range  = 30;     % 在主方向基础上左右扫描 +/- 30度
    PARA.scan_step         = 1.0;    % 扫描步长 1度

    %% 2. 预处理与全局 Hit/Miss 提取
    raw_times = [obs_data.time];
    if isempty(raw_times), analysis_results = []; return; end
    t0_global = raw_times(1);
    
    % 2.1 计算接收机位置 (基准)
    mid_idx = round(length(obs_data)/2);
    try
        [rec_pos, ~, sat_states_all] = calculate_receiver_position(obs_data, nav_data, mid_idx);
        % 若失败尝试重试...
        if any(isnan(rec_pos))
             for r = 100:100:length(obs_data)
                 [rec_pos, ~, sat_states_all] = calculate_receiver_position(obs_data, nav_data, r);
                 if ~any(isnan(rec_pos)), break; end
             end
        end
        [lat0, lon0, alt0] = ecef2geodetic(rec_pos(1), rec_pos(2), rec_pos(3));
    catch
        error('无法计算接收机位置');
    end

    % 2.2 全局扫描 Hit/Miss 点
    fprintf('    正在提取 Hit/Miss 特征点...\n');
    [all_hits, all_miss] = extract_features(obs_data, rec_pos, lat0, lon0, alt0, sat_states_all, PARA, t0_global, raw_times);
    
    if isempty(all_hits)
        fprintf('    [警告] 未检测到 Hit 点。\n');
        analysis_results = [];
        return;
    end

    %% 3. 动作聚类 (Time Clustering)
    % 按时间将 Hit 点切分为独立的笔画
    [~, sort_idx] = sort([all_hits.t_rel]);
    sorted_hits = all_hits(sort_idx);
    
    dt = diff([sorted_hits.t_rel]);
    split_indices = [0, find(dt > PARA.time_gap_split), length(sorted_hits)];
    
    analysis_results = struct('id', {}, 'p_start', {}, 'p_end', {}, 'direction', {});
    
    %% 4. 逐个动作进行线性扫描
    for i = 1:length(split_indices)-1
        curr_idx = (split_indices(i)+1) : split_indices(i+1);
        cluster_hits = sorted_hits(curr_idx);
        
        if length(cluster_hits) < 3, continue; end % 点太少
        
        seg_id = i;
        fprintf('\n--- 拟合动作 #%d (Hits: %d) ---\n', seg_id, length(cluster_hits));
        
        % 4.1 确定初始主方向 (Primary Direction)
        % 取前20%点的重心作为 Start，后20%点的重心作为 End
        n_pts = length(cluster_hits);
        n_head = max(1, round(n_pts * 0.3));
        
        pts_head = vertcat(cluster_hits(1:n_head).pos);
        pts_tail = vertcat(cluster_hits(end-n_head+1:end).pos);
        
        p_start_init = mean(pts_head, 1);
        p_end_init   = mean(pts_tail, 1);
        
        vec_init = p_end_init - p_start_init;
        len_init = norm(vec_init);
        if len_init < 0.01, vec_init = [1, 0]; end % 防止原地不动
        dir_init = vec_init / norm(vec_init);
        base_angle = atan2d(dir_init(2), dir_init(1));
        
        fprintf('    初始主方向: %.1f 度\n', base_angle);
        
        % 4.2 边界扫描与优化 (Boundary Scan)
        % 目标：找到一个角度 theta 和 偏移 offset，使得直线穿过 Hit 且避开 Miss
        
        best_score = -inf;
        best_line = struct('p1', p_start_init, 'p2', p_end_init);
        
        % 提取 Miss 障碍物 (只考虑附近的)
        cluster_center = mean(vertcat(cluster_hits.pos), 1);
        relevant_miss = [];
        if ~isempty(all_miss)
            m_pos = vertcat(all_miss.pos);
            dists = sqrt(sum((m_pos - cluster_center).^2, 2));
            relevant_miss = m_pos(dists < 1.5, :); % 只看1.5m内的障碍
        end
        
        % === 扫描循环 ===
        angles = linspace(base_angle - PARA.scan_angle_range, base_angle + PARA.scan_angle_range, 30);
        
        for theta = angles
            dir_scan = [cosd(theta), sind(theta)];
            
            % 将所有 Hit 点投影到该方向的法向量上，寻找最佳中心线
            % 法向量 normal = [-sin, cos]
            normal = [-dir_scan(2), dir_scan(1)];
            
            % 计算 Hit 点在法向量上的投影距离 (即横向偏差)
            h_pos = vertcat(cluster_hits.pos);
            offsets = h_pos * normal'; 
            
            % 直线的中心位置 (横向偏移均值)
            center_offset = median(offsets); % 使用中位数抗噪
            
            % 构造直线上的投影范围 (纵向)
            projs = h_pos * dir_scan';
            l_min = min(projs);
            l_max = max(projs);
            
            % 构造候选直线段
            p1_cand = (l_min * dir_scan) + (center_offset * normal);
            p2_cand = (l_max * dir_scan) + (center_offset * normal);
            
            % === 评分机制 ===
            score = 0;
            
            % 1. Hit 覆盖度 (点到直线的距离)
            dists_h = abs(offsets - center_offset);
            score = score - sum(dists_h.^2); % 距离平方和越小越好
            
            % 2. Miss 碰撞惩罚 (Boundary Check)
            if ~isempty(relevant_miss)
                % 计算 Miss 到直线的距离
                m_offsets = relevant_miss * normal';
                m_dists = abs(m_offsets - center_offset);
                
                % 投影范围检查 (Miss 是否在 线段长度范围内)
                m_projs = relevant_miss * dir_scan';
                in_segment = (m_projs > l_min - 0.1) & (m_projs < l_max + 0.1);
                
                % 碰撞检测
                collision = (m_dists < PARA.miss_radius) & in_segment;
                if any(collision)
                    score = score - 1000; % 严重惩罚
                end
                
                % 软约束：离障碍物越远越好
                min_m_dist = min(m_dists(in_segment));
                if ~isempty(min_m_dist)
                    score = score + min_m_dist * 10; 
                end
            end
            
            if score > best_score
                best_score = score;
                best_line.p1 = p1_cand;
                best_line.p2 = p2_cand;
            end
        end
        
        % 4.3 存储结果
        cnt = length(analysis_results) + 1;
        analysis_results(cnt).id = seg_id;
        analysis_results(cnt).p_start = best_line.p1;
        analysis_results(cnt).p_end   = best_line.p2;
        analysis_results(cnt).hit_points = cluster_hits;
        
        % 4.4 单次绘图
        visualize_scan_result(seg_id, cluster_hits, relevant_miss, best_line, PARA);
    end
    
    %% 5. 全局汇总绘图
    visualize_all_results(analysis_results, all_miss);
    fprintf('✅ 线性边界扫描分析完成。\n');
end

%% ================== 辅助函数 ==================

function [hits, misses] = extract_features(obs, rec_pos, lat0, lon0, alt0, sat_states, para, t0, times)
    hits = struct('sat', {}, 'pos', {}, 't_rel', {}, 'val', {});
    misses = struct('sat', {}, 'pos', {});
    
    % 获取所有可见卫星
    all_sats = {};
    for k=1:min(50, length(obs)), if ~isempty(obs(k).data), all_sats=[all_sats, fieldnames(obs(k).data)']; end; end
    unique_sats = unique(all_sats);
    
    for s_idx = 1:length(unique_sats)
        sat_id = unique_sats{s_idx};
        sys = sat_id(1);
        if ~ismember(sys, ['G','C','E','J','R']), continue; end
        
        cn0_seq = extract_processed_snr(obs, sat_id);
        if all(isnan(cn0_seq)), continue; end
        
        if isfield(sat_states, sat_id)
            sat_pos = sat_states.(sat_id).position;
        else
            continue;
        end
        
        [e, n, u] = ecef2enu(sat_pos(1)-rec_pos(1), sat_pos(2)-rec_pos(2), sat_pos(3)-rec_pos(3), lat0, lon0, alt0);
        vec = [e, n, u]; vec = vec/norm(vec);
        
        if asind(vec(3)) < para.min_elevation || vec(3) <= 0, continue; end
        pt = (para.gesture_height / vec(3)) * vec(1:2);
        if norm(pt) > 2.5, continue; end
        
        max_val = max(cn0_seq, [], 'omitnan');
        
        % 信号分类逻辑 (严谨版)
        if max_val > 18.0
            % Miss (High SNR)
            m_cnt = length(misses)+1;
            misses(m_cnt).sat = sat_id;
            misses(m_cnt).pos = pt;
        elseif max_val > 8.0 && max_val < 15.0
            % Hit (Processed 10.0)
            idx_hits = find(cn0_seq > 5.0);
            if ~isempty(idx_hits)
                diff_idx = [2; diff(idx_hits)];
                starts = idx_hits(diff_idx > 1);
                for k=1:length(starts)
                    h_cnt = length(hits)+1;
                    hits(h_cnt).sat = sat_id;
                    hits(h_cnt).pos = pt;
                    hits(h_cnt).t_rel = seconds(times(starts(k)) - t0);
                    hits(h_cnt).val = max_val;
                end
            end
        end
    end
end

function cn0 = extract_processed_snr(obs, sat)
    % 强制读取第一个字段 (对齐 align_duration)
    cn0 = NaN(length(obs), 1);
    field_name = '';
    for k=1:min(10, length(obs))
        if isfield(obs(k).data, sat) && isfield(obs(k).data.(sat), 'snr')
            f = fieldnames(obs(k).data.(sat).snr);
            if ~isempty(f), field_name = f{1}; break; end
        end
    end
    if isempty(field_name), return; end
    for k=1:length(obs)
        if isfield(obs(k).data, sat) && isfield(obs(k).data.(sat).snr, field_name)
            cn0(k) = obs(k).data.(sat).snr.(field_name);
        end
    end
end

function visualize_scan_result(id, hits, misses, line, para)
    figure('Name', sprintf('Gesture #%d Linear Scan', id), 'Color', 'w', 'Position', [100,100,600,600]);
    hold on; axis equal; grid on;
    xlabel('East (m)'); ylabel('North (m)');
    
    % Draw Miss
    if ~isempty(misses)
        plot(misses(:,1), misses(:,2), 'x', 'Color', [0.6 0.6 0.6]);
        theta = 0:0.1:2*pi;
        for i=1:size(misses,1)
            cx = misses(i,1) + para.miss_radius*cos(theta);
            cy = misses(i,2) + para.miss_radius*sin(theta);
            plot(cx, cy, '-', 'Color', [0.8 0.8 0.8]);
        end
    end
    
    % Draw Hits
    h_pos = vertcat(hits.pos);
    plot(h_pos(:,1), h_pos(:,2), 'bo', 'MarkerFaceColor', 'b');
    
    % Draw Best Line
    plot([line.p1(1), line.p2(1)], [line.p1(2), line.p2(2)], 'r-', 'LineWidth', 3);
    plot(line.p1(1), line.p1(2), 'g^', 'MarkerFaceColor', 'g', 'MarkerSize', 10);
    plot(line.p2(1), line.p2(2), 'rv', 'MarkerFaceColor', 'r', 'MarkerSize', 10);
    
    title(sprintf('Gesture #%d: Linear Boundary Scan', id));
end

function visualize_all_results(res, misses)
    if isempty(res), return; end
    figure('Name', 'Global Linear Trajectories', 'Color', 'w');
    hold on; axis equal; grid on;
    xlabel('East (m)'); ylabel('North (m)');
    
    if ~isempty(misses)
        mp = vertcat(misses.pos);
        plot(mp(:,1), mp(:,2), '.', 'Color', [0.8 0.8 0.8]);
    end
    
    cols = lines(length(res));
    for i=1:length(res)
        p1 = res(i).p_start;
        p2 = res(i).p_end;
        
        plot([p1(1), p2(1)], [p1(2), p2(2)], '-', 'Color', cols(i,:), 'LineWidth', 2);
        
        % Arrow head
        quiver(p1(1), p1(2), p2(1)-p1(1), p2(2)-p1(2), 0, 'Color', cols(i,:), 'MaxHeadSize', 0.5);
        
        mid = (p1+p2)/2;
        text(mid(1), mid(2), sprintf('#%d', res(i).id), 'BackgroundColor', 'w', 'EdgeColor', cols(i,:));
    end
    title('Final Recovered Linear Trajectories');
end