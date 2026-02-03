% % =========================================================================
% % run_gesture_analysis_boundary_trackV3.m (v6.0 动态稀疏度加权版)
% % 功能: 鲁棒手势感知 Step 2 - 模拟均匀分布 / 动态稀疏度加权
% % 描述:
% %   [核心假设]: 之前的重心偏移和畸变是由卫星分布不均(GDOP)造成的。
% %   [解决方案]: 假设理想情况下卫星应均匀分布。通过计算卫星在天空视图(Skyplot)
% %              上的局部密度，给予"稀疏区域"卫星更高的权重，给予"密集区域"
% %              卫星更低的权重，从而在算法层面"抹平"分布差异。
% %
% %   [算法步骤]:
% %     1. 计算瞬时卫星位置 (Azimuth, Elevation)。
% %     2. Sparsity Weighting: 对每一颗卫星，计算其与所有其他卫星的角距离。
% %        密度 = Sum(exp(-dist^2 / sigma^2))。
% %        稀疏权重 = 1 / 密度。
% %     3. Weighted Centroid: 使用 (信号能量 * 稀疏权重) 计算轨迹重心。
% % =========================================================================
% 
% function [traj_x, traj_y, traj_t, traj_e] = run_gesture_analysis_boundary_trackV3(obs_clean, nav_data, step1_res)
% 
% % --- 1. 数据解包 ---
% segments = step1_res.segments;
% volatility_matrix = step1_res.volatility_matrix;
% t_grid = step1_res.t_grid;
% valid_sats = step1_res.valid_sats;
% PARA = step1_res.PARA; 
% 
% if isempty(segments)
%     fprintf('⚠️ 警告: step1_res 中无有效分段，Step 2 终止。\n');
%     traj_x=[]; traj_y=[]; traj_t=[]; traj_e=[];
%     return;
% end
% 
% addpath(genpath('sky_plot')); 
% addpath(genpath('calculate_clock_bias_and_positon'));
% addpath(genpath('nav_parse'));
% 
% fprintf('--> [Step 2] 启动 v6.0 动态稀疏度加权 (Dynamic Sparsity Weighting)...\n');
% 
% %% ================= [Part 1] 参数设置 =================
% 
% % 1. 轨迹与几何参数
% TRAJ.gesture_height    = 0.30; 
% TRAJ.min_elevation     = 10;   
% TRAJ.min_action_dist   = 0.05; 
% 
% % 2. 身体参考参数
% ALG.body_pos           = [0.0, -1.0]; 
% 
% % 3. [核心] 稀疏度加权参数
% ALG.density_sigma      = 30;    % [度] 高斯核带宽。决定"多近才算拥挤"。
%                                 % 30度意味着两个卫星相距30度以内会显著互相削弱权重。
% ALG.sparsity_power     = 2.0;   % [强度] 稀疏权重的指数。
%                                 % 1.0 = 线性补偿; 2.0 = 强力补偿 (越稀疏越重要)
% 
% % 4. 原始点提取 (Top 15%)
% ALG.top_k_percent      = 0.15; 
% ALG.min_sat_vol        = 1.0; 
% 
% % 5. 平滑参数
% TRAJ.time_cluster_k    = 5;    
% TRAJ.traj_smooth_m     = 3;    
% 
% %% ================= [Part 2] 核心计算流程 =================
% 
% % 1. 构建掩膜
% num_samples = length(t_grid);
% is_active_mask = false(num_samples, 1);
% for k = 1:length(segments)
%     is_active_mask(segments(k).start_idx : segments(k).end_idx) = true;
% end
% 
% % 2. 预计算接收机位置 & 卫星位置缓存 (这是计算密度的基础)
% rec_pos_cache = NaN(num_samples, 3);
% sat_az_el_cache = cell(num_samples, 1); % 存储每时刻所有卫星的 [Az, El]
% 
% fprintf('    预计算卫星几何分布...\n');
% for t = 1:num_samples
%     if is_active_mask(t)
%         [~, epoch_idx] = min(abs([obs_clean.time] - t_grid(t)));
%         try 
%             [rp, ~, sat_states] = calculate_receiver_position(obs_clean, nav_data, epoch_idx); 
%             rec_pos_cache(t,:) = rp;
%             
%             % 提取所有可见卫星的 Az/El (不仅仅是波动卫星，而是所有存在的卫星)
%             % 因为密度是相对于"整个星座分布"而言的
%             sats_geometry = struct('id', {}, 'az', {}, 'el', {});
%             cnt = 0;
%             [lat0, lon0, alt0] = ecef2geodetic(rp(1), rp(2), rp(3));
%             
%             fns = fieldnames(sat_states);
%             for i = 1:length(fns)
%                 sid = fns{i};
%                 sat_p = sat_states.(sid).position;
%                 [e, n, u] = ecef2enu(sat_p(1)-rp(1), sat_p(2)-rp(2), sat_p(3)-rp(3), lat0, lon0, alt0);
%                 az = atan2d(e, n);
%                 el = asind(u / norm([e, n, u]));
%                 if el > 0
%                     cnt = cnt + 1;
%                     sats_geometry(cnt).id = sid;
%                     sats_geometry(cnt).az = az;
%                     sats_geometry(cnt).el = el;
%                 end
%             end
%             sat_az_el_cache{t} = sats_geometry;
%         catch
%         end
%     end
% end
% 
% % 3. 主循环: 稀疏度加权追踪
% K_cluster = TRAJ.time_cluster_k;
% track_results = struct('t_idx', {}, 'x', {}, 'y', {}, 'total_energy', {});
% track_cnt = 0;
% is_tracking_started = false; 
% 
% fprintf('    [Processing] 执行稀疏度加权 (Sigma=%d°, Power=%.1f)...\n', ALG.density_sigma, ALG.sparsity_power);
% 
% for t_start = 1 : K_cluster : num_samples
%     t_end = min(t_start + K_cluster - 1, num_samples);
%     range_indices = t_start : t_end;
%     
%     if ~any(is_active_mask(range_indices))
%         is_tracking_started = false; 
%         continue;
%     end
%     
%     % --- A. 候选点收集 ---
%     window_candidates = struct('x', {}, 'y', {}, 'opp_score', {}, 'weight', {});
%     wc_cnt = 0;
%     
%     for t = range_indices
%         if ~is_active_mask(t), continue; end
%         rec_pos = rec_pos_cache(t, :);
%         if any(isnan(rec_pos)), continue; end
%         [lat0, lon0, alt0] = ecef2geodetic(rec_pos(1), rec_pos(2), rec_pos(3));
%         
%         % 1. 计算当前时刻的稀疏度权重表 (Sparsity Map)
%         geom = sat_az_el_cache{t};
%         if isempty(geom), continue; end
%         
%         % 构建 Az-El 矩阵
%         num_vis = length(geom);
%         coords = zeros(num_vis, 2); 
%         for i=1:num_vis, coords(i,:) = [geom(i).az, geom(i).el]; end
%         
%         % 计算密度权重
%         sparsity_weights = zeros(num_vis, 1);
%         for i = 1:num_vis
%             % 计算该卫星与其他所有卫星的"角距离" (简化为欧氏距离 approximation)
%             % 对于更严谨的球面距离，可用 distance() 函数，但这里 Az/El 欧氏距离足够
%             dists = sqrt(sum((coords - coords(i,:)).^2, 2));
%             
%             % 高斯核密度估计
%             local_density = sum(exp(-(dists.^2) / (ALG.density_sigma^2)));
%             
%             % 权重 = 1 / 密度 (越稀疏，权重越大)
%             sparsity_weights(i) = (1 / local_density) ^ ALG.sparsity_power;
%         end
%         
%         % 建立 ID -> Weight 映射
%         sw_map = containers.Map({geom.id}, num2cell(sparsity_weights));
%         
%         
%         % 2. 遍历波动卫星计算投影
%         current_vols = volatility_matrix(t, :);
%         valid_energy_idx = find(current_vols > ALG.min_sat_vol); 
%         if isempty(valid_energy_idx), continue; end
%         
%         [~, epoch_idx] = min(abs([obs_clean.time] - t_grid(t)));
%         try [~, ~, sat_states] = calculate_receiver_position(obs_clean, nav_data, epoch_idx); catch, continue; end
%         
%         for k = 1:length(valid_energy_idx)
%             s_idx = valid_energy_idx(k); 
%             sid = valid_sats{s_idx};
%             if ~isfield(sat_states, sid), continue; end
%             
%             % 几何投影
%             sat_p = sat_states.(sid).position;
%             [e, n, u] = ecef2enu(sat_p(1)-rec_pos(1), sat_p(2)-rec_pos(2), sat_p(3)-rec_pos(3), lat0, lon0, alt0);
%             vec_u = [e, n, u] / norm([e, n, u]); 
%             if vec_u(3) <= 0 || (90 - acosd(vec_u(3))) < TRAJ.min_elevation, continue; end
%             
%             t_int = TRAJ.gesture_height / vec_u(3); 
%             pt_int = t_int * vec_u;
%             px = pt_int(1); py = pt_int(2);
%             
%             % 获取该卫星的稀疏权重
%             w_sparse = 1.0;
%             if isKey(sw_map, sid)
%                 w_sparse = sw_map(sid);
%             end
%             
%             % [核心公式] 最终权重 = 信号遮挡强度 * 几何稀疏度
%             % 信号强说明挡得严实，稀疏说明位置珍贵
%             w_final = current_vols(s_idx) * w_sparse;
%             
%             % 依然保留前沿排序逻辑
%             vec_body_to_pt = [px, py] - ALG.body_pos;
%             opp_score = norm(vec_body_to_pt);
%             
%             wc_cnt = wc_cnt + 1;
%             window_candidates(wc_cnt).x = px;
%             window_candidates(wc_cnt).y = py;
%             window_candidates(wc_cnt).opp_score = opp_score;
%             window_candidates(wc_cnt).weight = w_final; % 使用修正后的权重
%         end
%     end
%     
%     if wc_cnt == 0, continue; end
%     
%     % --- B. 前沿提取 (Top 15%) ---
%     % 注意：现在的 Top 15% 是基于"修正权重"加权后的结果
%     all_scores = [window_candidates.opp_score];
%     [~, sort_idx] = sort(all_scores, 'descend'); 
%     
%     num_to_pick = max(1, ceil(wc_cnt * ALG.top_k_percent));
%     selected_idx = sort_idx(1:num_to_pick);
%     
%     sum_w = 0; sum_wx = 0; sum_wy = 0;
%     for k = 1:length(selected_idx)
%         idx = selected_idx(k);
%         w = window_candidates(idx).weight;
%         sum_w = sum_w + w;
%         sum_wx = sum_wx + window_candidates(idx).x * w;
%         sum_wy = sum_wy + window_candidates(idx).y * w;
%     end
%     
%     if sum_w == 0, continue; end
%     center_x = sum_wx / sum_w;
%     center_y = sum_wy / sum_w;
%     
%     % --- C. 动作锁 ---
%     if ~is_tracking_started
%         if norm([center_x, center_y]) > TRAJ.min_action_dist
%             is_tracking_started = true; 
%         else
%             continue; 
%         end
%     end
%     
%     track_cnt = track_cnt + 1;
%     track_results(track_cnt).t_idx = mean(range_indices);
%     track_results(track_cnt).x = center_x;
%     track_results(track_cnt).y = center_y;
%     track_results(track_cnt).total_energy = sum_w;
% end
% 
% % 提取轨迹
% if track_cnt > 0
%     traj_x = [track_results.x]'; 
%     traj_y = [track_results.y]'; 
%     traj_t = [track_results.t_idx]'; 
%     traj_e = [track_results.total_energy]';
%     
%     % PCA 线性规整 (保留作为可选后处理，让线条更好看)
%     if length(traj_x) > 5
%          traj_x = smoothdata(traj_x, 'movmean', 3);
%          traj_y = smoothdata(traj_y, 'movmean', 3);
%     end
% else
%     traj_x = []; traj_y = []; traj_t = []; traj_e = [];
% end
% 
% %% ================= [Part 3] 绘图流程 =================
% if ~isempty(traj_x)
%     figure('Name', 'Dynamic Sparsity Track v6.0', 'Position', [100, 200, 600, 600], 'Color', 'w');
%     ax = axes; hold(ax, 'on'); grid(ax, 'on'); axis(ax, 'equal');
%     xlabel('East (m)'); ylabel('North (m)');
%     
%     % 画参考点
%     plot(ax, 0, 0, '^', 'MarkerSize', 12, 'MarkerFaceColor', 'k', 'DisplayName', 'Receiver');
%     plot(ax, ALG.body_pos(1), ALG.body_pos(2), 'bs', 'MarkerSize', 10, 'DisplayName', 'Body Ref');
%     
%     % 画轨迹
%     % 颜色代表能量，点大小代表"稀疏度加权后的影响力"
%     scatter(ax, traj_x, traj_y, 30, traj_e, 'filled', 'DisplayName', 'Weighted Pts');
%     plot(ax, traj_x, traj_y, 'b-', 'LineWidth', 2.5, 'DisplayName', 'Balanced Path');
%     colormap(ax, 'parula'); c = colorbar; c.Label.String = 'Weighted Energy';
%     
%     % 起终点
%     plot(ax, traj_x(1), traj_y(1), 'go', 'MarkerFaceColor', 'g', 'MarkerSize', 10, 'DisplayName', 'Start');
%     plot(ax, traj_x(end), traj_y(end), 'rs', 'MarkerFaceColor', 'r', 'MarkerSize', 10, 'DisplayName', 'End');
%     
%     title({'v6.0 动态稀疏度加权 (Simulated Uniformity)', ...
%            sprintf('Sigma=%d° | Power=%.1f | De-Clustering', ALG.density_sigma, ALG.sparsity_power)});
%     legend('Location', 'best');
%     
%     xlim([-0.8, 0.8]); ylim([-1.2, 0.8]);
% else
%     fprintf('⚠️ 本次未生成有效轨迹。\n');
% end
% 
% fprintf('✅ v6.0 分析完成 (动态稀疏度加权)。\n');
% end











% % =========================================================================
% % run_gesture_analysis_boundary_trackV3.m (极值跨度追踪版)
% % 功能: 鲁棒手势感知 Step 2 - 基于物理跨度的指尖追踪 (v3.6)
% % 描述:
% %   放弃传统的“重心法”，采用“物理跨度极值法”来还原轨迹。
% %   核心假设：
% %     1. 手臂遮挡形成长条形阴影。
% %     2. 阴影中距离身体(Body Ref)最远的波动点，物理上对应指尖位置。
% %   算法流程：
% %     1. 建立参考系：身体坐标 B=[0, -1.0]。
% %     2. 投影计算：将所有遮挡卫星投影到 H=0.3m 平面。
% %     3. 极值提取：在每一时间窗口内，寻找距离 B 最远的 Top-K% 点。
% %     4. 坐标回归：利用这些最远点生成轨迹，自然克服重心内缩问题。
% %
% % [调用格式]:
% %   [traj_x, traj_y, traj_t, traj_e] = run_gesture_analysis_boundary_trackV3(obs_clean, nav_data, step1_res);
% % =========================================================================
% 
% function [traj_x, traj_y, traj_t, traj_e] = run_gesture_analysis_boundary_trackV3(obs_clean, nav_data, step1_res)
% 
% % --- 1. 数据解包 (Unpack Data) ---
% segments = step1_res.segments;
% volatility_matrix = step1_res.volatility_matrix;
% t_grid = step1_res.t_grid;
% valid_sats = step1_res.valid_sats;
% PARA = step1_res.PARA; 
% 
% if isempty(segments)
%     fprintf('⚠️ 警告: step1_res 中无有效分段，Step 2 (Extremum-Track) 终止。\n');
%     traj_x=[]; traj_y=[]; traj_t=[]; traj_e=[];
%     return;
% end
% 
% addpath(genpath('sky_plot')); 
% addpath(genpath('calculate_clock_bias_and_positon'));
% addpath(genpath('nav_parse'));
% 
% fprintf('--> [Step 2] 启动极值跨度追踪 (Physical Span / Farthest Point Logic)...\n');
% 
% %% ================= [Part 1] 参数设置 =================
% 
% % 1. 轨迹与几何参数
% TRAJ.gesture_height    = 0.30;  % [物理] 手势平面高度 (与仿真一致)
% TRAJ.min_elevation     = 10;    % [物理] 最低仰角门限 (度)
% TRAJ.min_action_dist   = 0.05;  % [触发] 动作死区 (米)
% 
% % 2. 身体参考参数
% ALG.body_pos           = [0.0, -1.0]; % [关键] 身体参考点 (接收机正南方1米)
% 
% % 3. 极值提取参数
% % [核心逻辑]: 既然手臂有宽度且是线状分布，重心必然偏向身体。
% % 我们只取距离身体"最远"的那一小撮点，它们才代表指尖。
% ALG.top_k_percent      = 0.05;  % [激进] 仅取最远的 5% (模拟寻找最远波动点)
% ALG.min_sat_vol        = 1.0;   % [过滤] 能量阈值
% 
% % 4. 时域聚类与平滑
% TRAJ.time_cluster_k    = 3;     % [聚类] 时间窗口 (减少窗口以提高极值灵敏度)
% TRAJ.traj_smooth_m     = 5;     % [平滑] 增加平滑以抵消 Top-5% 的噪声
% 
% %% ================= [Part 2] 核心计算流程 =================
% 
% % 1. 构建激活掩膜 (Mask)
% num_samples = length(t_grid);
% is_active_mask = false(num_samples, 1);
% for k = 1:length(segments)
%     is_active_mask(segments(k).start_idx : segments(k).end_idx) = true;
% end
% 
% % 2. 预计算接收机位置缓存
% rec_pos_cache = NaN(num_samples, 3);
% for t = 1:num_samples
%     if is_active_mask(t)
%         [~, epoch_idx] = min(abs([obs_clean.time] - t_grid(t)));
%         try [rp, ~, ~] = calculate_receiver_position(obs_clean, nav_data, epoch_idx); rec_pos_cache(t,:) = rp; catch, end
%     end
% end
% 
% % 3. 主循环: 极值跨度追踪
% K_cluster = TRAJ.time_cluster_k;
% fprintf('    [Processing] 正在提取距离身体最远的波动点 (Top %.1f%%)...\n', ALG.top_k_percent*100);
% 
% track_results = struct('t_idx', {}, 'x', {}, 'y', {}, 'total_energy', {});
% track_cnt = 0;
% is_tracking_started = false; 
% 
% for t_start = 1 : K_cluster : num_samples
%     
%     t_end = min(t_start + K_cluster - 1, num_samples);
%     range_indices = t_start : t_end;
%     
%     if ~any(is_active_mask(range_indices))
%         is_tracking_started = false; 
%         continue;
%     end
%     
%     % --- A. 候选点收集 (Candidates) ---
%     window_candidates = struct('x', {}, 'y', {}, 'dist_from_body', {}, 'weight', {});
%     wc_cnt = 0;
%     
%     for t = range_indices
%         if ~is_active_mask(t), continue; end
%         
%         rec_pos = rec_pos_cache(t, :);
%         if any(isnan(rec_pos)), continue; end
%         [lat0, lon0, alt0] = ecef2geodetic(rec_pos(1), rec_pos(2), rec_pos(3));
%         
%         current_vols = volatility_matrix(t, :);
%         valid_energy_idx = find(current_vols > ALG.min_sat_vol); 
%         
%         if isempty(valid_energy_idx), continue; end
%         
%         [~, epoch_idx] = min(abs([obs_clean.time] - t_grid(t)));
%         try [~, ~, sat_states] = calculate_receiver_position(obs_clean, nav_data, epoch_idx); catch, continue; end
%         
%         for k = 1:length(valid_energy_idx)
%             s_idx = valid_energy_idx(k); 
%             sid = valid_sats{s_idx};
%             if ~isfield(sat_states, sid), continue; end
%             
%             sat_p = sat_states.(sid).position;
%             [e, n, u] = ecef2enu(sat_p(1)-rec_pos(1), sat_p(2)-rec_pos(2), sat_p(3)-rec_pos(3), lat0, lon0, alt0);
%             vec_u = [e, n, u] / norm([e, n, u]); 
%             
%             if vec_u(3) <= 0 || (90 - acosd(vec_u(3))) < TRAJ.min_elevation
%                 continue; 
%             end
%             
%             % 计算投影点 P (在 0.3m 高度)
%             t_int = TRAJ.gesture_height / vec_u(3); 
%             pt_int = t_int * vec_u;
%             px = pt_int(1); py = pt_int(2);
%             
%             % --- [核心逻辑] 计算到身体参考点的物理跨度 ---
%             % 我们假设点离身体越远，越接近指尖
%             vec_body_to_pt = [px, py] - ALG.body_pos;
%             dist_val = norm(vec_body_to_pt);
%             
%             wc_cnt = wc_cnt + 1;
%             window_candidates(wc_cnt).x = px;
%             window_candidates(wc_cnt).y = py;
%             window_candidates(wc_cnt).dist_from_body = dist_val;
%             window_candidates(wc_cnt).weight = current_vols(s_idx);
%         end
%     end
%     
%     if wc_cnt == 0, continue; end
%     
%     % --- B. 极值提取 (Extraction) ---
%     % 按照“距离身体的远近”排序
%     all_dists = [window_candidates.dist_from_body];
%     [~, sort_idx] = sort(all_dists, 'descend'); % 降序：最远的排前面
%     
%     % 仅取 Top-K% (锁定指尖)
%     num_to_pick = max(1, ceil(wc_cnt * ALG.top_k_percent));
%     selected_idx = sort_idx(1:num_to_pick);
%     
%     sum_w = 0; sum_wx = 0; sum_wy = 0;
%     for k = 1:length(selected_idx)
%         idx = selected_idx(k);
%         w = window_candidates(idx).weight;
%         sum_w = sum_w + w;
%         sum_wx = sum_wx + window_candidates(idx).x * w;
%         sum_wy = sum_wy + window_candidates(idx).y * w;
%     end
%     
%     if sum_w == 0, continue; end
%     
%     center_x = sum_wx / sum_w;
%     center_y = sum_wy / sum_w;
%     
%     % --- C. 动作触发锁 ---
%     if ~is_tracking_started
%         if norm([center_x, center_y]) > TRAJ.min_action_dist
%             is_tracking_started = true; 
%         else
%             continue; 
%         end
%     end
%     
%     track_cnt = track_cnt + 1;
%     track_results(track_cnt).t_idx = mean(range_indices);
%     track_results(track_cnt).x = center_x;
%     track_results(track_cnt).y = center_y;
%     track_results(track_cnt).total_energy = sum_w;
% end
% 
% % --- 3. 轨迹生成与平滑 ---
% if track_cnt > 0
%     traj_x = [track_results.x]'; 
%     traj_y = [track_results.y]'; 
%     traj_t = [track_results.t_idx]'; 
%     traj_e = [track_results.total_energy]';
%     
%     % [平滑策略]
%     % 由于取极值点 (Top 5%) 可能会带来高频噪声，这里使用中值滤波去除尖峰，再滑动平均
%     if length(traj_x) > 5
%         traj_x = medfilt1(traj_x, 3);
%         traj_y = medfilt1(traj_y, 3);
%         
%         traj_x = smoothdata(traj_x, 'movmean', TRAJ.traj_smooth_m);
%         traj_y = smoothdata(traj_y, 'movmean', TRAJ.traj_smooth_m);
%     end
% else
%     traj_x = []; traj_y = []; traj_t = []; traj_e = [];
% end
% 
% %% ================= [Part 3] 绘图流程 =================
% if ~isempty(traj_x)
%     figure('Name', 'Extremum Span Track v3.6', 'Position', [100, 200, 600, 600], 'Color', 'w');
%     ax = axes; hold(ax, 'on'); grid(ax, 'on'); axis(ax, 'equal');
%     xlabel('East (m)'); ylabel('North (m)');
%     
%     % 画参考点
%     plot(ax, 0, 0, '^', 'MarkerSize', 12, 'MarkerFaceColor', 'k', 'DisplayName', 'Receiver');
%     plot(ax, ALG.body_pos(1), ALG.body_pos(2), 'bs', 'MarkerSize', 10, 'DisplayName', 'Body Ref');
%     
%     % 画轨迹
%     % 用颜色表示能量，辅助判断轨迹可信度
%     scatter(ax, traj_x, traj_y, 30, traj_e, 'filled', 'DisplayName', 'Farthest Points');
%     plot(ax, traj_x, traj_y, 'b-', 'LineWidth', 2, 'DisplayName', 'Fingertip Path');
%     colormap(ax, 'parula'); c = colorbar; c.Label.String = 'Cluster Energy';
%     
%     % 起终点
%     plot(ax, traj_x(1), traj_y(1), 'go', 'MarkerFaceColor', 'g', 'MarkerSize', 10, 'DisplayName', 'Start');
%     plot(ax, traj_x(end), traj_y(end), 'rs', 'MarkerFaceColor', 'r', 'MarkerSize', 10, 'DisplayName', 'End');
%     
%     title({'v3.6 极值跨度追踪 (Farthest Point Track)', ...
%            sprintf('Height=%.1fm | Top %.1f%% Farthest From Body', TRAJ.gesture_height, ALG.top_k_percent*100)});
%     legend('Location', 'best');
%     
%     xlim([-0.8, 0.8]); ylim([-1.2, 0.8]);
% else
%     fprintf('⚠️ 本次未生成有效轨迹。\n');
% end
% 
% fprintf('✅ v3.6 极值跨度分析完成 (已返回指尖轨迹)。\n');
% end








% =========================================================================
% run_gesture_analysis_boundary_trackV3.m (对抗增强版)
% 功能: 鲁棒手势感知 Step 2 - 身体对抗/前沿追踪版 (v3.5)
% 描述:
%   针对“手臂遮挡重心偏移”问题进行了重构。
%   算法不再简单提取距离圆心最远的点，而是建立“身体-手部”对抗坐标系：
%   1. 锁定身体参考点 (Receiver 南方 1m)。
%   2. 计算各投影点相对于身体的方位矢量。
%   3. 提取在该矢量方向上投影值最大的 Top-K% 采样点。
%   4. 从而强制算法追踪“手臂的最远端”（即指尖），找回字母顶点。
% =========================================================================

function [traj_x, traj_y, traj_t, traj_e] = run_gesture_analysis_boundary_trackV3(obs_clean, nav_data, step1_res)

% --- 1. 数据解包 (Unpack Data) ---
segments = step1_res.segments;
volatility_matrix = step1_res.volatility_matrix;
t_grid = step1_res.t_grid;
valid_sats = step1_res.valid_sats;
PARA = step1_res.PARA; 

if isempty(segments)
    fprintf('⚠️ 警告: step1_res 中无有效分段，Step 2 (Boundary-Opposition) 终止。\n');
    traj_x=[]; traj_y=[]; traj_t=[]; traj_e=[];
    return;
end

addpath(genpath('sky_plot')); 
addpath(genpath('calculate_clock_bias_and_positon'));
addpath(genpath('nav_parse'));

fprintf('--> [Step 2] 启动对抗增强前沿追踪 (Opposition-Based Tip Tracking)...\n');

%% ================= [Part 1] 参数设置 =================

% 1. 轨迹与几何参数
TRAJ.gesture_height    = 0.30;  % [同步] 与仿真一致的高度 (0.3m)
TRAJ.min_elevation     = 10;    % [物理] 最低仰角门限 (度)
TRAJ.min_action_dist   = 0.05;  % [触发] 动作死区 (米)

% 2. 身体参考参数
ALG.body_pos           = [0.0, -1.0]; % [关键] 身体/肩膀参考点位置 (EN轴)

% 3. 抗干扰与边界算法参数
ALG.top_k_percent      = 0.15;  % [核心] 压缩至 15%，仅保留最前端的指尖信息
ALG.min_sat_vol        = 1.0;   % [过滤] 单星波动门限 (dB)

% 4. 时域聚类与平滑
TRAJ.time_cluster_k    = 5;     % [聚类] 时间窗口 (5点聚合)
TRAJ.traj_smooth_m     = 3;     % [平滑] 增加平滑度以应对 Top-K 的跳变

%% ================= [Part 2] 核心计算流程 =================

% 1. 构建激活掩膜 (Mask)
num_samples = length(t_grid);
is_active_mask = false(num_samples, 1);
for k = 1:length(segments)
    is_active_mask(segments(k).start_idx : segments(k).end_idx) = true;
end

% 2. 预计算接收机位置缓存
rec_pos_cache = NaN(num_samples, 3);
for t = 1:num_samples
    if is_active_mask(t)
        [~, epoch_idx] = min(abs([obs_clean.time] - t_grid(t)));
        try [rp, ~, ~] = calculate_receiver_position(obs_clean, nav_data, epoch_idx); rec_pos_cache(t,:) = rp; catch, end
    end
end

% 3. 主循环: 对抗边界追踪
K_cluster = TRAJ.time_cluster_k;
fprintf('    [Processing] 执行对抗投影过滤 (Body-Ref: [0, -1])...\n');

track_results = struct('t_idx', {}, 'x', {}, 'y', {}, 'total_energy', {});
track_cnt = 0;
is_tracking_started = false; 



for t_start = 1 : K_cluster : num_samples
    
    t_end = min(t_start + K_cluster - 1, num_samples);
    range_indices = t_start : t_end;
    
    if ~any(is_active_mask(range_indices))
        is_tracking_started = false; 
        continue;
    end
    
    % --- A. 候选点收集 ---
    window_candidates = struct('x', {}, 'y', {}, 'opp_score', {}, 'weight', {});
    wc_cnt = 0;
    
    for t = range_indices
        if ~is_active_mask(t), continue; end
        
        rec_pos = rec_pos_cache(t, :);
        if any(isnan(rec_pos)), continue; end
        [lat0, lon0, alt0] = ecef2geodetic(rec_pos(1), rec_pos(2), rec_pos(3));
        
        % 读取能量矩阵
        current_vols = volatility_matrix(t, :);
        valid_energy_idx = find(current_vols > ALG.min_sat_vol); 
        
        if isempty(valid_energy_idx), continue; end
        
        [~, epoch_idx] = min(abs([obs_clean.time] - t_grid(t)));
        try [~, ~, sat_states] = calculate_receiver_position(obs_clean, nav_data, epoch_idx); catch, continue; end
        
        for k = 1:length(valid_energy_idx)
            s_idx = valid_energy_idx(k); 
            sid = valid_sats{s_idx};
            if ~isfield(sat_states, sid), continue; end
            
            sat_p = sat_states.(sid).position;
            [e, n, u] = ecef2enu(sat_p(1)-rec_pos(1), sat_p(2)-rec_pos(2), sat_p(3)-rec_pos(3), lat0, lon0, alt0);
            vec_u = [e, n, u] / norm([e, n, u]); 
            
            if vec_u(3) <= 0 || (90 - acosd(vec_u(3))) < TRAJ.min_elevation
                continue; 
            end
            
            % 计算在 0.3m 平面上的投影点
            t_int = TRAJ.gesture_height / vec_u(3); 
            pt_int = t_int * vec_u;
            px = pt_int(1); py = pt_int(2);
            
            % --- [核心] 对抗得分计算 ---
            % 计算投影点相对于身体参考点的矢量
            vec_body_to_pt = [px, py] - ALG.body_pos;
            % 得分定义：距离身体越远的点，得分越高（越可能是指尖）
            opp_score = norm(vec_body_to_pt);
            
            wc_cnt = wc_cnt + 1;
            window_candidates(wc_cnt).x = px;
            window_candidates(wc_cnt).y = py;
            window_candidates(wc_cnt).opp_score = opp_score;
            window_candidates(wc_cnt).weight = current_vols(s_idx);
        end
    end
    
    if wc_cnt == 0, continue; end
    
    % --- B. 前沿提取 (Top-K Opposition) ---
    all_scores = [window_candidates.opp_score];
    [~, sort_idx] = sort(all_scores, 'descend'); 
    
    num_to_pick = max(1, ceil(wc_cnt * ALG.top_k_percent));
    selected_idx = sort_idx(1:num_to_pick);
    
    sum_w = 0; sum_wx = 0; sum_wy = 0;
    for k = 1:length(selected_idx)
        idx = selected_idx(k);
        w = window_candidates(idx).weight;
        sum_w = sum_w + w;
        sum_wx = sum_wx + window_candidates(idx).x * w;
        sum_wy = sum_wy + window_candidates(idx).y * w;
    end
    
    if sum_w == 0, continue; end
    % 计算加权中心
    center_x = sum_wx / sum_w;
    center_y = sum_wy / sum_w;
    
    % --- C. 动作触发锁 ---
    if ~is_tracking_started
        if norm([center_x, center_y]) > TRAJ.min_action_dist
            is_tracking_started = true; 
        else
            continue; 
        end
    end
    
    track_cnt = track_cnt + 1;
    track_results(track_cnt).t_idx = mean(range_indices);
    track_results(track_cnt).x = center_x;
    track_results(track_cnt).y = center_y;
    track_results(track_cnt).total_energy = sum_w;
end

% --- 3. 轨迹提取与平滑 ---
if track_cnt > 0
    traj_x = [track_results.x]'; 
    traj_y = [track_results.y]'; 
    traj_t = [track_results.t_idx]'; 
    traj_e = [track_results.total_energy]';
    
    % 物理偏置手动微调 (根据手臂宽度 40cm 估算，补偿重心内缩)
    traj_y = traj_y + 0.05; 
    
    if TRAJ.traj_smooth_m > 1
        traj_x = smoothdata(traj_x, 'movmean', TRAJ.traj_smooth_m);
        traj_y = smoothdata(traj_y, 'movmean', TRAJ.traj_smooth_m);
    end
else
    traj_x = []; traj_y = []; traj_t = []; traj_e = [];
end

%% ================= [Part 3] 绘图流程 =================
if ~isempty(traj_x)
    figure('Name', 'Opposition Boundary Track v3.5', 'Position', [100, 200, 600, 600], 'Color', 'w');
    ax = axes; hold(ax, 'on'); grid(ax, 'on'); axis(ax, 'equal');
    xlabel('East (m)'); ylabel('North (m)');
    
    % 画参考点
    plot(ax, 0, 0, '^', 'MarkerSize', 12, 'MarkerFaceColor', 'k', 'DisplayName', 'Receiver');
    plot(ax, ALG.body_pos(1), ALG.body_pos(2), 'bs', 'MarkerSize', 10, 'DisplayName', 'Body Ref');
    
    % 画轨迹
    plot(ax, traj_x, traj_y, 'b-', 'LineWidth', 2.5, 'DisplayName', 'Tip Trajectory');
    scatter(ax, traj_x, traj_y, 40, traj_e, 'filled', 'DisplayName', 'Leading Pts');
    
    % 起终点
    plot(ax, traj_x(1), traj_y(1), 'go', 'MarkerFaceColor', 'g', 'MarkerSize', 10, 'DisplayName', 'Start');
    plot(ax, traj_x(end), traj_y(end), 'rs', 'MarkerFaceColor', 'r', 'MarkerSize', 10, 'DisplayName', 'End');
    
    title({'v3.5 对抗前沿追踪 (Opposition Tip Track)', ...
           sprintf('Height=%.1fm | Top %.0f%% Leading Edge', TRAJ.gesture_height, ALG.top_k_percent*100)});
    legend('Location', 'best');
    
    xlim([-0.8, 0.8]); ylim([-1.2, 0.8]);
else
    fprintf('⚠️ 本次未生成有效轨迹。\n');
end

fprintf('✅ v3.5 对抗分析完成 (已返回指尖前沿轨迹数据)。\n');
end










% % =========================================================================
% % run_gesture_analysis_boundary_trackV3.m (函数版)
% % 功能: 鲁棒手势感知 Step 2 - 边界/前沿追踪版 (v3.3 Struct Input)
% % 描述:
% %   该函数是 Step 2 的一种实现，针对"远端手势" (如伸出手画图) 设计。
% %   它接收 Step 1 处理后的干净数据，不再计算所有遮挡点的重心，而是筛选
% %   距离身体中心"最远"的一批点 (Top-K%)，以此模拟追踪"指尖"的位置，
% %   有效避免了手掌或手臂遮挡信号将轨迹重心向身体侧拉扯的问题。
% %
% % [调用格式]:
% %   [traj_x, traj_y, traj_t, traj_e] = run_gesture_analysis_boundary_trackV3(obs_clean, nav_data, step1_res);
% %
% % [输入参数]:
% %   1. obs_clean (struct): Step 1 返回的清洗后观测数据 (用于定位)。
% %   2. nav_data (struct): 导航星历数据。
% %   3. step1_res (struct): Step 1 返回的结果包 (含能量矩阵、分段信息等)。
% %
% % [返回值说明]:
% %   1. traj_x / traj_y: [指尖/前沿] 的轨迹坐标 (米)。
% %   2. traj_t: 轨迹点时间索引。
% %   3. traj_e: 前沿聚类点的总能量。
% %
% % [核心算法]:
% %   1. Top-K Filtering: 仅保留距离圆心最远的 K% (如50%) 的点参与计算。
% %   2. Zenith Safe Mask: 强制剔除低仰角区域，防止身体核心区遮挡干扰。
% % =========================================================================
% 
% function [traj_x, traj_y, traj_t, traj_e] = run_gesture_analysis_boundary_trackV3(obs_clean, nav_data, step1_res)
% 
% % --- 1. 数据解包 (Unpack Data) ---
% segments = step1_res.segments;
% volatility_matrix = step1_res.volatility_matrix;
% t_grid = step1_res.t_grid;
% valid_sats = step1_res.valid_sats;
% PARA = step1_res.PARA; % 获取 Step 1 的参数
% 
% if isempty(segments)
%     fprintf('⚠️ 警告: step1_res 中无有效分段，Step 2 (Boundary) 终止。\n');
%     traj_x=[]; traj_y=[]; traj_t=[]; traj_e=[];
%     return;
% end
% 
% addpath(genpath('sky_plot')); 
% addpath(genpath('calculate_clock_bias_and_positon'));
% addpath(genpath('nav_parse'));
% 
% fprintf('--> [Step 2] 启动前沿边界追踪 (Function版 v3.3: Boundary/Top-K)...\n');
% 
% %% ================= [Part 1] 参数设置 =================
% 
% % 1. 轨迹与几何参数
% TRAJ.gesture_height    = 0.30;  % [物理] 手势平面高度 (米)
% TRAJ.min_elevation     = 10;    % [物理] 最低仰角门限 (度)
% TRAJ.min_action_dist   = 0.05;  % [触发] 动作死区 (米)
% 
% % 2. 过滤参数
% % [修复] 使用 Step 1 传入的 PARA 进行兼容，或直接定义局部参数
% PARA.min_sat_vol       = 1.0;   % [过滤] 单星波动门限 (dB)
% 
% % 3. 抗干扰与边界算法参数
% ALG.zenith_safe_deg    = 10;    % [掩膜] 安全仰角门限 (度): 剔除身体遮挡
% ALG.top_k_percent      = 0.5;   % [核心] 前沿追踪比例: 仅取最远的 50% 点
% 
% % 4. 时域聚类
% TRAJ.time_cluster_k    = 5;     % [聚类] 时间窗口
% TRAJ.traj_smooth_m     = 2;     % [平滑] 轨迹平滑窗口
% 
% %% ================= [Part 2] 核心计算流程 =================
% 
% % 1. 构建激活掩膜 (Mask)
% num_samples = length(t_grid);
% is_active_mask = false(num_samples, 1);
% for k = 1:length(segments)
%     is_active_mask(segments(k).start_idx : segments(k).end_idx) = true;
% end
% 
% % 2. 预计算接收机位置缓存
% rec_pos_cache = NaN(num_samples, 3);
% for t = 1:num_samples
%     if is_active_mask(t)
%         [~, epoch_idx] = min(abs([obs_clean.time] - t_grid(t)));
%         try [rp, ~, ~] = calculate_receiver_position(obs_clean, nav_data, epoch_idx); rec_pos_cache(t,:) = rp; catch, end
%     end
% end
% 
% % 3. 主循环: 边界追踪
% K_cluster = TRAJ.time_cluster_k;
% fprintf('--> 执行边界追踪 (Zenith<%d°, Top%.0f%%)...\n', ALG.zenith_safe_deg, ALG.top_k_percent*100);
% 
% track_results = struct('t_idx', {}, 'x', {}, 'y', {}, 'total_energy', {});
% track_cnt = 0;
% is_tracking_started = false; 
% 
% for t_start = 1 : K_cluster : num_samples
%     
%     t_end = min(t_start + K_cluster - 1, num_samples);
%     range_indices = t_start : t_end;
%     
%     if ~any(is_active_mask(range_indices))
%         is_tracking_started = false; 
%         continue;
%     end
%     
%     % --- A. 候选点收集 ---
%     window_candidates = struct('x', {}, 'y', {}, 'dist', {}, 'weight', {});
%     wc_cnt = 0;
%     sum_window_energy = 0;
%     
%     for t = range_indices
%         if ~is_active_mask(t), continue; end
%         
%         rec_pos = rec_pos_cache(t, :);
%         if any(isnan(rec_pos)), continue; end
%         [lat0, lon0, alt0] = ecef2geodetic(rec_pos(1), rec_pos(2), rec_pos(3));
%         
%         % 读取能量矩阵
%         current_vols = volatility_matrix(t, :);
%         valid_energy_idx = find(current_vols > PARA.min_sat_vol); 
%         
%         if isempty(valid_energy_idx), continue; end
%         
%         [~, epoch_idx] = min(abs([obs_clean.time] - t_grid(t)));
%         try [~, ~, sat_states] = calculate_receiver_position(obs_clean, nav_data, epoch_idx); catch, continue; end
%         
%         for k = 1:length(valid_energy_idx)
%             s_idx = valid_energy_idx(k); 
%             sid = valid_sats{s_idx};
%             if ~isfield(sat_states, sid), continue; end
%             
%             sat_p = sat_states.(sid).position;
%             [e, n, u] = ecef2enu(sat_p(1)-rec_pos(1), sat_p(2)-rec_pos(2), sat_p(3)-rec_pos(3), lat0, lon0, alt0);
%             dist = norm([e, n, u]);
%             vec_u = [e, n, u] / dist; 
%             
%             zen_deg = acosd(vec_u(3)); 
%             el_deg  = 90 - zen_deg;
%             
%             if vec_u(3) <= 0, continue; end 
%             
%             % --- [安全掩膜] 剔除身体遮挡 ---
%             if el_deg < ALG.zenith_safe_deg
%                 continue; 
%             end
%             
%             % 计算投影点
%             t_int = TRAJ.gesture_height / vec_u(3); 
%             pt_int = t_int * vec_u;
%             
%             % 半径门控
%             dist_from_center = norm(pt_int(1:2));
%             if dist_from_center > 5.0, continue; end
%             
%             wc_cnt = wc_cnt + 1;
%             window_candidates(wc_cnt).x = pt_int(1);
%             window_candidates(wc_cnt).y = pt_int(2);
%             window_candidates(wc_cnt).dist = dist_from_center;
%             
%             w = current_vols(s_idx);
%             window_candidates(wc_cnt).weight = w;
%             sum_window_energy = sum_window_energy + w;
%         end
%     end
%     
%     if wc_cnt == 0, continue; end
%     
%     % --- B. 边界提取 (Top-K Logic) ---
%     % 核心思想: 离圆心最远的点通常对应"指尖"位置
%     all_dists = [window_candidates.dist];
%     [~, sort_idx] = sort(all_dists, 'descend'); 
%     
%     num_to_pick = max(1, ceil(wc_cnt * ALG.top_k_percent));
%     selected_idx = sort_idx(1:num_to_pick);
%     
%     sum_w = 0; sum_wx = 0; sum_wy = 0;
%     for k = 1:length(selected_idx)
%         idx = selected_idx(k);
%         w = window_candidates(idx).weight;
%         sum_w = sum_w + w;
%         sum_wx = sum_wx + window_candidates(idx).x * w;
%         sum_wy = sum_wy + window_candidates(idx).y * w;
%     end
%     
%     if sum_w == 0, continue; end
%     center_x = sum_wx / sum_w;
%     center_y = sum_wy / sum_w;
%     
%     % --- C. 动作触发锁 ---
%     if ~is_tracking_started
%         if norm([center_x, center_y]) > TRAJ.min_action_dist
%             is_tracking_started = true; 
%         else
%             continue; 
%         end
%     end
%     
%     track_cnt = track_cnt + 1;
%     track_results(track_cnt).t_idx = mean(range_indices);
%     track_results(track_cnt).x = center_x;
%     track_results(track_cnt).y = center_y;
%     track_results(track_cnt).total_energy = sum_window_energy;
% end
% 
% % 轨迹提取
% if track_cnt > 0
%     traj_x = [track_results.x]'; 
%     traj_y = [track_results.y]'; 
%     traj_t = [track_results.t_idx]'; 
%     traj_e = [track_results.total_energy]';
%     
%     if TRAJ.traj_smooth_m > 1
%         traj_x = smoothdata(traj_x, 'movmean', TRAJ.traj_smooth_m);
%         traj_y = smoothdata(traj_y, 'movmean', TRAJ.traj_smooth_m);
%     end
% else
%     traj_x = []; traj_y = []; traj_t = []; traj_e = [];
% end
% 
% %% ================= [Part 3] 绘图流程 =================
% fprintf('\n--> 开始生成图表...\n');
% 
% if ~isempty(traj_x)
%     figure('Name', 'Reconstructed Boundary Trajectory v3.3', 'Position', [100, 200, 600, 600], 'Color', 'w');
%     ax = axes; hold(ax, 'on'); grid(ax, 'on'); axis(ax, 'equal');
%     xlabel('East (m)'); ylabel('North (m)');
%     
%     % 画接收机
%     plot(ax, 0, 0, '^', 'MarkerSize', 12, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'DisplayName', 'Receiver');
%     viscircles(ax, [0,0], TRAJ.min_action_dist, 'Color', [0.8 0.8 0.8], 'LineStyle', '--'); 
%     
%     % 画轨迹 (带能量颜色)
%     plot(ax, traj_x, traj_y, 'b-', 'LineWidth', 2, 'DisplayName', 'Finger Path');
%     scatter(ax, traj_x, traj_y, 40, traj_e, 'filled', 'DisplayName', 'Cluster Point');
%     c = colorbar; c.Label.String = 'Cluster Energy'; colormap(ax, 'parula');
%     
%     % 起终点
%     plot(ax, traj_x(1), traj_y(1), 'go', 'MarkerFaceColor', 'g', 'MarkerSize', 8, 'DisplayName', 'Start');
%     plot(ax, traj_x(end), traj_y(end), 'rs', 'MarkerFaceColor', 'r', 'MarkerSize', 8, 'DisplayName', 'End');
%     
%     title({'Boundary Track v3.3', ...
%            sprintf('Elev > %d° | Top %.0f%% | Cluster K=%d', ALG.zenith_safe_deg, ALG.top_k_percent*100, TRAJ.time_cluster_k)});
%     legend('Location', 'best');
%     
%     max_range = max(max(abs(traj_x)), max(abs(traj_y)));
%     if max_range < 0.5, max_range = 0.5; end
%     xlim([-max_range*1.2, max_range*1.2]); ylim([-max_range*1.2, max_range*1.2]);
% else
%     fprintf('⚠️ 本次未生成有效轨迹。\n');
% end
% 
% fprintf('✅ v3.3 边界追踪分析完成 (已返回轨迹数据)。\n');
% end