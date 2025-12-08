% =========================================================================
% run_gesture_analysis_robust.m
% 功能: 鲁棒手势感知全流程 (v10.2 - 参数更新版)
% 核心: 
%   1. 参数完全匹配用户指定要求 (GVI=5, Window=1.5, SpikeDur=3)。
%   2. 逻辑保持 v10.1 的修复状态 (无数据覆盖 Bug)。
% =========================================================================

%% ================= [Part 1] 环境检查与参数设置 =================
clearvars -except obs_data nav_data; 
if ~exist('obs_data', 'var'), error('请先加载 obs_data!'); end

addpath(genpath('sky_plot')); 
addpath(genpath('calculate_clock_bias_and_positon'));
addpath(genpath('nav_parse'));

fprintf('--> 启动鲁棒手势感知分析 (Parameters Updated)...\n');

% --- 参数设置 ---
% [Step 0: 信号预处理参数 (Smart Flattening)]
PARA.diff_lag_N         = 5;     % 趋势窗口 (N点): 用于识别缓坡台阶，防止被误削
PARA.noise_cutoff_db    = 1;     % 噪声阈值 (dB): 波动幅度 <= 此值时视为噪声被削平
PARA.spike_th_db        = 2;     % 毛刺幅度阈值 (dB): 跳变 > 此值才检查毛刺
PARA.spike_max_duration = 3;     % 毛刺最大持续点数 (1-5): 允许剔除宽毛刺
% [Step 1: 分段与检测参数]
PARA.smooth_window_sec = 1.5;   % 滑动平均窗口大小（秒）
PARA.gvi_threshold     = 5;     % GVI（波动指数）阈值
PARA.sampling_rate     = 25;    % 数据采样率（Hz）
PARA.merge_gap_sec     = 0.01;  % 合并间隙（秒）
PARA.min_duration_sec  = 0.4;   % 最小持续时间（秒）
% [Step 2: 3D 轨迹投影参数]
TRAJ.gesture_height    = 0.20;  % 手势物理高度（米）
TRAJ.min_elevation     = 15;    % 最低仰角阈值（度）
TRAJ.energy_th_ratio   = 0.2;   % 能量筛选比例
% [Step 2: 抗干扰与几何算法参数]
ALG.zenith_safe_deg    = 30;    % 天顶角安全区（度）
ALG.az_neighbor_dist   = 20;    % 方位角邻域（度）
ALG.density_penalty_k  = 1.0;   % 密度惩罚系数
ALG.ransac_iter        = 500;   % RANSAC 迭代次数
ALG.ransac_dist_th     = 0.20;  % RANSAC 内点阈值（米）


%% ================= [Part 2] 核心计算流程 =================

% ----------------- [Step 0 计算] 数据提取与智能削平 -----------------
fprintf('--> [Step 0] 提取数据并执行智能削平...\n');

% 1. 提取卫星ID
all_sat_ids = {};
for i = 1:min(100, length(obs_data)), if ~isempty(obs_data(i).data), all_sat_ids = [all_sat_ids, fieldnames(obs_data(i).data)']; end; end
unique_sat_ids = unique(all_sat_ids);
valid_sats = {};
for i = 1:length(unique_sat_ids)
    sid = unique_sat_ids{i};
    if ismember(sid(1), ['G','C','E','J']), valid_sats{end+1} = sid; end
end

% 2. 构建时间轴
raw_times = [obs_data.time];
t_grid = (min(raw_times) : seconds(1/PARA.sampling_rate) : max(raw_times))'; 
num_samples = length(t_grid);
num_sats = length(valid_sats);
t_grid_plot = t_grid + hours(8) - seconds(18); 

% 3. 提取原始数据 (DATA ENTRY POINT)
cn0_matrix = NaN(num_samples, num_sats);
for s_idx = 1:num_sats
    sat_id = valid_sats{s_idx};
    target_snr_code = '';
    for k = 1:min(50, length(obs_data))
        if isfield(obs_data(k).data, sat_id) && isfield(obs_data(k).data.(sat_id), 'snr')
            fields = fieldnames(obs_data(k).data.(sat_id).snr);
            if ismember('S1C', fields), target_snr_code = 'S1C'; elseif ismember('S2I', fields), target_snr_code = 'S2I'; else, if ~isempty(fields), target_snr_code = fields{1}; end; end
            if ~isempty(target_snr_code), break; end
        end
    end
    if isempty(target_snr_code), continue; end
    s_times = []; s_cn0 = [];
    for k = 1:length(obs_data)
        if isfield(obs_data(k).data, sat_id) && isfield(obs_data(k).data.(sat_id).snr, target_snr_code)
            val = obs_data(k).data.(sat_id).snr.(target_snr_code);
            if ~isnan(val) && val > 10, s_times = [s_times; obs_data(k).time]; s_cn0 = [s_cn0; val]; end
        end
    end
    if length(s_times) > 20
        [u_times, u_idx] = unique(s_times);
        cn0_matrix(:, s_idx) = interp1(u_times, s_cn0(u_idx), t_grid, 'linear', NaN);
    end
end

% 4. 智能削平算法 (直接修改 cn0_matrix)
% DATA FLOW: Raw CN0 -> Flattened CN0
N = PARA.diff_lag_N; Th = PARA.noise_cutoff_db; SpikeTh = PARA.spike_th_db; MaxDur = PARA.spike_max_duration;
for s = 1:num_sats
    raw_col = cn0_matrix(:, s); col = raw_col;
    last_val = NaN; skip_counter = 0; 
    for t = 1:num_samples
        curr_raw = raw_col(t);
        if skip_counter > 0, col(t) = last_val; skip_counter = skip_counter - 1; continue; end
        if isnan(curr_raw), last_val = NaN; continue; end
        if isnan(last_val), last_val = curr_raw; col(t) = last_val; continue; end
        
        % 宽毛刺检测
        is_spike = false; spike_width = 0; delta_curr = curr_raw - last_val;
        if abs(delta_curr) > SpikeTh
            for k = 1:MaxDur
                if t + k > num_samples, break; end
                if abs(raw_col(t+k) - last_val) <= Th, is_spike = true; spike_width = k; break; end
            end
        end
        if is_spike, col(t) = last_val; skip_counter = spike_width - 1; continue; end
        
        % 趋势判断
        if t > N, trend_diff = abs(curr_raw - raw_col(t-N)); else, trend_diff = 0; end
        if trend_diff > Th
            col(t) = curr_raw; last_val = curr_raw;
        else
            if abs(curr_raw - last_val) <= Th, col(t) = last_val; else, col(t) = curr_raw; last_val = curr_raw; end
        end
    end
    cn0_matrix(:, s) = col; % 更新矩阵为削平后的数据
end


% ----------------- [Step 1 计算] SG滤波与分段 -----------------
fprintf('--> [Step 1] SG滤波与手势分段 (基于削平数据)...\n');

% [数据流确认] 此时 cn0_matrix 已经是削平过的了，直接送入 SG
for s = 1:num_sats
    col = cn0_matrix(:, s); valid = ~isnan(col);
    if sum(valid) > 14
        idx = 1:length(col); filled = interp1(idx(valid), col(valid), idx, 'linear', 'extrap')';
        cn0_matrix(:, s) = sgolayfilt(filled, 2, 7); % 再次更新为 SG 滤波后的数据
    end
end

% 计算 GVI (基于 削平+SG 后的数据)
cn0_smooth = movmean(cn0_matrix, round(PARA.smooth_window_sec * PARA.sampling_rate), 1, 'omitnan');
volatility_matrix = abs(cn0_matrix - cn0_smooth);
gvi_curve_clean = movmean(sum(volatility_matrix, 2, 'omitnan'), 5);
is_active = gvi_curve_clean > PARA.gvi_threshold;

% 连通域
gap_pts = round(PARA.merge_gap_sec * PARA.sampling_rate);
pad_act = [1; is_active; 1];
g_starts = find(diff(pad_act)==-1); g_ends = find(diff(pad_act)==1)-1;
for i=1:length(g_starts), if (g_ends(i)-g_starts(i)+1) < gap_pts && g_starts(i)>1, is_active(g_starts(i):g_ends(i)-1) = 1; end; end
edges = diff([0; is_active; 0]); s_idxs = find(edges==1); e_idxs = find(edges==-1)-1;
min_dur = round(PARA.min_duration_sec * PARA.sampling_rate);
segments = struct('id', {}, 'start_idx', {}, 'end_idx', {}, 'peak_time', {}, 'peak_gvi', {}, 'peak_idx', {});
cnt = 0;
for i=1:length(s_idxs)
    if (e_idxs(i)-s_idxs(i)) >= min_dur
        cnt = cnt + 1;
        [m_v, m_i] = max(gvi_curve_clean(s_idxs(i):e_idxs(i)));
        segments(cnt).id = cnt; segments(cnt).start_idx = s_idxs(i); segments(cnt).end_idx = e_idxs(i);
        segments(cnt).peak_time = t_grid(s_idxs(i)+m_i-1); segments(cnt).peak_gvi = m_v;
        segments(cnt).peak_idx = s_idxs(i) + m_i - 1; 
    end
end

% 初始化结果容器
segment_fits = []; step2_vis_data = []; final_draw_data = []; seg_cnt = 0;

if cnt == 0
    fprintf('⚠️  提示: 本次未检测到有效手势片段 (GVI均低于 %d)。\n', PARA.gvi_threshold);
else
    fprintf('✅ [Step 1] 分段完成，共 %d 个片段。\n', cnt);
    
    % ----------------- [Step 2 计算] 鲁棒 3D 轨迹推演 -----------------
    fprintf('--> [Step 2] 开始鲁棒轨迹推演...\n');
    [~, anchor_ep_idx] = min(abs([obs_data.time] - segments(1).peak_time));
    try
        [ref_rec_pos, ~, ~] = calculate_receiver_position(obs_data, nav_data, anchor_ep_idx);
        if isempty(ref_rec_pos) || any(isnan(ref_rec_pos)), error('无法计算参考点位置'); end
        [ref_lat, ref_lon, ref_alt] = ecef2geodetic(ref_rec_pos(1), ref_rec_pos(2), ref_rec_pos(3));
    catch ME
        fprintf('⚠️ 参考点计算失败: %s\n', ME.message);
        cnt = 0; 
    end
    
    if cnt > 0
        segment_fits = struct('p_start', {}, 'p_end', {}, 't_center', {}, 'w_sum', {}, 'valid', {});
        step2_vis_data = struct('seg_id', {}, 'traj_az', {}, 'hits_data', {}, 'best_inliers', {}, 'p_start', {}, 'p_end', {});
        
        for i = 1:length(segments)
            seg = segments(i);
            idx_range = seg.start_idx : seg.end_idx;
            seg_times = t_grid(idx_range);
            % sub_vol 也是基于削平后的数据计算的
            sub_vol = volatility_matrix(idx_range, :);
            
            [~, ep_idx] = min(abs([obs_data.time] - seg.peak_time));
            try
                [~, ~, sat_states] = calculate_receiver_position(obs_data, nav_data, ep_idx);
            catch
                continue;
            end
            
            hits_data = struct('pos', {}, 'w_final', {}, 't_off', {}, 'id', {}, 'zen_deg', {});
            raw_sats = struct('pos', {}, 'az', {}, 'zen', {}, 'energy', {}, 't_off', {}, 'id', {});
            raw_cnt = 0; hit_cnt = 0;
            
            for s = 1:length(valid_sats)
                sid = valid_sats{s}; if ~isfield(sat_states, sid), continue; end
                sat_p = sat_states.(sid).position;
                [e, n, u] = ecef2enu(sat_p(1)-ref_rec_pos(1), sat_p(2)-ref_rec_pos(2), sat_p(3)-ref_rec_pos(3), ref_lat, ref_lon, ref_alt);
                dist = norm([e, n, u]); vec_u = [e, n, u]/dist;
                zen_deg = acosd(vec_u(3)); az_deg  = atan2d(e, n); if az_deg<0, az_deg=az_deg+360; end
                
                if vec_u(3) <= 0 || (90-zen_deg) < TRAJ.min_elevation, continue; end
                t_int = TRAJ.gesture_height / vec_u(3); pt_int = t_int * vec_u;
                if norm(pt_int(1:2)) > 3.0, continue; end 
                
                % energy 也是基于削平数据的 sub_vol
                energy = sum(sub_vol(:, s), 'omitnan'); [~, mx_i] = max(sub_vol(:, s));
                
                raw_cnt = raw_cnt + 1;
                raw_sats(raw_cnt).pos = pt_int; raw_sats(raw_cnt).az = az_deg; raw_sats(raw_cnt).zen = zen_deg;
                raw_sats(raw_cnt).energy = energy; raw_sats(raw_cnt).t_off = seconds(seg_times(mx_i) - seg_times(1));
                raw_sats(raw_cnt).id = sid;
            end
            
            if raw_cnt < 2, continue; end
            all_e = [raw_sats.energy]; th_e = max(all_e) * TRAJ.energy_th_ratio;
            hit_candidates = raw_sats(all_e > th_e);
            if length(hit_candidates) < 3, continue; end
            
            for k = 1:length(hit_candidates)
                cand = hit_candidates(k); w_base = cosd(cand.zen); penalty = 1.0;
                if cand.zen > ALG.zenith_safe_deg 
                    n_neighbors = 0;
                    for j = 1:length(hit_candidates)
                        if k==j, continue; end
                        if abs(cand.az - hit_candidates(j).az) < ALG.az_neighbor_dist || abs(cand.az - hit_candidates(j).az) > 360-ALG.az_neighbor_dist, n_neighbors=n_neighbors+1; end
                    end
                    penalty = 1.0 / (1.0 + ALG.density_penalty_k * n_neighbors);
                end
                hit_cnt = hit_cnt + 1;
                hits_data(hit_cnt).pos = cand.pos; hits_data(hit_cnt).t_off = cand.t_off; hits_data(hit_cnt).id = cand.id;
                hits_data(hit_cnt).zen_deg = cand.zen; hits_data(hit_cnt).w_final = w_base * penalty;
            end
            
            pts_xy = vertcat(hits_data.pos); pts_xy = pts_xy(:, 1:2); weights = [hits_data.w_final]';
            best_score = -1; best_inliers = false(hit_cnt, 1);
            num_pts = size(pts_xy, 1);
            for iter = 1:ALG.ransac_iter
                if num_pts >= 2, sample_idx = randsample(num_pts, 2, true, weights); else, sample_idx=[1,2]; end
                p1 = pts_xy(sample_idx(1), :); p2 = pts_xy(sample_idx(2), :);
                vec = p2 - p1; if norm(vec) < 1e-3, continue; end 
                vec = vec / norm(vec); normal = [-vec(2), vec(1)]; 
                dists = abs((pts_xy - p1) * normal');
                is_inlier = dists < ALG.ransac_dist_th; current_score = sum(weights(is_inlier));
                if current_score > best_score, best_score = current_score; best_inliers = is_inlier; end
            end
            if sum(best_inliers) < 2, continue; end
            
            inlier_pts = pts_xy(best_inliers, :); inlier_w = weights(best_inliers);
            w_sum = sum(inlier_w); mean_w = sum(inlier_pts .* inlier_w) / w_sum;
            centered = inlier_pts - mean_w; weighted_centered = centered .* sqrt(inlier_w);
            [U_svd, ~, ~] = svd(weighted_centered' * weighted_centered);
            dir_final = U_svd(:, 1)'; 
            
            t_inliers = [hits_data(best_inliers).t_off]'; proj_vals = centered * dir_final';
            corr_v = corr(proj_vals, t_inliers);
            if ~isnan(corr_v) && corr_v < 0, dir_final = -dir_final; end
            
            traj_az = atan2d(dir_final(1), dir_final(2)); if traj_az < 0, traj_az = traj_az + 360; end
            fprintf('   Seg #%d: 鲁棒方向 %.1f° (Corr: %.2f)\n', i, traj_az, corr_v);
            
            dists_along_axis = centered * dir_final';
            d_min = min(dists_along_axis); d_max = max(dists_along_axis);
            p_start_real = mean_w + d_min * dir_final;
            p_end_real   = mean_w + d_max * dir_final;
            
            seg_cnt = seg_cnt + 1;
            segment_fits(seg_cnt).p_start = p_start_real;
            segment_fits(seg_cnt).p_end   = p_end_real;
            segment_fits(seg_cnt).t_center = seg.peak_time;
            
            step2_vis_data(seg_cnt).seg_id = i; step2_vis_data(seg_cnt).traj_az = traj_az;
            step2_vis_data(seg_cnt).hits_data = hits_data; step2_vis_data(seg_cnt).best_inliers = best_inliers;
            step2_vis_data(seg_cnt).p_start = p_start_real; step2_vis_data(seg_cnt).p_end = p_end_real;
        end
        
        % ----------------- [Step 3 计算] 轨迹重构 -----------------
        fprintf('--> [Step 3] 生成最终手势矢量链...\n');
        if seg_cnt > 0
            [~, sort_idx] = sort([segment_fits.t_center]); sorted_segs = segment_fits(sort_idx);
            current_pen_pos = sorted_segs(1).p_start; 
            for k = 1:seg_cnt
                raw_seg = sorted_segs(k);
                vec_motion = raw_seg.p_end - raw_seg.p_start;
                new_start = current_pen_pos; new_end = current_pen_pos + vec_motion;
                final_draw_data(k).start = new_start; final_draw_data(k).end = new_end;
                final_draw_data(k).vec = vec_motion; final_draw_data(k).id = k;
                current_pen_pos = new_end; 
            end
        end
    end
end


%% ================= [Part 3] 统一绘图流程 =================
fprintf('\n--> 开始生成所有图表...\n');

% --- [图表 1] Step 1 GVI 分段详情图 ---
if ~isempty(segments)
    figure('Name', 'Step 1: GVI Segmentation', 'Position', [50, 100, 1000, 600], 'Color', 'w');
    subplot(2,1,1);
    plot(t_grid_plot, gvi_curve_clean, 'k-', 'LineWidth', 1); hold on;
    yline(PARA.gvi_threshold, 'b--', 'Threshold');
    title(sprintf('全局波动指数 (GVI) 与阈值 (Th=%d)', PARA.gvi_threshold));
    ylabel('GVI'); xlabel('Time (BJT)');
    datetick('x', 'HH:MM:ss', 'keepticks', 'keeplimits'); grid on; axis tight;
    subplot(2,1,2);
    plot(t_grid_plot, gvi_curve_clean, 'Color', [0.8 0.8 0.8]); hold on;
    yline(PARA.gvi_threshold, 'b--');
    yl = ylim;
    for i = 1:cnt
        idx = segments(i).start_idx : segments(i).end_idx;
        t_s = t_grid_plot(segments(i).start_idx); t_e = t_grid_plot(segments(i).end_idx);
        patch([t_s t_e t_e t_s], [yl(1) yl(1) yl(2) yl(2)], 'r', 'FaceAlpha', 0.1, 'EdgeColor', 'none');
        plot(t_grid_plot(idx), gvi_curve_clean(idx), 'r-', 'LineWidth', 1.5);
        text(t_grid_plot(segments(i).peak_idx), segments(i).peak_gvi, sprintf('#%d', i), 'Color', 'r', 'FontWeight', 'bold', 'VerticalAlignment', 'bottom');
    end
    title(sprintf('检测到的手势片段: %d 个 (若为0则表示全是噪声)', cnt));
    xlabel('Time (BJT)'); ylabel('GVI');
    datetick('x', 'HH:MM:ss', 'keepticks', 'keeplimits'); grid on; axis tight;
end

% --- [图表 2 & 3] ---
if cnt > 0
    for k = 1:length(step2_vis_data)
        d = step2_vis_data(k);
        fig_name = sprintf('Seg #%d Analysis (Az=%.1f)', d.seg_id, d.traj_az);
        f = figure('Name', fig_name, 'Position', [100+(k-1)*30, 200, 500, 500], 'Color', 'w');
        ax = axes('Parent', f); hold(ax, 'on'); grid(ax, 'on'); axis(ax, 'equal');
        xlabel('East (m)'); ylabel('North (m)');
        title({sprintf('手势片段 #%d 几何分析', d.seg_id), '红色=Inliers | 灰色=Outliers | 蓝色=拟合向量'});
        plot(ax, 0, 0, '^', 'MarkerSize', 10, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k'); 
        for j = 1:length(d.hits_data)
            pt = d.hits_data(j).pos; w = d.hits_data(j).w_final;
            is_in = d.best_inliers(j);
            col = [0.7 0.7 0.7]; if is_in, col = 'r'; end
            ms = 5 + w * 15; 
            plot(ax, pt(1), pt(2), 'o', 'MarkerFaceColor', col, 'MarkerEdgeColor', 'k', 'MarkerSize', ms);
        end
        plot(ax, [d.p_start(1), d.p_end(1)], [d.p_start(2), d.p_end(2)], 'b-', 'LineWidth', 2);
        quiver(ax, d.p_start(1), d.p_start(2), d.p_end(1)-d.p_start(1), d.p_end(2)-d.p_start(2), 'Color', 'b', 'LineWidth', 2, 'MaxHeadSize', 0.5, 'AutoScale', 'off');
        xlim([-2.5 2.5]); ylim([-2.5 2.5]);
    end
    
    if seg_cnt > 0
        f_final = figure('Name', 'Final Gesture Vector Map', 'Position', [300, 100, 800, 800], 'Color', 'w');
        ax_f = axes('Parent', f_final); hold(ax_f, 'on'); grid(ax_f, 'on'); axis(ax_f, 'equal');
        xlabel('East (m)'); ylabel('North (m)');
        title({'最终手势矢量重构图', '首尾无缝拼接 (保留原始方向和长度)'});
        plot(ax_f, 0, 0, '^', 'MarkerSize', 12, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'DisplayName', 'Receiver');
        colors = lines(seg_cnt);
        start_p = final_draw_data(1).start;
        plot(ax_f, start_p(1), start_p(2), 'kp', 'MarkerSize', 12, 'MarkerFaceColor', 'g', 'DisplayName', 'Start');
        for k = 1:seg_cnt
            d = final_draw_data(k); col = colors(k, :);
            plot(ax_f, [d.start(1), d.end(1)], [d.start(2), d.end(2)], '-', 'Color', col, 'LineWidth', 4, 'DisplayName', sprintf('Stroke #%d', k));
            quiver(ax_f, d.start(1), d.start(2), d.vec(1), d.vec(2), 'Color', col, 'LineWidth', 2, 'MaxHeadSize', 0.4, 'AutoScale', 'off', 'HandleVisibility', 'off');
            text(ax_f, d.start(1), d.start(2), sprintf('%d', k), 'Color', 'w', 'FontWeight', 'bold', 'FontSize', 10, 'BackgroundColor', col, 'HorizontalAlignment','center');
        end
        axis(ax_f, 'tight'); xl = xlim(ax_f); yl = ylim(ax_f);
        xlim(ax_f, xl + [-0.5 0.5]); ylim(ax_f, yl + [-0.5 0.5]); legend('Location', 'bestoutside');
    end
end
fprintf('✅ 所有分析完成。\n');









% 原流程：原始数据插值 -> SG 滤波 -> GVI 计算 -> 分段...


% % =========================================================================
% % run_gesture_analysis_robust.m
% % 功能: 鲁棒手势感知全流程 (v8.0 - 计算与绘图分离版)
% % 结构:
% %   [Part 1] 参数设置 & 环境准备
% %   [Part 2] 核心计算流程 (Step 1 -> Step 2 -> Step 3)
% %   [Part 3] 统一绘图流程 (所有 Figure 在最后生成)
% % =========================================================================
% 
% %% ================= [Part 1] 环境检查与参数设置 =================
% clearvars -except obs_data nav_data; 
% if ~exist('obs_data', 'var'), error('请先加载 obs_data!'); end
% 
% addpath(genpath('sky_plot')); 
% addpath(genpath('calculate_clock_bias_and_positon'));
% addpath(genpath('nav_parse'));
% 
% fprintf('--> 启动鲁棒手势感知分析 (Logic-Plot Separation)...\n');
% 
% % --- 参数设置 (严格保留用户指定参数) ---
% % [Step 1: 分段与检测参数]
% PARA.smooth_window_sec = 1.8;   % 滑动平均窗口大小（秒）：用于计算信号的平滑基线，从而提取波动量。
% PARA.gvi_threshold     = 8;    % GVI（波动指数）阈值：信号波动强度超过此值才会被标记为潜在的手势活动。
% PARA.sampling_rate     = 25;    % 数据采样率（Hz）：用于将时间转换为采样点数，需与接收机设置一致。
% PARA.merge_gap_sec     = 0.01;  % 合并间隙（秒）：两个动作片段间隔小于此值时，视为同一个连续动作合并（防止动作断裂）。
% PARA.min_duration_sec  = 0.4;   % 最小持续时间（秒）：持续时间短于此值的波动被视为偶然噪声并丢弃。
% 
% % [Step 2: 3D 轨迹投影参数]
% TRAJ.gesture_height    = 0.20;  % 手势物理高度（米）：假设手在接收机天线上方 0.3 米的平面上运动（用于视线投影计算）。
% TRAJ.min_elevation     = 15;    % 最低仰角阈值（度）：仰角低于此值（即天顶角 > 75度）的卫星因多径严重或投影误差大而被剔除。
% TRAJ.energy_th_ratio   = 0.2;   % 能量筛选比例：只有波动能量 > (最大能量 * 0.4) 的卫星才会被判定为 Hit（有效触点）。
% 
% % [Step 2: 抗干扰与几何算法参数]
% ALG.zenith_safe_deg    = 30;    % 天顶角安全区（度）：天顶角小于此值（即头顶区域）的卫星，被认为不受手臂干扰，不进行密度惩罚。
% ALG.az_neighbor_dist   = 20;    % 方位角邻域（度）：定义“邻居”的范围。如果侧面卫星在此范围内有其他 Hit 点，视为聚集干扰。
% ALG.density_penalty_k  = 1.0;   % 密度惩罚系数：用于打压侧面聚集的卫星权重。值越大，对手臂干扰的抑制越狠。
% ALG.ransac_iter        = 500;   % RANSAC 迭代次数：随机采样拟合直线的尝试次数，次数越多越容易找到最优解。
% ALG.ransac_dist_th     = 0.20;  % RANSAC 内点阈值（米）：投影点到拟合直线的距离小于此值时，被判定为“Inlier”（有效手势点）。
% 
% 
% %% ================= [Part 2] 核心计算流程 =================
% 
% % ----------------- [Step 1 计算] 数据提取、滤波与分段 -----------------
% fprintf('--> [Step 1 计算] 提取全星座数据...\n');
% all_sat_ids = {};
% for i = 1:min(100, length(obs_data)), if ~isempty(obs_data(i).data), all_sat_ids = [all_sat_ids, fieldnames(obs_data(i).data)']; end; end
% unique_sat_ids = unique(all_sat_ids);
% valid_sats = {};
% for i = 1:length(unique_sat_ids)
%     sid = unique_sat_ids{i};
%     if ismember(sid(1), ['G','C','E','J']), valid_sats{end+1} = sid; end
% end
% 
% raw_times = [obs_data.time];
% t_grid = (min(raw_times) : seconds(1/PARA.sampling_rate) : max(raw_times))'; 
% num_samples = length(t_grid);
% num_sats = length(valid_sats);
% t_grid_plot = t_grid + hours(8) - seconds(18); % 北京时间轴(含闰秒修正)
% 
% cn0_matrix = NaN(num_samples, num_sats);
% for s_idx = 1:num_sats
%     sat_id = valid_sats{s_idx};
%     target_snr_code = '';
%     for k = 1:min(50, length(obs_data))
%         if isfield(obs_data(k).data, sat_id) && isfield(obs_data(k).data.(sat_id), 'snr')
%             fields = fieldnames(obs_data(k).data.(sat_id).snr);
%             if ismember('S1C', fields), target_snr_code = 'S1C'; elseif ismember('S2I', fields), target_snr_code = 'S2I'; else, if ~isempty(fields), target_snr_code = fields{1}; end; end
%             if ~isempty(target_snr_code), break; end
%         end
%     end
%     if isempty(target_snr_code), continue; end
%     s_times = []; s_cn0 = [];
%     for k = 1:length(obs_data)
%         if isfield(obs_data(k).data, sat_id) && isfield(obs_data(k).data.(sat_id).snr, target_snr_code)
%             val = obs_data(k).data.(sat_id).snr.(target_snr_code);
%             if ~isnan(val) && val > 10, s_times = [s_times; obs_data(k).time]; s_cn0 = [s_cn0; val]; end
%         end
%     end
%     if length(s_times) > 20
%         [u_times, u_idx] = unique(s_times);
%         cn0_matrix(:, s_idx) = interp1(u_times, s_cn0(u_idx), t_grid, 'linear', NaN);
%     end
% end
% 
% % SG 预滤波
% for s = 1:num_sats
%     col = cn0_matrix(:, s); valid = ~isnan(col);
%     if sum(valid) > 14
%         idx = 1:length(col); filled = interp1(idx(valid), col(valid), idx, 'linear', 'extrap')';
%         cn0_matrix(:, s) = sgolayfilt(filled, 2, 7);
%     end
% end
% 
% % 分段逻辑
% cn0_smooth = movmean(cn0_matrix, round(PARA.smooth_window_sec * PARA.sampling_rate), 1, 'omitnan');
% volatility_matrix = abs(cn0_matrix - cn0_smooth);
% gvi_curve_clean = movmean(sum(volatility_matrix, 2, 'omitnan'), 5);
% is_active = gvi_curve_clean > PARA.gvi_threshold;
% gap_pts = round(PARA.merge_gap_sec * PARA.sampling_rate);
% pad_act = [1; is_active; 1];
% g_starts = find(diff(pad_act)==-1); g_ends = find(diff(pad_act)==1)-1;
% for i=1:length(g_starts), if (g_ends(i)-g_starts(i)+1) < gap_pts && g_starts(i)>1, is_active(g_starts(i):g_ends(i)-1) = 1; end; end
% edges = diff([0; is_active; 0]); s_idxs = find(edges==1); e_idxs = find(edges==-1)-1;
% min_dur = round(PARA.min_duration_sec * PARA.sampling_rate);
% segments = struct('id', {}, 'start_idx', {}, 'end_idx', {}, 'peak_time', {}, 'peak_gvi', {}, 'peak_idx', {});
% cnt = 0;
% for i=1:length(s_idxs)
%     if (e_idxs(i)-s_idxs(i)) >= min_dur
%         cnt = cnt + 1;
%         [m_v, m_i] = max(gvi_curve_clean(s_idxs(i):e_idxs(i)));
%         segments(cnt).id = cnt; segments(cnt).start_idx = s_idxs(i); segments(cnt).end_idx = e_idxs(i);
%         segments(cnt).peak_time = t_grid(s_idxs(i)+m_i-1); segments(cnt).peak_gvi = m_v;
%         segments(cnt).peak_idx = s_idxs(i) + m_i - 1; % [关键修复]
%     end
% end
% fprintf('✅ [Step 1] 分段完成，共 %d 个片段。\n', cnt);
% 
% 
% % ----------------- [Step 2 计算] 鲁棒 3D 轨迹推演 -----------------
% fprintf('--> [Step 2 计算] 开始鲁棒轨迹推演 (锁定全局参考系)...\n');
% 
% if isempty(segments), error('未检测到有效手势片段'); end
% [~, anchor_ep_idx] = min(abs([obs_data.time] - segments(1).peak_time));
% try
%     [ref_rec_pos, ~, ~] = calculate_receiver_position(obs_data, nav_data, anchor_ep_idx);
%     if isempty(ref_rec_pos) || any(isnan(ref_rec_pos)), error('无法计算参考点位置'); end
%     [ref_lat, ref_lon, ref_alt] = ecef2geodetic(ref_rec_pos(1), ref_rec_pos(2), ref_rec_pos(3));
% catch ME
%     error('参考点计算失败: %s', ME.message);
% end
% 
% % 初始化结果容器
% segment_fits = struct('p_start', {}, 'p_end', {}, 't_center', {}, 'w_sum', {}, 'valid', {});
% % [新增] 用于 Step 2 绘图的数据缓存容器
% step2_vis_data = struct('seg_id', {}, 'traj_az', {}, 'hits_data', {}, 'best_inliers', {}, 'p_start', {}, 'p_end', {});
% 
% seg_cnt = 0;
% 
% for i = 1:length(segments)
%     seg = segments(i);
%     idx_range = seg.start_idx : seg.end_idx;
%     seg_times = t_grid(idx_range);
%     sub_vol = volatility_matrix(idx_range, :);
%     
%     [~, ep_idx] = min(abs([obs_data.time] - seg.peak_time));
%     try
%         [~, ~, sat_states] = calculate_receiver_position(obs_data, nav_data, ep_idx);
%     catch
%         continue;
%     end
%     
%     hits_data = struct('pos', {}, 'w_final', {}, 't_off', {}, 'id', {}, 'zen_deg', {});
%     raw_sats = struct('pos', {}, 'az', {}, 'zen', {}, 'energy', {}, 't_off', {}, 'id', {});
%     raw_cnt = 0; hit_cnt = 0;
%     
%     % 收集卫星数据
%     for s = 1:length(valid_sats)
%         sid = valid_sats{s};
%         if ~isfield(sat_states, sid), continue; end
%         sat_p = sat_states.(sid).position;
%         [e, n, u] = ecef2enu(sat_p(1)-ref_rec_pos(1), sat_p(2)-ref_rec_pos(2), sat_p(3)-ref_rec_pos(3), ref_lat, ref_lon, ref_alt);
%         dist = norm([e, n, u]); vec_u = [e, n, u]/dist;
%         zen_deg = acosd(vec_u(3)); az_deg  = atan2d(e, n); if az_deg<0, az_deg=az_deg+360; end
%         
%         if vec_u(3) <= 0 || (90-zen_deg) < TRAJ.min_elevation, continue; end
%         t_int = TRAJ.gesture_height / vec_u(3);
%         pt_int = t_int * vec_u;
%         if norm(pt_int(1:2)) > 3.0, continue; end 
%         energy = sum(sub_vol(:, s), 'omitnan'); [~, mx_i] = max(sub_vol(:, s));
%         
%         raw_cnt = raw_cnt + 1;
%         raw_sats(raw_cnt).pos = pt_int; raw_sats(raw_cnt).az = az_deg; raw_sats(raw_cnt).zen = zen_deg;
%         raw_sats(raw_cnt).energy = energy; raw_sats(raw_cnt).t_off = seconds(seg_times(mx_i) - seg_times(1));
%         raw_sats(raw_cnt).id = sid;
%     end
%     
%     if raw_cnt < 2, continue; end
%     all_e = [raw_sats.energy]; th_e = max(all_e) * TRAJ.energy_th_ratio;
%     hit_candidates = raw_sats(all_e > th_e);
%     if length(hit_candidates) < 3, continue; end
%     
%     % 抗干扰权重
%     for k = 1:length(hit_candidates)
%         cand = hit_candidates(k);
%         w_base = cosd(cand.zen); 
%         penalty = 1.0;
%         if cand.zen > ALG.zenith_safe_deg 
%             n_neighbors = 0;
%             for j = 1:length(hit_candidates)
%                 if k == j, continue; end
%                 d_az = abs(cand.az - hit_candidates(j).az);
%                 if d_az > 180, d_az = 360 - d_az; end
%                 if d_az < ALG.az_neighbor_dist, n_neighbors = n_neighbors + 1; end
%             end
%             penalty = 1.0 / (1.0 + ALG.density_penalty_k * n_neighbors);
%         end
%         hit_cnt = hit_cnt + 1;
%         hits_data(hit_cnt).pos = cand.pos; hits_data(hit_cnt).t_off = cand.t_off; hits_data(hit_cnt).id = cand.id;
%         hits_data(hit_cnt).zen_deg = cand.zen; hits_data(hit_cnt).w_final = w_base * penalty;
%     end
%     
%     % RANSAC
%     pts_xy = vertcat(hits_data.pos); pts_xy = pts_xy(:, 1:2); weights = [hits_data.w_final]';
%     best_score = -1; best_inliers = false(hit_cnt, 1);
%     num_pts = size(pts_xy, 1);
%     for iter = 1:ALG.ransac_iter
%         if num_pts >= 2, sample_idx = randsample(num_pts, 2, true, weights); else, sample_idx=[1,2]; end
%         p1 = pts_xy(sample_idx(1), :); p2 = pts_xy(sample_idx(2), :);
%         vec = p2 - p1;
%         if norm(vec) < 1e-3, continue; end 
%         vec = vec / norm(vec); normal = [-vec(2), vec(1)]; 
%         dists = abs((pts_xy - p1) * normal');
%         is_inlier = dists < ALG.ransac_dist_th;
%         current_score = sum(weights(is_inlier));
%         if current_score > best_score, best_score = current_score; best_inliers = is_inlier; end
%     end
%     if sum(best_inliers) < 2, continue; end
%     
%     % PCA & Correlation
%     inlier_pts = pts_xy(best_inliers, :); inlier_w = weights(best_inliers);
%     w_sum = sum(inlier_w); mean_w = sum(inlier_pts .* inlier_w) / w_sum;
%     centered = inlier_pts - mean_w; weighted_centered = centered .* sqrt(inlier_w);
%     [U_svd, ~, ~] = svd(weighted_centered' * weighted_centered);
%     dir_final = U_svd(:, 1)'; 
%     
%     t_inliers = [hits_data(best_inliers).t_off]'; proj_vals = centered * dir_final';
%     corr_v = corr(proj_vals, t_inliers);
%     if ~isnan(corr_v) && corr_v < 0, dir_final = -dir_final; end
%     
%     traj_az = atan2d(dir_final(1), dir_final(2)); if traj_az < 0, traj_az = traj_az + 360; end
%     fprintf('   Seg #%d: 鲁棒方向 %.1f° (Corr: %.2f)\n', i, traj_az, corr_v);
%     
%     % 计算线段
%     dists_along_axis = centered * dir_final';
%     d_min = min(dists_along_axis); d_max = max(dists_along_axis);
%     p_start_real = mean_w + d_min * dir_final;
%     p_end_real   = mean_w + d_max * dir_final;
%     
%     seg_cnt = seg_cnt + 1;
%     segment_fits(seg_cnt).p_start = p_start_real;
%     segment_fits(seg_cnt).p_end   = p_end_real;
%     segment_fits(seg_cnt).t_center = seg.peak_time;
%     
%     % [绘图数据保存]
%     step2_vis_data(seg_cnt).seg_id = i;
%     step2_vis_data(seg_cnt).traj_az = traj_az;
%     step2_vis_data(seg_cnt).hits_data = hits_data;
%     step2_vis_data(seg_cnt).best_inliers = best_inliers;
%     step2_vis_data(seg_cnt).p_start = p_start_real;
%     step2_vis_data(seg_cnt).p_end = p_end_real;
% end
% 
% 
% % ----------------- [Step 3 计算] 轨迹重构 (相对位移拼接) -----------------
% fprintf('--> [Step 3 计算] 生成最终手势矢量链 (Stitching)...\n');
% final_draw_data = struct('start', {}, 'end', {}, 'vec', {}, 'id', {});
% if seg_cnt > 0
%     [~, sort_idx] = sort([segment_fits.t_center]);
%     sorted_segs = segment_fits(sort_idx);
%     current_pen_pos = sorted_segs(1).p_start; % 初始下笔点
%     
%     for k = 1:seg_cnt
%         raw_seg = sorted_segs(k);
%         vec_motion = raw_seg.p_end - raw_seg.p_start;
%         new_start = current_pen_pos;
%         new_end   = current_pen_pos + vec_motion;
%         
%         final_draw_data(k).start = new_start;
%         final_draw_data(k).end = new_end;
%         final_draw_data(k).vec = vec_motion;
%         final_draw_data(k).id = k;
%         
%         current_pen_pos = new_end; % 移动笔触
%     end
% end
% 
% 
% %% ================= [Part 3] 统一绘图流程 =================
% fprintf('\n--> 开始生成所有图表...\n');
% 
% % --- [图表 1] Step 1 GVI 分段详情图 ---
% % 功能: 展示全时段的信号波动情况，以及算法识别出的手势起止片段。
% % 红色背景块表示识别出的动作区间，波峰数字对应片段ID。
% if ~isempty(segments)
%     figure('Name', 'Step 1: GVI Segmentation', 'Position', [50, 100, 1000, 600], 'Color', 'w');
%     subplot(2,1,1);
%     plot(t_grid_plot, gvi_curve_clean, 'k-', 'LineWidth', 1); hold on;
%     yline(PARA.gvi_threshold, 'b--', 'Threshold');
%     title(sprintf('全局波动指数 (GVI) 与阈值 (Th=%d)', PARA.gvi_threshold));
%     ylabel('GVI'); xlabel('Time (BJT)');
%     datetick('x', 'HH:MM:ss', 'keepticks', 'keeplimits'); grid on; axis tight;
% 
%     subplot(2,1,2);
%     plot(t_grid_plot, gvi_curve_clean, 'Color', [0.8 0.8 0.8]); hold on;
%     yline(PARA.gvi_threshold, 'b--');
%     yl = ylim;
%     for i = 1:cnt
%         idx = segments(i).start_idx : segments(i).end_idx;
%         t_s = t_grid_plot(segments(i).start_idx);
%         t_e = t_grid_plot(segments(i).end_idx);
%         patch([t_s t_e t_e t_s], [yl(1) yl(1) yl(2) yl(2)], 'r', 'FaceAlpha', 0.1, 'EdgeColor', 'none');
%         plot(t_grid_plot(idx), gvi_curve_clean(idx), 'r-', 'LineWidth', 1.5);
%         text(t_grid_plot(segments(i).peak_idx), segments(i).peak_gvi, sprintf('#%d', i), ...
%             'Color', 'r', 'FontWeight', 'bold', 'VerticalAlignment', 'bottom');
%     end
%     title('检测到的手势片段 (红色高亮)');
%     xlabel('Time (BJT)'); ylabel('GVI');
%     datetick('x', 'HH:MM:ss', 'keepticks', 'keeplimits'); grid on; axis tight;
% end
% 
% % --- [图表 2] Step 2 单片段几何分析图 (循环生成) ---
% % 功能: 逐个展示每个手势片段的微观几何结构。
% % 红色点 = RANSAC 筛选出的内点 (Inliers)，代表手部信号。
% % 灰色点 = 被剔除的噪声/手臂干扰 (Outliers)。
% % 蓝色箭头 = 最终拟合出的该片段运动向量。
% 
% % for k = 1:length(step2_vis_data)
% %     d = step2_vis_data(k);
% %     fig_name = sprintf('Seg #%d Analysis (Az=%.1f)', d.seg_id, d.traj_az);
% %     f = figure('Name', fig_name, 'Position', [100+(k-1)*30, 200, 500, 500], 'Color', 'w');
% %     ax = axes('Parent', f); hold(ax, 'on'); grid(ax, 'on'); axis(ax, 'equal');
% %     xlabel('East (m)'); ylabel('North (m)');
% %     title({sprintf('手势片段 #%d 几何分析', d.seg_id), '红色=Inliers | 灰色=Outliers | 蓝色=拟合向量'});
% %     
% %     plot(ax, 0, 0, '^', 'MarkerSize', 10, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k'); % 接收机
% %     
% %     for j = 1:length(d.hits_data)
% %         pt = d.hits_data(j).pos; w = d.hits_data(j).w_final;
% %         is_in = d.best_inliers(j);
% %         col = [0.7 0.7 0.7]; if is_in, col = 'r'; end
% %         ms = 5 + w * 15; 
% %         plot(ax, pt(1), pt(2), 'o', 'MarkerFaceColor', col, 'MarkerEdgeColor', 'k', 'MarkerSize', ms);
% %     end
% %     
% %     plot(ax, [d.p_start(1), d.p_end(1)], [d.p_start(2), d.p_end(2)], 'b-', 'LineWidth', 2);
% %     quiver(ax, d.p_start(1), d.p_start(2), ...
% %            d.p_end(1)-d.p_start(1), d.p_end(2)-d.p_start(2), ...
% %            'Color', 'b', 'LineWidth', 2, 'MaxHeadSize', 0.5, 'AutoScale', 'off');
% %     xlim([-2.5 2.5]); ylim([-2.5 2.5]);
% % end
% 
% 
% % --- [图表 3] Step 3 最终手势矢量重构图 ---
% % 功能: 将所有片段按时间顺序首尾相接，还原完整的手势轨迹。
% % 忽略了片段之间的绝对位置跳变，只保留“方向”和“长度”信息。
% % 实线 = 动作行程；箭头 = 运动方向。
% if seg_cnt > 0
%     f_final = figure('Name', 'Final Gesture Vector Map', 'Position', [300, 100, 800, 800], 'Color', 'w');
%     ax_f = axes('Parent', f_final);
%     hold(ax_f, 'on'); grid(ax_f, 'on'); axis(ax_f, 'equal');
%     xlabel('East (m)'); ylabel('North (m)');
%     title({'最终手势矢量重构图', '首尾无缝拼接 (保留原始方向和长度)'});
%     
%     plot(ax_f, 0, 0, '^', 'MarkerSize', 12, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'DisplayName', 'Receiver');
%     colors = lines(seg_cnt);
%     
%     % 画起点
%     start_p = final_draw_data(1).start;
%     plot(ax_f, start_p(1), start_p(2), 'kp', 'MarkerSize', 12, 'MarkerFaceColor', 'g', 'DisplayName', 'Start');
% 
%     for k = 1:seg_cnt
%         d = final_draw_data(k);
%         col = colors(k, :);
%         
%         plot(ax_f, [d.start(1), d.end(1)], [d.start(2), d.end(2)], ...
%                    '-', 'Color', col, 'LineWidth', 4, 'DisplayName', sprintf('Stroke #%d', k));
%         quiver(ax_f, d.start(1), d.start(2), d.vec(1), d.vec(2), ...
%                'Color', col, 'LineWidth', 2, 'MaxHeadSize', 0.4, 'AutoScale', 'off', 'HandleVisibility', 'off');
%         text(ax_f, d.start(1), d.start(2), sprintf('%d', k), ...
%             'Color', 'w', 'FontWeight', 'bold', 'FontSize', 10, 'BackgroundColor', col, 'HorizontalAlignment','center');
%     end
%     axis(ax_f, 'tight');
%     xl = xlim(ax_f); yl = ylim(ax_f);
%     xlim(ax_f, xl + [-0.5 0.5]); ylim(ax_f, yl + [-0.5 0.5]);
%     legend('Location', 'bestoutside');
% else
%     fprintf('⚠️ 没有检测到有效线段，跳过重构。\n');
% end
% 
% fprintf('✅ 所有分析完成。\n');

