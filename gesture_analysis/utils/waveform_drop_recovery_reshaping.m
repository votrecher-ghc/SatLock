% % =========================================================================
% % waveform_drop_recovery_parametric.m
% % 功能: 参数化跌落检测 (Parametric Drop Detection)
% % 描述:
% %   基于“跌落-恢复”模型检测手势。
% %   严格执行用户指令：
% %   1. 动作持续时间：只检查最小持续时间 (min_dur_pts)。
% %   2. 深度阈值：只要在连通域内【有一个采样点】深度超过阈值即可。
% %
% % [调用格式]:
% %   [obs_waveform_param, step1_res_shaped] = waveform_drop_recovery_parametric(obs_clean, step1_res);
% %
% % [核心参数]:
% %   ALG.min_dur_pts:  最小持续点数 (默认 5, 约0.2s)
% %   ALG.min_depth_db: 最小深度阈值 (默认 2.0dB)
% % =========================================================================
% 
% function [obs_waveform_param, step1_res_shaped] = waveform_drop_recovery_parametric(obs_clean, step1_res)
% 
% fprintf('--> [Event] 启动参数化跌落检测 (Parametric Drop Detection)...\n');
% 
% %% 1. 初始化
% obs_waveform_param = obs_clean; 
% step1_res_shaped = step1_res;
% 
% cn0_mat = step1_res.cn0_clean_matrix; 
% t_grid  = step1_res.t_grid;
% valid_sats = step1_res.valid_sats;
% [num_samples, num_sats] = size(cn0_mat);
% 
% %% 2. 核心过滤参数 (Filter Parameters)
% % [参数1] 最小持续时间 (采样点数)
% % 只要动作持续超过这个点数，就认为时长达标
% ALG.min_dur_pts   = 20;      
% 
% % [参数2] 最小深度阈值 (dB)
% % 只要动作过程中有任意一个点跌落超过此值，就认为深度达标
% ALG.min_depth_db  = 2.0;    
% 
% % 输出电平
% ALG.output_weight = 10.0;
% 
% fprintf('    策略: Duration >= %d pts, Any Point Depth > %.1fdB\n', ...
%     ALG.min_dur_pts, ALG.min_depth_db);
% 
% %% 3. 逐卫星处理
% shaped_matrix = zeros(size(cn0_mat));
% 
% for s = 1:num_sats
%     raw_col = cn0_mat(:, s);
%     if sum(~isnan(raw_col)) < 10, continue; end
%     
%     % --- Step A: 确定基准线 (Baseline) ---
%     % 使用全局众数 (Mode) 作为最稳健的“原 sn0”
%     valid_data = raw_col(~isnan(raw_col));
%     baseline = mode(round(valid_data)); 
%     
%     % 扩展为向量
%     baseline_curve = repmat(baseline, num_samples, 1);
%     
%     % --- Step B: 提取潜在区间 (Connectivity) ---
%     % 定义：凡是低于基准线的，都算作“跌落过程”
%     % 这里不加任何缓冲，严格执行 < Baseline
%     is_below = raw_col < baseline_curve;
%     is_below(isnan(is_below)) = 0;
%     
%     % 标记连通域 (Labeling)
%     [label_mat, num_events] = bwlabel(is_below);
%     
%     final_mask = false(num_samples, 1);
%     
%     % --- Step C: 参数过滤 (Parametric Filtering) ---
%     for k = 1:num_events
%         idx = find(label_mat == k);
%         if isempty(idx), continue; end
%         
%         % [条件1] 动作持续时间检查
%         % 只需要最小持续时间
%         if length(idx) < ALG.min_dur_pts
%             continue; % 太短，丢弃
%         end
%         
%         % [条件2] 深度阈值检查
%         % 只要在这一段时间内【有一个采样点】深度超过了阈值就可以
%         current_vals = raw_col(idx);
%         current_base = baseline_curve(idx);
%         drops = current_base - current_vals;
%         
%         if max(drops) <= ALG.min_depth_db
%             continue; % 全程都很浅，没有一个点超过阈值，丢弃
%         end
%         
%         % --- 满足所有条件，保留 ---
%         % 稍微向外扩展1个点，保证包住边缘
%         s_idx = max(1, idx(1)-1);
%         e_idx = min(num_samples, idx(end)+1);
%         final_mask(s_idx:e_idx) = true;
%     end
%     
%     % 赋值
%     shaped_matrix(:, s) = double(final_mask) * ALG.output_weight;
% end
% 
% % 更新结果
% step1_res_shaped.volatility_matrix = shaped_matrix;
% 
% %% 4. 数据回填
% fprintf('--> [Injection] 回填参数化检测数据...\n');
% sat_map = containers.Map(valid_sats, 1:num_sats);
% 
% for k = 1:length(obs_waveform_param)
%     t_now = obs_waveform_param(k).time;
%     [~, t_idx] = min(abs(t_grid - t_now));
%     
%     epoch_sats = fieldnames(obs_waveform_param(k).data);
%     for i = 1:length(epoch_sats)
%         sid = epoch_sats{i};
%         if isKey(sat_map, sid)
%             col_idx = sat_map(sid);
%             val_to_inject = shaped_matrix(t_idx, col_idx);
%             
%             snr_struct = obs_waveform_param(k).data.(sid).snr;
%             fds = fieldnames(snr_struct);
%             if ~isempty(fds)
%                 target_code = fds{1};
%                 obs_waveform_param(k).data.(sid).snr.(target_code) = val_to_inject;
%             end
%         end
%     end
% end
% 
% fprintf('✅ 参数化跌落检测完成。\n');
% 
% end














% =========================================================================
% waveform_drop_recovery_reshaping.m (修改版)
% 功能: 参数化跌落检测 (Parametric Drop Detection) & 并集重分段
% 描述:
%   1. 基于“跌落-恢复”模型对每颗卫星进行单星去噪和整形 (0/10)。
%   2. [新增] 基于整形后数据的"并集" (Union) 重新计算动作起止时间，
%      确保轨迹推演覆盖动作的全程 (从第一颗星进入到最后一颗星离开)。
%
% [调用格式]:
%   [obs_waveform_param, step1_res_shaped] = waveform_drop_recovery_parametric(obs_clean, step1_res);
%
% [核心参数]:
%   ALG.min_dur_pts:  最小持续点数 (默认 20)
%   ALG.min_depth_db: 最小深度阈值 (默认 2.0dB)
% =========================================================================

function [obs_waveform_param, step1_res_shaped] = waveform_drop_recovery_reshaping(obs_clean, step1_res)

fprintf('--> [Event] 启动参数化跌落检测 (Parametric Drop Detection)...\n');

%% 1. 初始化
obs_waveform_param = obs_clean; 
step1_res_shaped = step1_res;

cn0_mat = step1_res.cn0_clean_matrix; 
t_grid  = step1_res.t_grid;
valid_sats = step1_res.valid_sats;
[num_samples, num_sats] = size(cn0_mat);

%% 2. 核心过滤参数 (Filter Parameters)
% [参数1] 最小持续时间 (采样点数)
% 只要动作持续超过这个点数，就认为时长达标
ALG.min_dur_pts   = 20;      

% [参数2] 最小深度阈值 (dB)
% 只要动作过程中有任意一个点跌落超过此值，就认为深度达标
ALG.min_depth_db  = 2.0;    

% 输出电平
ALG.output_weight = 10.0;

fprintf('    策略: Duration >= %d pts, Any Point Depth > %.1fdB\n', ...
    ALG.min_dur_pts, ALG.min_depth_db);

%% 3. 逐卫星处理 (生成整形矩阵)
shaped_matrix = zeros(size(cn0_mat));

for s = 1:num_sats
    raw_col = cn0_mat(:, s);
    if sum(~isnan(raw_col)) < 10, continue; end
    
    % --- Step A: 确定基准线 (Baseline) ---
    % 使用全局众数 (Mode) 作为最稳健的“原 sn0”
    valid_data = raw_col(~isnan(raw_col));
    baseline = mode(round(valid_data)); 
    
    % 扩展为向量
    baseline_curve = repmat(baseline, num_samples, 1);
    
    % --- Step B: 提取潜在区间 (Connectivity) ---
    % 定义：凡是低于基准线的，都算作“跌落过程”
    % 这里不加任何缓冲，严格执行 < Baseline
    is_below = raw_col < baseline_curve;
    is_below(isnan(is_below)) = 0;
    
    % 标记连通域 (Labeling)
    [label_mat, num_events] = bwlabel(is_below);
    
    final_mask = false(num_samples, 1);
    
    % --- Step C: 参数过滤 (Parametric Filtering) ---
    for k = 1:num_events
        idx = find(label_mat == k);
        if isempty(idx), continue; end
        
        % [条件1] 动作持续时间检查
        if length(idx) < ALG.min_dur_pts
            continue; % 太短，丢弃
        end
        
        % [条件2] 深度阈值检查
        % 只要在这一段时间内【有一个采样点】深度超过了阈值就可以
        current_vals = raw_col(idx);
        current_base = baseline_curve(idx);
        drops = current_base - current_vals;
        
        if max(drops) <= ALG.min_depth_db
            continue; % 全程都很浅，没有一个点超过阈值，丢弃
        end
        
        % --- 满足所有条件，保留 ---
        % 稍微向外扩展1个点，保证包住边缘
        s_idx = max(1, idx(1)-1);
        e_idx = min(num_samples, idx(end)+1);
        final_mask(s_idx:e_idx) = true;
    end
    
    % 赋值
    shaped_matrix(:, s) = double(final_mask) * ALG.output_weight;
end

% 保存整形矩阵
step1_res_shaped.volatility_matrix = shaped_matrix;

%% 4. [新增] 核心修改：基于并集的重分段 (Union Segmentation)
% 目的：根据清洗后的干净数据，重新定义动作的物理起止时间。
% 逻辑：只要有任意一颗卫星是 Hit (10)，就算动作开始；全都是 Miss (0)，才算结束。

fprintf('    [Refinement] 正在基于整形数据执行并集重分段 (Union Strategy)...\n');

% 4.1 计算全局激活掩膜 (Union Mask)
% 对每一行求最大值，只要 > 0 即为 active
is_global_active = max(shaped_matrix, [], 2) > 0;

% 4.2 提取上升沿和下降沿
% 前后补0以处理边界情况
padded_active = [0; is_global_active; 0];
edges = diff(padded_active);
s_indices = find(edges == 1);       % 上升沿 (Start)
e_indices = find(edges == -1) - 1;  % 下降沿 (End)

% 4.3 构建新的分段结构体
num_new_segs = length(s_indices);
new_segments = struct('id', {}, 'start_idx', {}, 'end_idx', {});

for i = 1:num_new_segs
    new_segments(i).id = i;
    new_segments(i).start_idx = s_indices(i);
    new_segments(i).end_idx   = e_indices(i);
    
    % 计算持续时间 (秒)
    if ~isempty(t_grid)
        dur_sec = seconds(t_grid(e_indices(i)) - t_grid(s_indices(i)));
    else
        dur_sec = NaN;
    end
    
    % (可选) 可以在这里补充 peak_time 等信息，但对于轨迹推演不是必须的
end

% 4.4 覆盖输出结果
% 注意：我们修改的是 step1_res_shaped 里的 segments，step1_res 保持不变
step1_res_shaped.segments = new_segments;

fprintf('    ✅ 分段已更新: 从 Step1 的 %d 段优化为 %d 段 (基于物理并集)。\n', ...
    length(step1_res.segments), num_new_segs);

%% 5. 数据回填 (Data Injection)
fprintf('--> [Injection] 回填参数化检测数据...\n');
sat_map = containers.Map(valid_sats, 1:num_sats);

for k = 1:length(obs_waveform_param)
    t_now = obs_waveform_param(k).time;
    % 寻找最近的时间点
    [~, t_idx] = min(abs(t_grid - t_now));
    
    epoch_sats = fieldnames(obs_waveform_param(k).data);
    for i = 1:length(epoch_sats)
        sid = epoch_sats{i};
        if isKey(sat_map, sid)
            col_idx = sat_map(sid);
            val_to_inject = shaped_matrix(t_idx, col_idx);
            
            % 注入到 SNR 字段
            snr_struct = obs_waveform_param(k).data.(sid).snr;
            fds = fieldnames(snr_struct);
            if ~isempty(fds)
                target_code = fds{1};
                obs_waveform_param(k).data.(sid).snr.(target_code) = val_to_inject;
            end
        end
    end
end

fprintf('✅ 参数化跌落检测及分段优化完成。\n');

end

