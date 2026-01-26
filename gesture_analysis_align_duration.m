% =========================================================================
% gesture_analysis_align_duration.m (全周期版 V3)
% 功能: 强制对齐卫星动作时长 & 全周期覆盖
%
% 描述:
%   1. [对齐]: 针对 0/10 数据，检测每颗卫星的独立脉冲，将其时长强制对齐为 fixed_dur_pts。
%   2. [分段]: 放弃传统的起止时间判定，直接将分段设置为【整个数据周期】。
%      (依赖后续轨迹算法自动跳过值为0的区域)
%
% 输入:
%   obs_in          : 输入观测数据 (包含 0/10 的 SNR)
%   step1_res_in    : 输入结果结构体 (包含 volatility_matrix 和 segments)
%
% 输出:
%   obs_out         : 对齐及回填后的观测数据
%   step1_res_out   : 更新后的结果结构体 (segments 变为 1~End)
% =========================================================================

function [obs_out, step1_res_out] = gesture_analysis_align_duration(obs_in, step1_res_in)

    %% 1. 参数定义
    PARA.fixed_duration_pts = 80; % 强制对齐长度 (约 0.8秒)

    fprintf('--> [特征预处理] 启动强制时长对齐 (Fixed Dur: %d pts) & 全周期覆盖...\n', PARA.fixed_duration_pts);

    %% 2. 初始化
    obs_out = obs_in;
    step1_res_out = step1_res_in;
    
    mat = step1_res_in.volatility_matrix;
    % 依然需要原始分段来定位“去哪里找脉冲”，避免全图瞎找
    old_segments = step1_res_in.segments; 
    [num_samples, num_sats] = size(mat);
    
    % 创建新矩阵 (全0初始化)
    aligned_matrix = zeros(size(mat));

    %% 3. 遍历原始分段与卫星 (执行对齐)
    % 逻辑保持不变：在 Step 1 指示的大致范围内，精修每颗卫星的波形
    
    for k = 1:length(old_segments)
        seg_s = old_segments(k).start_idx;
        seg_e = old_segments(k).end_idx;
        
        for s = 1:num_sats
            snippet = mat(seg_s:seg_e, s);
            
            % 寻找所有独立的波动段 (Multi-Pulse)
            is_hit = snippet > 5;
            if ~any(is_hit), continue; end
            
            [labels, num_events] = bwlabel(is_hit);
            
            for ev = 1:num_events
                idx = find(labels == ev);
                local_start = idx(1);
                
                % 转换为全局坐标
                global_start = seg_s + local_start - 1;
                
                % 计算强制结束时间
                global_end = global_start + PARA.fixed_duration_pts - 1;
                global_end = min(global_end, num_samples);
                
                % 填入 10.0
                aligned_matrix(global_start:global_end, s) = 10.0;
            end
        end
    end

    % 更新矩阵
    step1_res_out.volatility_matrix = aligned_matrix;

    %% 4. [修改] 重置分段为全周期 (Full Cycle Segmentation)
    % 逻辑：不再计算哪里开始哪里结束，直接把大门全部打开。
    % 既然数据是干净的，轨迹算法读到 0 自然会跳过，读到 10 自然会计算。
    
    fprintf('    [Override] 正在将动作分段设置为整个数据周期 (1 - %d)...\n', num_samples);
    
    % 创建唯一的全长分段
    % id=1, start=1, end=num_samples
    new_segments = struct('id', 1, 'start_idx', 1, 'end_idx', num_samples);
    
    step1_res_out.segments = new_segments;

    %% 5. 数据回填 (Data Injection)
    fprintf('    正在将对齐后的数据回填至观测结构体...\n');
    
    t_grid = step1_res_in.t_grid;
    valid_sats = step1_res_in.valid_sats;
    sat_map = containers.Map(valid_sats, 1:num_sats);
    
    for k = 1:length(obs_out)
        [~, t_idx] = min(abs(t_grid - obs_out(k).time));
        epoch_sats = fieldnames(obs_out(k).data);
        for i = 1:length(epoch_sats)
            sid = epoch_sats{i};
            if isKey(sat_map, sid)
                col_idx = sat_map(sid);
                val_to_inject = aligned_matrix(t_idx, col_idx);
                
                snr_struct = obs_out(k).data.(sid).snr;
                fds = fieldnames(snr_struct);
                if ~isempty(fds)
                    obs_out(k).data.(sid).snr.(fds{1}) = val_to_inject;
                end
            end
        end
    end

    fprintf('✅ 处理完成：数据已强制对齐，分段已设为全周期。\n');
end