% =========================================================================
% waveform_synchronization.m
% 功能: 多星动作时域强制对齐 (Time-Domain Synchronization)
% 描述:
%   解决由于卫星几何分布导致的手势动作持续时间不一致的问题。
%   (例如: 卫星A跌落4秒，卫星B跌落2秒，会导致后续聚点不稳定)
%
%   [核心逻辑]:
%   1. 基于 GVI (Global Volatility Index) 确定手势的"最大时空包络"。
%   2. 提取满足最小持续时间 (min_action_pts) 的全局动作窗口。
%   3. 检查每颗卫星在该窗口内是否"参与" (有过任何波动)。
%   4. 对参与的卫星，将其信号在该窗口内强制拉伸填满，实现起止时间严格对齐。
%
% [调用格式]:
%   [obs_sync, res_sync] = waveform_synchronization(obs_waveform, step1_res_shaped);
%
% [输入参数]:
%   1. obs_waveform (struct): 
%      上一级 (跌落检测) 输出的波形数据，包含初步整形的方波。
%   2. step1_res_shaped (struct): 
%      上一级输出的结果包，必须包含 .volatility_matrix (矩阵形式的波形)。
%
% [返回值]:
%   1. obs_sync (struct): 
%      同步对齐后的观测数据 (方波已拉伸)。
%   2. res_sync (struct): 
%      同步对齐后的结果包 (矩阵已更新)。
% =========================================================================

function [obs_sync, res_sync] = waveform_synchronization(obs_waveform, step1_res_shaped)

fprintf('--> [Sync] 启动多星动作时域同步 (Waveform Synchronization)...\n');

%% 1. 初始化
obs_sync = obs_waveform;
res_sync = step1_res_shaped;

% 提取核心矩阵 (此时应该是 0/10 的方波矩阵)
vol_mat   = step1_res_shaped.volatility_matrix;
t_grid    = step1_res_shaped.t_grid;
valid_sats = step1_res_shaped.valid_sats;
[num_samples, num_sats] = size(vol_mat);

%% 2. 同步参数 (Synchronization Parameters)
% [参数] 最小动作持续时间 (采样点数)
% 定义一个有效的"全局动作"至少应该持续多久。
% 小于此长度的 GVI 波动将被视为噪声忽略。
SYNC.min_action_pts = 20;  

% [参数] 填充值 (方波高电平)
SYNC.fill_value     = 10.0;

fprintf('    同步策略: Global Window Min Duration = %d pts\n', SYNC.min_action_pts);

%% 3. 计算 GVI 并提取全局窗口
% GVI: 所有卫星波形的叠加，代表了动作的"总包络"
gvi_curve = sum(vol_mat, 2);

% 二值化 GVI (只要有任意一颗星在动，GVI > 0)
is_active = gvi_curve > 0;

% 连通域分析，提取动作片段
[label_mat, num_events] = bwlabel(is_active);

fprintf('    检测到 %d 个潜在动作片段 (GVI Segments)\n', num_events);

%% 4. 逐窗口执行强制对齐
sync_mat = vol_mat; % 复制一份用于修改

for k = 1:num_events
    % 获取当前窗口的索引
    idx = find(label_mat == k);
    
    % --- 过滤: 忽略过短的碎片 ---
    if length(idx) < SYNC.min_action_pts
        continue; 
    end
    
    % 获取窗口的时间范围 (用于日志)
    t_start_idx = idx(1);
    t_end_idx   = idx(end);
    
    % --- 遍历每一颗卫星，检查参与性 ---
    for s = 1:num_sats
        % 提取该卫星在当前窗口内的数据
        sat_segment = vol_mat(idx, s);
        
        % [核心逻辑] 参与性检测
        % 只要该卫星在这个时间段内有过"动静" (最大值 > 0)，
        % 就认为它属于这个动作事件，必须强制对齐。
        if max(sat_segment) > 0
            % 强制拉伸: 将该窗口内的所有点填满
            sync_mat(idx, s) = SYNC.fill_value;
        end
    end
end

% 更新结果矩阵
res_sync.volatility_matrix = sync_mat;

%% 5. 数据回填 (Data Injection)
fprintf('--> [Injection] 正在回填同步后的波形数据...\n');

sat_map = containers.Map(valid_sats, 1:num_sats);

for k = 1:length(obs_sync)
    t_now = obs_sync(k).time;
    % 寻找最近的时间网格点
    [~, t_idx] = min(abs(t_grid - t_now));
    
    epoch_sats = fieldnames(obs_sync(k).data);
    for i = 1:length(epoch_sats)
        sid = epoch_sats{i};
        if isKey(sat_map, sid)
            col_idx = sat_map(sid);
            
            % 取出同步后的值
            val_sync = sync_mat(t_idx, col_idx);
            
            % 覆盖 SNR 数据
            snr_struct = obs_sync(k).data.(sid).snr;
            fds = fieldnames(snr_struct);
            if ~isempty(fds)
                target_code = fds{1};
                obs_sync(k).data.(sid).snr.(target_code) = val_sync;
            end
        end
    end
end

fprintf('✅ 同步完成。所有参与卫星的动作起止时间已强制对齐。\n');

end