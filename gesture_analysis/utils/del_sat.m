% =========================================================================
% del_sat.m
% 功能: 从观测数据中剔除指定的卫星 (Blacklist Removal)
%
% 描述:
%   遍历观测数据结构体，检查并删除在 black_list 中指定的卫星数据。
%   用于去除那些数据混乱、无法滤波或已知有问题的卫星。
%
% 输入:
%   obs_in      : 原始观测数据结构体 (obs_data)
%   black_list  : 要删除的卫星ID列表 (Cell Array, e.g., {'G20', 'C01'})
%
% 输出:
%   obs_out     : 剔除指定卫星后的观测数据
% =========================================================================

function obs_out = del_sat(obs_in, black_list)

    fprintf('--> [数据清洗] 正在剔除黑名单卫星...\n');
    
    % 1. 输入检查
    if ischar(black_list)
        % 如果用户只输了一个字符串 'G01'，转为 cell
        black_list = {black_list};
    elseif isstring(black_list)
        % 如果用户输的是 string 数组，转为 cell
        black_list = cellstr(black_list);
    end
    
    if isempty(black_list)
        fprintf('    警告: 黑名单为空，未做任何修改。\n');
        obs_out = obs_in;
        return;
    end

    % 2. 打印要剔除的列表
    fprintf('    目标剔除: %s\n', strjoin(black_list, ', '));

    % 3. 遍历每一历元进行剔除
    obs_out = obs_in;
    total_epochs = length(obs_out);
    removed_count = 0; % 统计删除了多少次(调试用)

    for k = 1:total_epochs
        % 获取当前历元下所有的卫星ID (Fields)
        current_sats = fieldnames(obs_out(k).data);
        
        % 找出当前历元中，既存在于卫星列表又存在于黑名单中的卫星 (交集)
        % intersect 返回的是我们要删除的字段名
        sats_to_remove = intersect(current_sats, black_list);
        
        % 如果交集不为空，说明这就需要删
        if ~isempty(sats_to_remove)
            % rmfield 可以一次性删除多个字段，效率较高
            obs_out(k).data = rmfield(obs_out(k).data, sats_to_remove);
            removed_count = removed_count + length(sats_to_remove);
        end
    end

    fprintf('✅ 剔除完成。共移除了 %d 个卫星数据点。\n', removed_count);

end