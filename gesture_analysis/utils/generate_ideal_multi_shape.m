% =========================================================================
% generate_ideal_multi_shape_v7_switch.m
% 功能: 通用仿真生成器 (支持 A, B, M, Star, D, L, N, X, Z)
% 特性: 
%   1. [Switch Case]: 集成所有字母坐标，通过 TARGET_LETTER 切换。
%   2. [Base V4]: 严格保留了 v4 版的手臂遮挡、长间隔和绘图逻辑。
% =========================================================================

clearvars -except obs_data nav_data;

% --- [用户设置] 请在这里修改要生成的形状 ---
TARGET_LETTER = 'Star';  % 可选: 'A', 'B', 'M', 'Star', 'D', 'L', 'N', 'X', 'Z'

% --- [Part 0] 环境与数据检查 ---
if ~exist('obs_data', 'var') || ~exist('nav_data', 'var')
    error('错误: 工作区缺少数据模板！请先加载一份真实的 obs_data 和 nav_data。');
end

addpath(genpath('calculate_clock_bias_and_positon'));
addpath(genpath('nav_parse'));

fprintf('--> 启动仿真 V7 (Switch版): 目标 [%s]...\n', TARGET_LETTER);

%% ================= [Part 1] 仿真参数设置 =================

% 1. 信号参数
SIM.baseline_db    = 45;       % [基线] 平稳时的信号强度 (dB)
SIM.drop_depth_db  = 15;       % [波动] 遮挡时的下降深度 (dB)
SIM.noise_sigma    = 0.02;     % [噪声] 微小抖动

% 2. 几何参数
SIM.gesture_height = 0.20;     % 手势平面高度 (米)
SIM.arm_width      = 0.12;     % [关键] 手臂有效遮挡宽度 (米)
SIM.body_pos       = [0.0, -1.0]; % [关键] 身体位置 (接收机正南方1米)

% 3. 形状定义 (Switch Case)
switch TARGET_LETTER
    case 'A'
        P1 = [-0.40, -0.50]; P2 = [ 0.00,  0.50]; P3 = [ 0.40, -0.50];
        P4 = [-0.20, -0.10]; P5 = [ 0.20, -0.10];
        STAGES = {
            P1, P2, 1.5, true;   % /
            P2, P2, 3.0, false;  % Interval
            P2, P3, 1.5, true;   % \
            P3, P4, 3.0, false;  % Interval
            P4, P5, 1.5, true    % -
        };

    case 'B'
        P1 = [-0.40, -1]; P2 = [-0.40,  1]; P3 = [ 1.5,  0.40]; 
        P4 = [-0.40,  0.00]; P5 = [ 1.5,  0.00]; P6 = [-0.40, -1]; 
        STAGES = {
            P1, P2, 1.5, true; P2, P2, 3.0, false;
            P2, P3, 1.5, true; P3, P3, 3.0, false;
            P3, P4, 3.0, true; P4, P4, 3.0, false;
            P4, P5, 1.5, true; P5, P5, 3.0, false;
            P5, P6, 1.5, true
        };

    case 'M'
        P1 = [-0.40, -0.40]; P2 = [-0.40,  0.40]; P3 = [ 0.00,  0.00];
        P4 = [ 0.40,  0.40]; P5 = [ 0.40, -0.40];
        STAGES = {
            P1, P2, 1.5, true;  P2, P2, 3.0, false;
            P2, P3, 1.5, true;  P3, P3, 3.0, false;
            P3, P4, 1.5, true;  P4, P4, 3.0, false;
            P4, P5, 1.5, true
        };

    case 'Star'
        P1 = [-0.30, -0.45]; P2 = [ 0.00,  0.55]; P3 = [ 0.30, -0.45];
        P4 = [-0.48,  0.15]; P5 = [ 0.48,  0.15];
        STAGES = {
            P1, P2, 1.5, true; P2, P2, 3.0, false;
            P2, P3, 1.5, true; P3, P3, 3.0, false;
            P3, P4, 1.5, true; P4, P4, 3.0, false;
            P4, P5, 1.5, true; P5, P5, 3.0, false;
            P5, P1, 1.5, true
        };

    case 'L'
        P_TL=[-0.3,0.4]; P_BL=[-0.3,-0.4]; P_BR=[0.3,-0.4];
        STAGES = {
            P_TL, P_BL, 1.5, true;  % |
            P_BL, P_BL, 3.0, false; % Interval
            P_BL, P_BR, 1.5, true   % _
        };

    case 'D'
        P1=[-0.3,0.5]; P2=[-0.3,-0.5]; P3=[0.0,0.5]; P4=[0.3,0.0]; P5=[0.0,-0.5];
        STAGES = {
            P1, P2, 1.5, true;   % |
            P2, P1, 3.0, false;  % Interval
            P1, P3, 1.0, true;   % -
            P3, P4, 1.0, true;   % )
            P4, P5, 1.0, true;   % )
            P5, P2, 1.0, true    % -
        };

    case 'X'
        P_TL=[-0.3,0.4]; P_TR=[0.3,0.4]; P_BL=[-0.3,-0.4]; P_BR=[0.3,-0.4];
        STAGES = {
            P_TL, P_BR, 1.5, true;  % \
            P_BR, P_TR, 3.0, false; % Interval
            P_TR, P_BL, 1.5, true   % /
        };

    case 'Z'
        P_TL=[-0.3,0.4]; P_TR=[0.3,0.4]; P_BL=[-0.3,-0.4]; P_BR=[0.7,-0.4];
        STAGES = {
            P_TL, P_TR, 1.5, true;  % -
            P_TR, P_TR, 3.0, false; % Interval
            P_TR, P_BL, 1.5, true;  % /
            P_BL, P_BL, 3.0, false; % Interval
            P_BL, P_BR, 1.5, true   % _
        };

    case 'N'
        P_BL=[-0.3,-0.4]; P_TL=[-0.3,0.4]; P_BR=[0.3,-0.4]; P_TR=[0.3,0.4];
        STAGES = {
            P_BL, P_TL, 1.5, true;  % |
            P_TL, P_TL, 3.0, false; % Interval
            P_TL, P_BR, 1.5, true;  % \
            P_BR, P_BR, 3.0, false; % Interval
            P_BR, P_TR, 1.5, true   % |
        };

    otherwise
        error('未定义的字母: %s', TARGET_LETTER);
end

% 计算总时长
total_sim_duration = 0;
for k=1:size(STAGES, 1), total_sim_duration = total_sim_duration + STAGES{k,3}; end

%% ================= [Part 2] 数据准备 =================
ideal_obs = obs_data;
raw_times = [ideal_obs.time];
num_samples = length(raw_times);

start_idx = round(num_samples * 0.3); 
sampling_rate = 25; 
end_idx = start_idx + round(total_sim_duration * sampling_rate);

if end_idx > num_samples
    error('原始数据太短，请更换更长的数据文件。');
end

% 提取有效卫星
all_sat_ids = {}; for i=1:min(100, length(ideal_obs)), if ~isempty(ideal_obs(i).data), all_sat_ids = [all_sat_ids, fieldnames(ideal_obs(i).data)']; end; end
unique_sat_ids = unique(all_sat_ids); valid_sats = {}; for i=1:length(unique_sat_ids), sid=unique_sat_ids{i}; if ismember(sid(1),['G','C','E','J']), valid_sats{end+1}=sid; end; end

% 计算接收机参考位置
fprintf('    计算接收机参考位置...\n');
rec_pos_acc = [0,0,0]; count = 0;
for t = start_idx : 10 : end_idx
    try [rp,~,~]=calculate_receiver_position(obs_data, nav_data, t); rec_pos_acc=rec_pos_acc+rp; count=count+1; catch, end
end
if count==0, error('无法计算接收机位置'); end
rec_pos_mean = rec_pos_acc / count;
[lat0, lon0, alt0] = ecef2geodetic(rec_pos_mean(1), rec_pos_mean(2), rec_pos_mean(3));

%% ================= [Part 3] 核心循环: 信号注入与轨迹记录 =================
fprintf('    正在注入信号并记录真值轨迹...\n');

% --- 用于记录 Ground Truth 的数组 ---
num_sim_pts = end_idx - start_idx + 1;
gt_trace_x = NaN(num_sim_pts, 1);
gt_trace_y = NaN(num_sim_pts, 1);
gt_pen_down = false(num_sim_pts, 1);

progr_step = floor(num_sim_pts/10);

for t_idx = 1 : num_samples
    
    if mod(t_idx-start_idx, progr_step) == 0 && t_idx >= start_idx && t_idx <= end_idx
        fprintf('    进度: %.0f%%\n', (t_idx-start_idx)/(end_idx-start_idx)*100);
    end
    
    % --- 1. 计算当前手的位置 & 笔触状态 ---
    current_hand_pos = [NaN, NaN]; % 默认无效
    is_pen_down = false;
    
    if t_idx >= start_idx && t_idx <= end_idx
        dt = (t_idx - start_idx) / sampling_rate;
        elapsed = 0;
        for k = 1:size(STAGES, 1)
            dur = STAGES{k,3};
            if dt <= (elapsed + dur)
                local_prog = max(0, min(1, (dt - elapsed) / dur));
                current_hand_pos = STAGES{k,1} + (STAGES{k,2} - STAGES{k,1}) * local_prog;
                is_pen_down = STAGES{k,4};
                
                % [记录 Ground Truth]
                trace_idx = t_idx - start_idx + 1;
                gt_trace_x(trace_idx) = current_hand_pos(1);
                gt_trace_y(trace_idx) = current_hand_pos(2);
                gt_pen_down(trace_idx) = is_pen_down;
                break;
            end
            elapsed = elapsed + dur;
        end
    end
    
    % --- 2. 遍历卫星 & 手臂遮挡计算 ---
    if isempty(ideal_obs(t_idx).data), continue; end
    try [~, ~, sat_states] = calculate_receiver_position(obs_data, nav_data, t_idx); catch, continue; end
    
    for s = 1:length(valid_sats)
        sid = valid_sats{s};
        if ~isfield(ideal_obs(t_idx).data, sid), continue; end
        
        sim_val = SIM.baseline_db + randn() * SIM.noise_sigma;
        
        % 仅在"落笔"且位置有效时计算遮挡
        if is_pen_down && ~isnan(current_hand_pos(1)) && isfield(sat_states, sid)
            sat_p = sat_states.(sid).position;
            [e, n, u] = ecef2enu(sat_p(1)-rec_pos_mean(1), sat_p(2)-rec_pos_mean(2), sat_p(3)-rec_pos_mean(3), lat0, lon0, alt0);
            
            if u > 0
                scale = SIM.gesture_height / u;
                P = [scale * e, scale * n]; % 投影点
                A = current_hand_pos; B = SIM.body_pos; % 手臂线段端点
                
                vec_AB = B - A; vec_AP = P - A;
                len_sq = sum(vec_AB.^2);
                if len_sq > 0
                    t_proj = max(0, min(1, dot(vec_AP, vec_AB) / len_sq));
                    dist_to_arm = norm(P - (A + t_proj * vec_AB));
                    if dist_to_arm < SIM.arm_width
                        sim_val = SIM.baseline_db - SIM.drop_depth_db + randn() * SIM.noise_sigma;
                    end
                end
            end
        end
        
        snr_struct = ideal_obs(t_idx).data.(sid).snr;
        fields = fieldnames(snr_struct);
        for f = 1:length(fields), ideal_obs(t_idx).data.(sid).snr.(fields{f}) = sim_val; end
    end
end
obs_data = ideal_obs; % 更新工作区

%% ================= [Part 4] Ground Truth 可视化 =================
fprintf('--> 生成 Ground Truth 轨迹图...\n');
figure('Name', sprintf('Ideal Ground Truth Trajectory [%s]', TARGET_LETTER), 'Position', [100, 100, 600, 600], 'Color', 'w');
ax = axes; hold(ax, 'on'); grid(ax, 'on'); axis(ax, 'equal');
xlabel('East (m)'); ylabel('North (m)');

% 1. 画接收机位置和身体参考点
plot(ax, 0, 0, '^', 'MarkerSize', 12, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'DisplayName', 'Receiver');
plot(ax, SIM.body_pos(1), SIM.body_pos(2), 'bs', 'MarkerSize', 10, 'MarkerFaceColor', 'b', 'DisplayName', 'Body/Shoulder Ref');

% 2. 画轨迹 (区分落笔和抬手)
plot_x_up = gt_trace_x; plot_y_up = gt_trace_y;
plot_x_up(gt_pen_down) = NaN; plot_y_up(gt_pen_down) = NaN;
plot_x_down = gt_trace_x; plot_y_down = gt_trace_y;
plot_x_down(~gt_pen_down) = NaN; plot_y_down(~gt_pen_down) = NaN;

% 绘制抬手移动 (灰色虚线)
plot(ax, plot_x_up, plot_y_up, 'k--', 'LineWidth', 1, 'Color', [0.6 0.6 0.6], 'DisplayName', 'Pen Up (Interval Move)');
% 绘制落笔书写 (蓝色实线 - 核心轨迹)
plot(ax, plot_x_down, plot_y_down, 'b-', 'LineWidth', 3, 'DisplayName', 'Pen Down (Target Stroke)');

% 3. 标记起点和终点
first_idx = find(~isnan(gt_trace_x), 1, 'first');
if ~isempty(first_idx)
     plot(ax, gt_trace_x(first_idx), gt_trace_y(first_idx), 'go', 'MarkerFaceColor', 'g', 'MarkerSize', 10, 'DisplayName', 'Start');
end
last_idx = find(~isnan(gt_trace_x), 1, 'last');
 if ~isempty(last_idx)
     plot(ax, gt_trace_x(last_idx), gt_trace_y(last_idx), 'rs', 'MarkerFaceColor', 'r', 'MarkerSize', 10, 'DisplayName', 'End');
end

title({sprintf('Ideal Ground Truth Trajectory (Letter %s)', TARGET_LETTER), 'Blue: Signal Drop Regions, Gray: Intervals'});
legend('Location', 'best');

% 自动缩放视角
valid_pts_x = gt_trace_x(~isnan(gt_trace_x));
valid_pts_y = gt_trace_y(~isnan(gt_trace_y));
if ~isempty(valid_pts_x)
    max_range = max([max(abs(valid_pts_x)), max(abs(valid_pts_y)), 0.6]);
else
    max_range = 0.6;
end
xlim([-max_range*1.2, max_range*1.2]);
ylim([-max_range*1.2, max_range*1.2] + SIM.body_pos(2)/2); 

fprintf('✅ 仿真与绘图全部完成。请运行分析脚本进行测试。\n');