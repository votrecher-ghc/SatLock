% =========================================================================
% gesture_test.m
% 两层式流程（保持原调用结构思路）：
%   第一层：数据处理（解析 + 仿真注入 + 预处理）
%   第二层：轨迹推演（纯数据驱动）
%
% 说明：
%   1) 第二层推演不使用 groundtruth，不使用目标字母；
%   2) groundtruth 仅用于“仿真验证对比绘图”；
%   3) 默认感知范围：左右 50cm、上下 50cm。
% =========================================================================

clear;
clc;
close all;

repo_dir = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(repo_dir));

%% ================= [用户配置] =================
obs_filepath = fullfile(repo_dir, 'data', '1_8', 'A_1_8_1.obs');
nav_filepath = fullfile(repo_dir, 'data', '1_8', '2026_1_8.nav');
target_letter = 'Star';            % 仅用于第一层仿真注入与对比绘图

% 感知范围配置（总跨度）
span_cfg = struct();
span_cfg.max_span_x = 0.50;     % 左右总跨度 50cm
span_cfg.max_span_y = 0.50;     % 上下总跨度 50cm

% 按你的最新要求：不做量化对比，但旧算法图要全部画出来
run_legacy_baselines = true;
show_gt_compare_plot = true;

%% ================= [第一层] 数据处理 =================
sim_cfg = struct();
sim_cfg.enable = true;
sim_cfg.target_letter = target_letter;
sim_cfg.max_span_x = span_cfg.max_span_x;
sim_cfg.max_span_y = span_cfg.max_span_y;

[obs_data, nav_data, ~] = prepare_input_data(obs_filepath, nav_filepath, sim_cfg);
[obs_clean, step1_res, obs_waveform, step1_res_shaped, obs_aligned, step1_aligned] = ...
    run_preprocess_pipeline(obs_data); %#ok<ASGLU>

% 结果输出目录（保存每种算法的 GT 对比图）
stamp = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
out_dir = fullfile(repo_dir, 'gesture_analysis', 'results', ['gesture_test_compare_', stamp]);
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

alg_results = struct('name', {}, 'x', {}, 'y', {}, 't', {}, 'conf', {});

%% ================= [可选] 原有算法链路 =================
if run_legacy_baselines
    % 以下为你原有算法链路，仅用于可视化查看效果，不参与对比打分
    [x1, y1, t1, e1] = run_gesture_analysis_continuous_track(obs_waveform, nav_data, step1_res_shaped);
    alg_results(end+1) = struct('name', 'continuous_track', 'x', x1, 'y', y1, 't', t1, 'conf', e1); %#ok<SAGROW>

    [x2, y2, t2, e2] = run_gesture_analysis_continuous_track_line(obs_aligned, nav_data, step1_aligned);
    alg_results(end+1) = struct('name', 'continuous_track_line', 'x', x2, 'y', y2, 't', t2, 'conf', e2); %#ok<SAGROW>

    [x3, y3, t3, e3] = run_gesture_analysis_boundary_trackV3(obs_waveform, nav_data, step1_res_shaped);
    alg_results(end+1) = struct('name', 'boundary_track_v3', 'x', x3, 'y', y3, 't', t3, 'conf', e3); %#ok<SAGROW>

    final_traj = run_gesture_analysis_physics_engine(obs_waveform, nav_data, step1_res_shaped);
    if isstruct(final_traj) && isfield(final_traj, 'x') && isfield(final_traj, 'y')
        if isfield(final_traj, 't')
            t4 = final_traj.t;
        else
            t4 = [];
        end
        alg_results(end+1) = struct('name', 'physics_engine', ...
            'x', final_traj.x, 'y', final_traj.y, 't', t4, 'conf', []); %#ok<SAGROW>
    end

    % 保留原 inverse-beam 图，方便观察你之前改动过的版本效果
    cfg_old_inv = struct();
    cfg_old_inv.debug = struct('verbose', true, 'plot', true);
    cfg_old_inv.track = struct('lambda_smooth', 10.0, 'final_smooth_pts', 3, 'max_jump_m', 0.22);
    [x5, y5, t5, c5, ~] = run_gesture_analysis_inverse_beam(obs_waveform, nav_data, step1_res_shaped, cfg_old_inv);
    alg_results(end+1) = struct('name', 'inverse_beam_legacy', 'x', x5, 'y', y5, 't', t5, 'conf', c5); %#ok<SAGROW>
end

%% ================= [第二层] 新算法轨迹推演 =================
cfg_l2 = struct();
cfg_l2.debug = struct('verbose', true, 'plot', true); % 直接弹窗画图
cfg_l2.model = struct('max_hand_radius', 0.40);
cfg_l2.grid = struct('x_min', -0.35, 'x_max', 0.35, 'y_min', -0.35, 'y_max', 0.35, 'step', 0.015);
cfg_l2.track = struct( ...
    'lambda_smooth', 12.0, ...
    'final_smooth_pts', 3, ...
    'max_jump_m', 0.18, ...
    'enforce_piecewise_linear', false, ...
    'template_snap_enable', false);

[traj_x, traj_y, traj_t, traj_conf, dbg_l2] = ...
    run_gesture_analysis_data_driven(obs_waveform, nav_data, step1_res_shaped, cfg_l2);
alg_results(end+1) = struct('name', 'data_driven_new', ...
    'x', traj_x, 'y', traj_y, 't', traj_t, 'conf', traj_conf); %#ok<SAGROW>

fprintf('\n[第二层] 新算法推演完成。\n');
fprintf('  输出点数: %d\n', numel(traj_x));
if isstruct(dbg_l2) && isfield(dbg_l2, 'uses_groundtruth')
    fprintf('  推演是否使用GT: %d\n', dbg_l2.uses_groundtruth);
end

%% ================= [对比绘图] 仅 GT vs 新算法 =================
if show_gt_compare_plot
    plot_all_algorithms_vs_gt( ...
        alg_results, step1_res_shaped.t_grid, target_letter, span_cfg, out_dir, true);
    fprintf('已保存各算法对比图到: %s\n', out_dir);
end
