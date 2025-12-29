% ============== test.m (最终调用脚本) ==============
clear;
clc; 
close all;
%%
obs_filepath = 'fingure_little_A_12_12_3.obs'; 
nav_filepath = 'arounds_12_12_1.nav'; 
% --- 2. 解析文件 ---
fprintf('--> 正在解析观测文件: %s\n', obs_filepath);
obs_data = parse_rinex_obs(obs_filepath);
fprintf('--> 正在解析导航文件: %s\n', nav_filepath);
nav_data = parse_rinex_nav_multi_gnss(nav_filepath);


obs_data = generate_ideal_multi_shape(obs_data, nav_data, 'Star');
% fprintf('\n✅ 文件解析全部完成\n\n');
% calculate_and_plot_all_skyplot(obs_data, nav_data);

%%
%旧方法

%手势切片
% [segments, volatility, t_grid, sats, PARA] = step1_segmentation_GVI(obs_data);

%切片手势识别方向
% step2_direction_estimation(obs_data, nav_data, segments, volatility, t_grid, sats);

%%

%可视化感知范围
% visualize_sensing_range(obs_data, nav_data);
%%

%模拟GNSS欺骗
% obs_replay = simulate_gnss_spoofing(obs_data, nav_data, 'REPLAY');

%%
%旧方法

%测试aheadN预处理算法
% [clean_snr, sat_ids, t_axis] = test_snr_flattening(obs_data, nav_data);

%aheadN预处理算法处理后利用切片识别方向算法，带手臂去除
% [draw_vecs1, segments1] = run_gesture_analysis_robust_aheadN(obs_data, nav_data);


% [analysis_results, segments] = run_gesture_analysis_3d_pipeline(obs_data, nav_data);
%%
%测试baseline预处理算法
% [clean_snr_matrix, sat_ids, t_axis] = test_snr_baseline_algorithm(obs_data, nav_data);

%baseline预处理算法处理后利用切片识别方向算法，带手臂去除
% [draw_vecs, segments] = run_gesture_analysis_robust_baseline(obs_data, nav_data);



%连续手势识别，重心聚类，会导致重心偏移 
%12.17修为仰角加权，极度抑制低仰角卫星


%第二步的轨迹推演需要这一步的参数，将obs_clean单独拎出来方便进行下一步特征归一化，
% obs_clean就是经过baseline处理过的os文件
[obs_clean, step1_res] = gesture_analysis_baseline_gvi(obs_data);

% run_gesture_analysis_continuous_track(obs_clean, nav_data, step1_res);
run_gesture_analysis_continuous_track_line(obs_clean, nav_data, step1_res);

    %还需要加入方向约束
run_gesture_analysis_boundary_trackV3(obs_clean, nav_data, step1_res);
run_gesture_analysis_boundary_track(obs_clean, nav_data, step1_res);