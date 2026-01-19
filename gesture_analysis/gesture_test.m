% ============== test.m (最终调用脚本) ==============
clear;
clc; 
close all;
%%
obs_filepath = 'A_1_13_1.obs'; 
nav_filepath = '2026_1_13.nav'; 
% --- 2. 解析文件 ---
fprintf('--> 正在解析观测文件: %s\n', obs_filepath);
obs_data = parse_rinex_obs(obs_filepath);
fprintf('--> 正在解析导航文件: %s\n', nav_filepath);
nav_data = parse_rinex_nav_multi_gnss(nav_filepath);

% obs_data = generate_ideal_multi_shape(obs_data, nav_data, 'Star');

% fprintf('\n✅ 文件解析全部完成\n\n');
% calculate_and_plot_all_skyplot(obs_data, nav_data);

%%

%模拟GNSS欺骗
% obs_replay = simulate_gnss_spoofing(obs_data, nav_data, 'REPLAY');

%%

%第二步的轨迹推演需要这一步的参数，将obs_clean单独拎出来方便进行下一步特征归一化，
% obs_clean就是经过baseline处理过的os文件
plot_sn(obs_data);
[obs_clean, step1_res] = gesture_analysis_baseline_gvi(obs_data);
% plot_sn(obs_clean);
[obs_waveform, step1_res_shaped] = waveform_reshaping(obs_data, obs_clean, step1_res);
% plot_sn(obs_waveform);

%%

run_gesture_analysis_continuous_track(obs_waveform, nav_data, step1_res_shaped);
run_gesture_analysis_continuous_track_line(obs_waveform, nav_data, step1_res_shaped);
