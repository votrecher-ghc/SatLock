% ============== test.m (最终调用脚本) ==============
clear;
clc; 
close all;
%%
obs_filepath = 'fingure_C_12_8_20_1.obs'; 
nav_filepath = 'arounds_12_8_1.nav'; 
% --- 2. 解析文件 ---
fprintf('--> 正在解析观测文件: %s\n', obs_filepath);
obs_data = parse_rinex_obs(obs_filepath);
fprintf('--> 正在解析导航文件: %s\n', nav_filepath);
nav_data = parse_rinex_nav_multi_gnss(nav_filepath);

obs_data = generate_ideal_multi_shape(obs_data, nav_data, 'Star');


% remove_set = {'J08','S28','C49','S29','S37','R14','R15','R17','R23','R24'}; 
% obs_data = del_sat(obs_data, remove_set);

% fprintf('\n✅ 文件解析全部完成\n\n');
calculate_and_plot_all_skyplot(obs_data, nav_data);

%模拟GNSS欺骗
% obs_replay = simulate_gnss_spoofing(obs_data, nav_data, 'REPLAY');

%%
% plot_sn(obs_data);

[obs_clean, step1_res] = gesture_analysis_baseline_gvi(obs_data);
% plot_sn(obs_clean);

[obs_waveform, step1_res_shaped] = waveform_drop_recovery_reshaping(obs_clean, step1_res);
% plot_sn(obs_waveform);

[obs_aligned, step1_aligned] = gesture_analysis_align_duration(obs_waveform, step1_res_shaped);
% plot_sn(obs_aligned);


%%

run_gesture_analysis_continuous_track(obs_waveform, nav_data, step1_res_shaped);
run_gesture_analysis_continuous_track_line(obs_aligned, nav_data, step1_aligned);
run_gesture_analysis_boundary_trackV3(obs_waveform, nav_data, step1_res_shaped);

%%
% run_gesture_analysis_temporal_trajectory(obs_waveform, nav_data);