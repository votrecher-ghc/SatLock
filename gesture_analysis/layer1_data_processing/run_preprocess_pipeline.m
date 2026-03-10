function [obs_clean, step1_res, obs_waveform, step1_res_shaped, obs_aligned, step1_aligned] = ...
    run_preprocess_pipeline(obs_data)
% RUN_PREPROCESS_PIPELINE
% 第一层：数据处理流水线
% 包含：GVI 分段 -> 波形整形 -> 时长对齐

fprintf('--> [第一层] 执行 GVI 基础预处理...\n');
[obs_clean, step1_res] = gesture_analysis_baseline_gvi(obs_data);

fprintf('--> [第一层] 执行波形整形...\n');
[obs_waveform, step1_res_shaped] = waveform_drop_recovery_reshaping(obs_clean, step1_res);

fprintf('--> [第一层] 执行时长对齐...\n');
[obs_aligned, step1_aligned] = gesture_analysis_align_duration(obs_waveform, step1_res_shaped);
end

