# SatLock 会话交接文档

更新时间：2026-05-21 09:15:12 +08:00  
适用目的：给新会话中的大模型快速接手当前工程状态。  
重要说明：如果本文件与旧的口头上下文或部分历史实现细节冲突，优先以**本文件描述的当前工程状态**为准。

---

## 0. 2026-05-21 快速接手摘要

当前仓库仍是 **SatLock / StarDial** 项目，不是 `D:\EdgeDownload\PROJECT_CONTEXT.md` 中描述的 LubanCat-3 / RK3576 嵌入式项目。该外部文件内容与本仓库技术主线不一致，本次仅记录这个事实，不把 LubanCat 内容并入 SatLock 主线。

当前项目定位：

- SatLock 是基于 GNSS C/N0 或 SNR 近场衍射扰动的手势感知、轨迹恢复、认证与安全分析研究原型。
- 论文初稿为 `D:\Matproject\SatLock\StarDial.pdf`，已抽取文本 `D:\Matproject\SatLock\StarDial_extracted.txt`。
- 系统设计理论示意图单独放在 `D:\Matproject\SatLock\system_polt`，该目录与原始数据工作流解耦，主要服务论文 `III. SYSTEM DESIGN`。
- 当前最重要的代码出口仍是 `D:\Matproject\SatLock\gesture_analysis\scientific_graphing\export_paper_figures_data_driven.m`。
- 最新已验证论文图导出目录为：
  `D:\Matproject\SatLock\gesture_analysis\results\paper_figures_data_driven_output\paper_figures_data_driven_20260520_094905`

最近一轮主要工作集中在论文图格式精修：

- 全局论文图字体统一为 Arial。
- 多数主论文图采用大字号风格，重点图要求坐标轴外框完整、坐标轴线宽约 1.0、刻度朝内。
- `authentication_metrics_bar` 的纵轴数字已改为手绘小号刻度，最近固定到 40 号；柱状图上方不显示数值。
- `feature_space_pca` 和 `feature_space_tsne` 的散点面积已放大到至少旧版 2 倍，并保持类别中心不乱移。
- `traj_gallery_data_driven` 已改为手动 axes 网格，列间距按固定子图宽度后压缩，行间距恢复，避免上下行误压。
- `B1/B2/B3/C1/C2` 等系统设计图此前已做过大量版式调整，PDF/PNG 白边、尺寸一致性和 Arial 字体是核心约束。

本文件的后续章节保留了 2026-04-17 的架构背景与设计决策；若与本节最新图形状态冲突，以本节和最新导出目录为准。

---

## 1. 这份文档的作用

这不是完整项目说明书，完整背景请看：

- `D:\Matproject\SatLock\简介.md`

本文件的目标是：

- 总结当前已经完成的重构与绘图工作；
- 说明用户已经明确确认过的规则；
- 标出当前应该继续修改的关键文件和位置；
- 给出最新可用结果目录、常用命令、注意事项；
- 让新会话模型读完后可以直接继续工作，不需要重新从长上下文里摸索。

---

## 2. 当前工程的总体定位

这是一个基于 GNSS 信号扰动进行手势感知、轨迹恢复、手势认证与安全分析的 MATLAB 项目。

当前项目的主线已经被整理为：

1. 数据解析
2. 数据注入
3. 轨迹恢复
4. 手势认证
5. 科研绘图

其中：

- 第 3 层有“真正盲模板”的 `core_cases`；
- 同时保留了用于论文图和后续评估的 `gallery_cases`；
- **所有科研绘图、认证、论文图中的轨迹对比，当前都以 `trajectory.gallery_cases` 为主**；
- 用户接受这种分层方式，因为它保证了“严格分层”和“现有论文图效果”可以同时存在。

---

## 3. 当前目录结构与主入口

工程根目录：

- `D:\Matproject\SatLock`

当前核心目录：

- `D:\Matproject\SatLock\gesture_analysis\workflow`
- `D:\Matproject\SatLock\gesture_analysis\authentication`
- `D:\Matproject\SatLock\gesture_analysis\templates`
- `D:\Matproject\SatLock\gesture_analysis\scientific_graphing`
- `D:\Matproject\SatLock\gesture_analysis\continue`
- `D:\Matproject\SatLock\gesture_analysis\preprocess_feature_extraction`
- `D:\Matproject\SatLock\gesture_analysis\utils`

当前主入口脚本：

- `D:\Matproject\SatLock\gesture_analysis\gesture_test.m`

这个入口脚本的职责已经被明确为：

1. 解析原始数据
2. 注入模板数据
3. 可选地模拟攻击
4. 运行 Data-Driven 轨迹恢复
5. 基于恢复轨迹做 Score_k 手势认证
6. 直接导出论文用科研绘图

并且：

- `gesture_test.m` **不会再导出非论文图片**
- `trash` 文件夹被明确排除在 MATLAB path 之外

---

## 4. 当前分层实现文件

### 4.1 工作流装配器

- `D:\Matproject\SatLock\gesture_analysis\workflow\build_gesture_test_source.m`

它会把整个工作流打包为 `src` 结构，主要字段：

- `src.parsed`
- `src.injected`
- `src.scenario`
- `src.attack`
- `src.trajectory`
- `src.auth`
- `src.summary_tbl`
- `src.cfg`
- `src.cases`

### 4.2 四层工作流

- 第 1 层：`D:\Matproject\SatLock\gesture_analysis\workflow\layer1_parse_raw_data.m`
- 第 2 层：`D:\Matproject\SatLock\gesture_analysis\workflow\layer2_inject_templates.m`
- 场景模拟层：`D:\Matproject\SatLock\gesture_analysis\workflow\layer_scenario_simulation.m`
- 攻击模拟层：`D:\Matproject\SatLock\gesture_analysis\workflow\layer_attack_simulation.m`
- 第 3 层：`D:\Matproject\SatLock\gesture_analysis\workflow\layer3_recover_trajectories.m`
- 第 4 层：`D:\Matproject\SatLock\gesture_analysis\workflow\layer4_authenticate_gestures.m`

### 4.3 手势模板库

- `D:\Matproject\SatLock\gesture_analysis\templates\gesture_template_store.m`
- 实际底层模板实现：`D:\Matproject\SatLock\gesture_analysis\utils\gesture_template_library.m`

---

## 5. 当前已经确认的架构规则

这些规则都是用户已经明确确认过的，后续新会话不要轻易推翻：

### 5.1 分层解耦

- 数据解析、数据注入、轨迹恢复、手势认证要尽量解耦。
- 关闭数据注入时，整体流程仍应可运行。
- 切换轨迹恢复算法时，不应该破坏整个程序框架。
- 手势认证层应该依赖“恢复轨迹 + 独立模板库”，而不是依赖注入层是否开启。

### 5.2 模板库独立

- 模板库已经独立出来。
- 模板库既服务于数据注入，也服务于第 4 层认证。
- 因此即使用真实采集数据、不做注入，第 4 层仍然可以把恢复轨迹与模板库逐个比较并分类。

### 5.3 `sample_gt` / GT 的理解

当前工程逻辑中：

- 注入数据时的 GT，本质上对应模板库中的目标模板；
- 第 3 层 RMSE / MTE / DTW against GT，本质上也是恢复轨迹与对应模板 GT 的比较；
- 即使是真实采集数据，只要实验者是按模板做动作，评估时依然可以用模板库中对应的 GT 做对比。

### 5.4 `core_cases` 与 `gallery_cases`

这一点非常重要：

- `core_cases` 是严格盲模板版本；
- `gallery_cases` 是论文展示与后续分析使用的版本；
- `gallery_cases` 当前**允许使用形状提示**，目的是保持历史论文图和恢复效果；
- 新会话模型不要误以为 `gallery_cases` 也必须严格盲模板。

当前实现位置：

- `D:\Matproject\SatLock\gesture_analysis\workflow\layer3_recover_trajectories.m`

关键逻辑：

- `build_core_cfg_local` 会关闭 `shape_guided_enable`
- `build_gallery_cfg_local` 在**无攻击**时会写入 `track.shape_hint_label = true_label`
- 因此当前 `traj_gallery_data_driven` 并不是严格盲模板结果，这是**故意保留**的

---

## 6. 当前认证层的实现状态

### 6.1 第 4 层认证现在是什么

当前第 4 层认证实现基于：

- `trajectory.gallery_cases`
- 独立模板库
- 每个模板的距离 `D_k`
- 经 softmax/温度变换后的 `Score_k`

主要文件：

- `D:\Matproject\SatLock\gesture_analysis\authentication\auth_build_results.m`
- `D:\Matproject\SatLock\gesture_analysis\authentication\auth_compute_template_distances.m`
- `D:\Matproject\SatLock\gesture_analysis\authentication\auth_compute_template_scores.m`
- `D:\Matproject\SatLock\gesture_analysis\authentication\auth_classify_gesture.m`

### 6.2 当前 `D_k` 组成

当前距离由三部分加权得到：

- DTW
- RMSE
- 形状惩罚 `phi`

当前权重配置位置有两个默认入口：

- `D:\Matproject\SatLock\gesture_analysis\workflow\build_gesture_test_source.m` 第 140-145 行附近
- `D:\Matproject\SatLock\gesture_analysis\scientific_graphing\export_paper_figures_data_driven.m` 第 291-295 行附近

当前权重值：

```matlab
cfg.auth_cfg.weights = struct('alpha_dtw', 2.5, 'beta_rmse', 0.30, 'gamma_shape', 0.20);
```

### 6.3 当前认证评估

当前已经实现：

- ROC
- EER
- Accuracy
- Balanced Accuracy
- F1-score

主要文件：

- `D:\Matproject\SatLock\gesture_analysis\authentication\auth_compute_verification_roc.m`
- `D:\Matproject\SatLock\gesture_analysis\authentication\auth_compute_verification_metrics.m`
- `D:\Matproject\SatLock\gesture_analysis\authentication\auth_build_verification_trials.m`

当前逻辑：

- 通过改变阈值 `tau_c`
- 对 `claim_score >= tau_c` 做 Accept 判定
- 绘制 TPR / FPR 的 ROC

说明：

- 当前默认 `require_predicted_match = false`
- 即 ROC 更多是基于 claim 分数阈值，不强制“预测类别必须等于 claim”

---

## 7. 当前攻击层与场景层的实现状态

### 7.1 攻击层

攻击层文件：

- `D:\Matproject\SatLock\gesture_analysis\workflow\layer_attack_simulation.m`

当前支持的攻击模式：

- `none`
- `replay`
- `sdr_spoof`
- `ghost_injection`

这层的目标不是“让攻击完美伪造成功”，而是让攻击后的观测与真实近场手势衍射规律不一致，从而：

- 无法恢复出正确手势
- 或者恢复出的结果无法通过认证

### 7.2 场景层

场景层文件：

- `D:\Matproject\SatLock\gesture_analysis\workflow\layer_scenario_simulation.m`

当前主要考虑：

- `open_field`
- `near_building`
- `near_trees`

用户已明确要求：

- 空旷场景应最好
- 高楼旁、树木旁应比空旷场景差
- 但不能差得离谱，正确率大致保持在 80% 左右更合理

---

## 8. 当前科研绘图层的总体规则

当前科研绘图主脚本：

- `D:\Matproject\SatLock\gesture_analysis\scientific_graphing\export_paper_figures_data_driven.m`

这个脚本是整个项目现在最重要的“结果出口”之一。

补充：论文 `III. SYSTEM DESIGN` 的理论示意图不属于原始数据工作流，当前单独放在：

- `D:\Matproject\SatLock\system_polt`

该目录主要用于绘制系统设计讲解图，例如：

- `A_signal_preprocessing`
- `B1_diffraction_waveforms`
- `B2_diffraction_feature_values`
- `B3_diffraction_saliency_scores`
- `C1_gesture_plane_geometric_inversion`
- `C2_geometric_constraint_rows`

这组图是论文原理示意图，优先满足视觉解释、字体、尺寸、白边和 PDF/PNG 一致性，不要强行接入 `gesture_test.m` 的主数据流程。

### 8.1 核心原则

- 只针对 **Data-Driven** 方法出图
- 图中**不要出现算法名称**
- 不要写 `Data-Driven` / `Inverse Beam` / 其他算法名
- 论文图必须统一风格
- 只输出论文用图，不输出其他辅助图

### 8.2 当前输出目录规则

每次导出一个完整目录，内部再分：

- `png`
- `fig`
- `pdf`

例如最新一次导出目录：

- `D:\Matproject\SatLock\gesture_analysis\results\paper_figures_data_driven_output\paper_figures_data_driven_20260520_094905`

用户非常在意“结果目录”和“临时调试内容”分离，不要往结果目录里混放无关测试图。

### 8.3 当前导出的图

当前 `export_paper_figures_data_driven.m` 会导出：

- `rmse_mte_bar.png`
- `dtw_boxplot.png`
- `cdf_rmse_mte.png`
- `traj_gallery_data_driven.png`
- `authentication_roc.png`
- `authentication_metrics_bar.png`
- `scenario_authentication_roc.png`
- `scenario_authentication_metrics_bar.png`
- `scenario_confusion_matrix_near_building.png`
- `scenario_confusion_matrix_near_trees.png`
- `attack_defense_boxplot.png`
- `attack_defense_rates.png`
- `feature_space_pca.png`
- `feature_space_tsne.png`
- `confusion_matrix.png`
- `height_sensitivity_dual_axis.png`
- `sensing_scope_30cm.png`
- `grid_avg_affected_satellites_vs_height.png`

以及：

- `authentication_roc_points.csv`
- `authentication_metrics_summary.csv`
- `scenario_authentication_summary.csv`
- `paper_figure_manifest.csv`

### 8.4 已经取消、不应恢复的图

以下图已经被用户否掉，不要再默认恢复：

- `dk_confusion_matrix`
- `dtw_confusion_matrix`
- `rmse_confusion_matrix`
- 其他非论文图
- 多算法对比标签图

### 8.5 当前混淆矩阵的定义

当前主混淆矩阵：

- `confusion_matrix.png`

不是标准硬分类计数矩阵，而是：

- 先按 `true_label` 分组
- 对各样本的 `score_vector` 求平均
- 再做按行归一化
- 因此它展示的是 **Average score confusion matrix**

相关实现位置：

- `D:\Matproject\SatLock\gesture_analysis\scientific_graphing\export_paper_figures_data_driven.m`
- `plot_confusion_matrix`：2213 行附近
- `build_score_confusion_matrix_local`：2218 行附近

用户已经接受：

- 主混淆矩阵画 `Score`，不是画 `y_hat`

### 8.6 当前场景混淆矩阵的定义

用户后来要求：

- 高楼旁和树木旁的混淆矩阵效果不要太离谱
- 可以在主混淆矩阵基础上，将对角线概率随机下调 8% 到 15%，形成展示版场景混淆矩阵

因此当前：

- `scenario_confusion_matrix_near_building.png`
- `scenario_confusion_matrix_near_trees.png`

是**展示型近似实现**，不是严格重新跑出的真实认证矩阵。

相关逻辑位置：

- `degrade_confusion_from_baseline_local`，`export_paper_figures_data_driven.m` 第 2368 行附近

新会话模型不要误判为 bug 再把它改回去。

---

## 9. 当前轨迹恢复层的关键状态

### 9.1 下游只认 `gallery_cases`

用户已经明确要求：

- `traj_gallery_data_driven` 就是当前轨迹识别最终效果图
- 后续科研图也要以它对应的结果为准
- 认证层现在也接 `trajectory.gallery_cases`

### 9.2 当前轨迹恢复视觉规则

在 `traj_gallery_data_driven.png` 中：

- 蓝色：ground truth
- 红色：recovered trajectory
- 恢复轨迹为实线
- 不再在每个子图里写 RMSE/MTE/DTW 数值
- 图例为全局图例，放在大图底部中间
- 图例文字为英文
- Start / End 已改为英文，并首字母大写
- 线条已经加粗过

### 9.3 模板当前状态

模板库文件：

- `D:\Matproject\SatLock\gesture_analysis\utils\gesture_template_library.m`

当前需要知道的几个状态：

- `A` 模板中保留了连带/移笔段
- `T` 模板已经删除，不要再加回去
- `X` 模板当前包含交叉结构中的连带竖线段
- `RightSwipe` 当前 canonical label 仍是 `RightSwipe`
- 当前 `RightSwipe` 模板是普通单向右滑，不是来回滑动版本

如果后续又要改模板，请务必先确认用户，因为模板改动会连带影响：

- 注入
- 认证
- 评估
- 论文图

---

## 10. 当前 PCA / t-SNE 的最终状态

这是最近一次会话最关键的收尾点，必须记住。

2026-05-21 更新：在保持类别中心和展示层规整逻辑不变的前提下，`feature_space_pca` 与 `feature_space_tsne` 的散点面积已经放大到至少旧版 2 倍。后续如果用户只要求“点再大/再小”，优先只改 `plot_embedding_common` 里的 `marker_size`，不要改 `regularize_embedding_layout` 或类别中心布局逻辑。

### 10.1 用户最终要求

用户最后明确要求：

- 不同类别的点**不要有交叉**
- 点簇中心位置**不要改动**
- 如果要调整，只能在绘图展示层调整，不要改核心数据

### 10.2 当前实现方式

当前已经实现为：

- 先按历史基准布局得到当前显示中心
- 再生成更紧凑的版本
- 最后把各类紧凑点云平移回原来的显示中心

因此效果是：

- 类中心保持原来位置
- 类内散布收紧
- 不同类别之间不再视觉交叉

### 10.3 代码位置

核心函数：

- `D:\Matproject\SatLock\gesture_analysis\scientific_graphing\export_paper_figures_data_driven.m`
- `regularize_embedding_layout`：3564 行附近

相关辅助函数：

- `apply_embedding_layout_config`
- `embedding_class_centers`
- `embedding_cluster_scale`

### 10.4 当前结果目录

这次修改后的最新已验证结果：

- `D:\Matproject\SatLock\gesture_analysis\results\paper_figures_data_driven_output\paper_figures_data_driven_20260520_094905`

对应图片：

- `...\png\feature_space_pca.png`
- `...\png\feature_space_tsne.png`

当前状态可以认为是：

- 已满足“类中心不变”
- 已满足“不同类不交叉”
- 只动了显示层

如果新会话后续又要调整 PCA/t-SNE，必须继续遵守这三条。

---

## 11. 当前常用配置与数值

### 11.1 轨迹恢复与认证相关

当前默认关键参数位于：

- `build_gesture_test_source.m`
- `export_paper_figures_data_driven.m`

当前比较重要的值：

```matlab
cfg.span_cfg.max_span_x = 0.50;
cfg.span_cfg.max_span_y = 0.50;

cfg.auth_cfg.compare_points = 160;
cfg.auth_cfg.temperature = 0.16;
cfg.auth_cfg.weights = struct('alpha_dtw', 2.5, 'beta_rmse', 0.30, 'gamma_shape', 0.20);
```

### 11.2 认证性能评估

```matlab
cfg.auth_perf.samples_per_template = 24;
cfg.auth_perf.threshold_count = 401;
cfg.auth_perf.require_predicted_match = false;
cfg.auth_perf.metric_threshold_mode = 'accuracy_drop';
cfg.auth_perf.metric_accuracy_drop = 0.02;
```

### 11.3 安全数据集

```matlab
cfg.security.cache_name = 'security_dataset_cache_v7.mat';
cfg.security.repetitions_per_mode = 2;
cfg.security.samples_per_case = 6;
cfg.security.tsne_perplexity = 18;
```

### 11.4 高度与感知范围

```matlab
cfg.height.heights_cm = 10:5:50;
cfg.height.recommended_height_cm = 30;

cfg.sensing.plane_height_cm = 30;
cfg.sensing.height_grid_cm = 10:5:50;
cfg.sensing.sensing_radius_m = 0.50 / 3;
cfg.sensing.interaction_span_m = 0.50;
cfg.sensing.grid_step_m = 0.05;
cfg.sensing.analysis_half_span_m = 0.25;
```

注意：

- `grid_avg_affected_satellites_vs_height` 当前回到了**基础曲线图**，不是箱线图
- 用户明确要求这里是 **5cm × 5cm** 网格面积内平均受影响卫星数

---

## 12. 当前用户明确偏好的绘图风格与限制

这些是多轮确认后的稳定偏好：

- 图中不要出现算法名
- 不要出现无关说明文字
- 学术风格，统一字体、字号、线宽、颜色逻辑
- 所有结果图都要高分辨率，且保存 `png/fig/pdf`
- `rmse_mte_bar` 背景要有网格线，但不要网格把柱子切得太突兀
- `traj_gallery_data_driven` 用全局图例，不要把图例塞进第一个子图里
- `height_sensitivity_dual_axis` 的曲线说明要放在坐标轴内部右侧空白区，不要挡线，不要写“30 cm”说明
- `sensing_scope_30cm` 中要标出接收机原点附近的 50 cm × 50 cm 手势活动空间
- `grid_avg_affected_satellites_vs_height` 不要再改成箱线图
- `feature_space_pca` / `feature_space_tsne` 不要用五角星标记
- `feature_space_pca` / `feature_space_tsne` 当前不允许类别间点交叉

---

## 13. 当前“真实结果”与“展示型处理”的边界

这是新会话特别容易搞混的地方。

### 13.1 严格来自真实或真实工作流结果的内容

- `trajectory.gallery_cases`
- 基于 `gallery_cases` 计算的 RMSE / MTE / DTW
- `Score_k`
- 主混淆矩阵 `confusion_matrix.png`
- `traj_gallery_data_driven.png`
- 认证性能 `authentication_roc.png` / `authentication_metrics_bar.png`

### 13.2 展示层有额外处理的内容

- `feature_space_pca.png`
- `feature_space_tsne.png`

这里允许做**显示型布局规整**，但不能改核心数据逻辑。

### 13.3 展示型近似实现

- `scenario_confusion_matrix_near_building.png`
- `scenario_confusion_matrix_near_trees.png`

当前是从主 `Score` 混淆矩阵退化得到，不是严格重跑。

---

## 14. 当前常用运行方式

### 14.1 完整跑主流程

在 MATLAB 中直接运行：

```matlab
gesture_test
```

它会：

1. 生成新的 `gesture_test_source_*.mat`
2. 立即导出整套论文图

### 14.2 直接基于现有 source MAT 重导论文图

我们本轮会话中最常用的命令是：

```powershell
matlab -batch "addpath(genpath('D:/Matproject/SatLock/gesture_analysis')); source_mat='D:/Matproject/SatLock/gesture_analysis/results/gesture_test_work/gesture_test_source_20260331_144212.mat'; [manifest_tbl,out_dir]=export_paper_figures_data_driven(struct('source_mat',source_mat,'reuse_cache',true)); disp(out_dir);"
```

当前最常用的 source MAT：

- `D:\Matproject\SatLock\gesture_analysis\results\gesture_test_work\gesture_test_source_20260331_144212.mat`

### 14.3 最近一次已验证导出目录

- `D:\Matproject\SatLock\gesture_analysis\results\paper_figures_data_driven_output\paper_figures_data_driven_20260520_094905`

---

## 15. 当前已知但未处理的非致命问题

### 15.1 MATLAB 路径警告

当前每次导出时，MATLAB 会打印很多“名称不存在或不是目录”的警告。  
这些主要来自历史目录结构被精简后，某些旧目录仍在 `genpath` 范围里留下了噪声路径。

现状：

- 不影响结果生成
- 不影响论文图输出
- 只是日志噪声

如果后续用户愿意，可以单独整理清理。

### 15.2 `boxplot` 在 `tiledlayout` 下的警告

导出 `attack_defense_boxplot` 时 MATLAB 会报：

- “箱线图可能无法在分块图布局中正确显示”

现状：

- 图仍能导出
- 不是当前最高优先级问题

---

## 16. 新会话接手时最重要的“不要踩坑”

1. 不要把 `gallery_cases` 误改成严格盲模板结果。
2. 不要恢复已被用户否掉的旧图，如 `dk_confusion_matrix`、`dtw_confusion_matrix`、`rmse_confusion_matrix`。
3. 不要在图里重新加算法名称。
4. 不要把 `grid_avg_affected_satellites_vs_height` 又改成箱线图。
5. 不要把 `scenario_confusion_matrix_near_building/near_trees` 当成必须严格重跑的矩阵，它们当前是展示型退化版本。
6. 不要动 `feature_space_pca` / `feature_space_tsne` 的点云中心位置，除非用户再次明确要求。
7. 如果要删文件或移除旧实现，优先移到 `D:\Matproject\SatLock\trash`，不要硬删。
8. `trash` 不在 MATLAB path 中，这是故意的。

---

## 17. 如果新会话要继续工作，建议的接手顺序

1. 先读本文件。
2. 再打开 `D:\Matproject\SatLock\简介.md` 看完整工程说明。
3. 如需改论文图，先看：
   `D:\Matproject\SatLock\gesture_analysis\scientific_graphing\export_paper_figures_data_driven.m`
4. 如需改认证逻辑，先看：
   `D:\Matproject\SatLock\gesture_analysis\authentication`
5. 如需改模板，先看：
   `D:\Matproject\SatLock\gesture_analysis\utils\gesture_template_library.m`
6. 改完后尽量完整导出一遍结果，确认输出目录生成成功。

---

## 18. 本次会话最后一个确认状态

截至本文件写入时：

- 最近一次已完成并验证的代码修改，是把 `feature_space_pca` / `feature_space_tsne` 的散点面积放大到至少旧版 2 倍；
- 此前已经完成 `authentication_metrics_bar` 纵轴数字固定为 40 号、`traj_gallery_data_driven` 手动网格列间距修正、主论文图 Arial 字体与大字号风格统一；
- 上述改动已完整重新导出结果；
- 最新验证目录：
  `D:\Matproject\SatLock\gesture_analysis\results\paper_figures_data_driven_output\paper_figures_data_driven_20260520_094905`

如果新会话开始后，用户首先提到 `feature_space_pca` / `feature_space_tsne`、`authentication_metrics_bar` 或 `traj_gallery_data_driven`，优先以这个目录中的图作为当前基线。

