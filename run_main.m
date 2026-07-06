

clc; clear; close all;


% ATL03 原始 .mat 文件所在文件夹（支持通配符批量处理）
data_path = 'C:\Users\TianyuXue\Desktop\Np_correction';

% 中间结果保存路径（profile 和 deconv 文件）
save_path = 'C:\Users\TianyuXue\Desktop\Np_correction\Output\';

% 最终雪深结果保存路径
result_path = 'C:\Users\TianyuXue\Desktop\Np_correction\Result\';

%（可选）对应本次处理的 ATL09 文件路径：
% 填空则保持原来的 Ta=0.9 粗略估计；填入后将用 ATL09 提取/匹配逐廓线 Ta
atl09_file = 'C:\Users\TianyuXue\Desktop\Np_correction\ATL09_read\ATL09_20230102055504_01941801_007_01.mat';
atl09_Ta_dataset = '';
atl09_gt = '';
atl09_radius_m = 2000;  % 2km 
Ta_is_twoway = false;
atl09_mat_lat_var = 'initial_lat';
atl09_mat_lon_var = 'initial_lon';
atl09_mat_Ta_var  = 'column_od_asr';

% 数据文件名匹配规则
file_pattern = 'DataDay_h5data_strg3_ATL03_*.mat';

%% ============================================================

if ~exist(save_path,   'dir'); mkdir(save_path);   end
if ~exist(result_path, 'dir'); mkdir(result_path); end

files = dir(fullfile(data_path, file_pattern));
if isempty(files)
    error('未找到匹配文件，请检查 data_path 和 file_pattern');
end
fprintf('共找到 %d 个 ATL03 文件\n\n', length(files));

% 汇总结果
all_lat  = [];
all_snow = [];
all_snowA = [];

for i = 1:length(files)
    file_choose = fullfile(files(i).folder, files(i).name);
    fprintf('========== [%d/%d] %s ==========\n', i, length(files), files(i).name);

    
    %% 第一步：廓线生成
    
    fprintf('  [步骤1] 廓线生成 (f_profile_processing)...\n');
    try
        Dayflag = 'Ngt';
        % 只计算 lat [-65, -50] 范围内的廓线（与ATL09重叠区域）
        lat_range = [-65, -50];
        f_profile_processing(file_choose, Dayflag, save_path, lat_range);
        % 推算输出文件名
        [~, stem, ~] = fileparts(file_choose);
        profile_file = fullfile(save_path, [stem, '_profile.mat']);
    catch ME
        fprintf('  [跳过] 廓线生成失败: %s\n', ME.message);
        close all; continue;
    end

    if ~exist(profile_file, 'file')
        fprintf('  [跳过] profile 文件未生成\n');
        close all; continue;
    end

    %% ----------------------------------------------------------
    %% 第二步：退卷积（系统响应去除）
   
    %% ----------------------------------------------------------
    fprintf('  [步骤2] 退卷积 (f_deconvolution)...\n');
    try
        file_deconv = f_deconvolution(profile_file);
    catch ME
        fprintf('  [跳过] 退卷积失败: %s\n', ME.message);
        close all; continue;
    end

    if isempty(file_deconv) || ~exist(file_deconv, 'file')
        fprintf('  [跳过] deconv 文件未生成\n');
        close all; continue;
    end

    %% ----------------------------------------------------------
    %% 第三步：PMT修正 + 雪深反演
    %% ----------------------------------------------------------
    fprintf('  [步骤3] PMT修正 + 雪深反演 (snow_reterival_corrected)...\n');
    try
        [snow, lat, snowA, step_info] = snow_reterival_corrected( ...
            file_deconv, ...           
           'Ta',      0.9, ...   % 单程大气透过率（学长给出）
           'atl09_file', atl09_file, ...
           'atl09_Ta_dataset', atl09_Ta_dataset, ...
           'atl09_gt', atl09_gt, ...
           'atl09_radius_m', atl09_radius_m, ...
           'Ta_is_twoway', Ta_is_twoway, ...
           'atl09_mat_lat_var', atl09_mat_lat_var, ...
           'atl09_mat_lon_var', atl09_mat_lon_var, ...
           'atl09_mat_Ta_var',  atl09_mat_Ta_var, ...
           'CF',      0.52, ...   % 系统在轨标定因子（论文Antarctic标定）
           'beam',    'strong', ...
           'verbose', true);
    catch ME
        fprintf('  [跳过] 雪深反演失败: %s\n', ME.message);
        close all; continue;
    end

    %% ----------------------------------------------------------
    %% 第四步：保存单文件结果
    %% ----------------------------------------------------------
    valid = ~isnan(snow);
    if sum(valid) == 0
        fprintf('  [提示] 无收敛结果，跳过保存\n');
        close all; continue;
    end

    [~, dstem, ~] = fileparts(file_deconv);
    out_file = fullfile(result_path, [dstem, '_snow.mat']);
    save(out_file, 'snow', 'lat', 'snowA', 'step_info');
    fprintf('  [保存] %s\n', out_file);

    % --- 每条轨道单独保存雪深图 ---
    fig_track = figure('Color','w','Visible','off','Position',[50 200 1400 380]);
    [lat_s, si_t] = sort(lat(valid));
    snow_s = snow(valid);
    snow_s = snow_s(si_t);
    plot(lat_s, snow_s, 'Color',[0 0.4470 0.7410], 'LineWidth',0.8);
    hold on;
    yline(mean(snow_s),'--','Color',[0.85 0.33 0.10],'LineWidth',1.2, ...
          'Label',sprintf('均值 %.3f m', mean(snow_s)), ...
          'LabelHorizontalAlignment','left','FontSize',9);
    % 叠加反照率（右轴）
    yyaxis right
    snowA_s = snowA(valid); snowA_s = snowA_s(si_t);
    plot(lat_s, snowA_s, 'Color',[0.47 0.67 0.19], 'LineWidth',0.6, 'LineStyle','--');
    ylabel('反照率', 'FontSize',10, 'Color',[0.47 0.67 0.19]);
    set(gca,'YColor',[0.47 0.67 0.19]);
    yyaxis left
    set(gca,'XDir','reverse','FontSize',10,'XMinorTick','on','YMinorTick','on');
    grid on; box on;
    y_mg = 0.02;
    ylim([max(0, min(snow_s)-y_mg), max(snow_s)+y_mg]);
    xlim([min(lat_s)-0.05, max(lat_s)+0.05]);
    xlabel('纬度 (°)','FontSize',11);
    ylabel('雪深 (m)','FontSize',11);
    title(strrep(dstem,'_',' '),'FontSize',11,'FontWeight','normal','Interpreter','none');
    fig_jpg = fullfile(result_path, [dstem, '_snow.jpg']);
    saveas(fig_track, fig_jpg);
    close(fig_track);
    fprintf('  [图像] %s\n', fig_jpg);

    all_lat   = [all_lat,   lat(valid)];
    all_snow  = [all_snow,  snow(valid)];
    all_snowA = [all_snowA, snowA(valid)];

    close all;
end

%% ============================================================
%% 汇总绘图
%% ============================================================
if isempty(all_snow)
    fprintf('\n[结束] 无有效结果\n');
    return;
end

fprintf('\n========== 全部处理完成 ==========\n');
fprintf('有效廓线总数 : %d\n',      length(all_snow));
fprintf('雪深范围     : [%.3f, %.3f] m\n', min(all_snow), max(all_snow));
fprintf('雪深均值     : %.3f m\n',  mean(all_snow));
fprintf('反照率均值   : %.3f\n',    mean(all_snowA));

% --- 按纬度排序 ---
[all_lat_sorted, si] = sort(all_lat);
all_snow_sorted      = all_snow(si);

% --- 宽幅单图：雪深随纬度变化 ---
figure('Color','w','Position',[50 200 1600 420]);

plot(all_lat_sorted, all_snow_sorted, ...
     'Color',[0 0.4470 0.7410], 'LineWidth', 0.8);
hold on;

% 均值参考线
ybar = mean(all_snow_sorted);
yline(ybar, '--', sprintf('均值 %.3f m', ybar), ...
      'Color',[0.8500 0.3250 0.0980], 'LineWidth',1.2, ...
      'LabelHorizontalAlignment','left', 'FontSize',10);

set(gca, 'XDir','reverse');
set(gca, 'FontSize', 11);
set(gca, 'XMinorTick','on', 'YMinorTick','on');
grid on; box on;

% 纵轴放大细节
y_margin = 0.02;
ylim([max(0, min(all_snow_sorted) - y_margin), ...
      max(all_snow_sorted) + y_margin]);

% 横轴紧贴数据范围
lat_margin = 0.05;
xlim([min(all_lat_sorted) - lat_margin, ...
      max(all_lat_sorted) + lat_margin]);

xlabel('纬度 (°)', 'FontSize', 12);
ylabel('雪深 (m)',  'FontSize', 12);
title('ICESat-2 雪深沿轨分布', 'FontSize', 13, 'FontWeight','normal');

saveas(gcf, fullfile(result_path, 'snow_depth_vs_lat.jpg'));
fprintf('结果图已保存至 %s\n', fullfile(result_path,'snow_depth_vs_lat.jpg'));
