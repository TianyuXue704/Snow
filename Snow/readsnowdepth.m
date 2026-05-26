
clear; close all
%% 1. 路径设置
input_dir = 'C:\Users\Jinghong Xu\Desktop\OIB_raw\data\oib_raw\OIB_data\2019';
output_dir = fullfile(input_dir, 'Processed_MAT');

if ~exist(output_dir, 'dir'), mkdir(output_dir); end

%% 2. 获取文件列表
file_list = dir(fullfile(input_dir, '*.txt'));
num_files = length(file_list);
summary_info = cell(num_files, 5); 

fprintf('开始处理 %d 个文件，经度将自动转换为 [-180, 180]...\n', num_files);

%% 3. 循环处理
for i = 1:num_files
    curr_file_name = file_list(i).name;
    curr_full_path = fullfile(input_dir, curr_file_name);
    
    try
        data = readtable(curr_full_path, 'VariableNamingRule', 'preserve');
        
        % 过滤无效值
        idx = (data.lat ~= -99999) & (data.lon ~= -99999) & (data.snow_depth ~= -99999);
        lat = data.lat(idx);
        lon = data.lon(idx);
        snow_depth = data.snow_depth(idx);
        thickness = data.thickness(idx);
        
        if isempty(lat), continue; end
        
        % --- 关键步骤：经度转换 ---
        % 如果经度大于 180，则减去 360 转换到负数区间（西经）
        lon(lon > 180) = lon(lon > 180) - 360;
        
        % 计算新的范围
        lat_range = [min(lat), max(lat)];
        lon_range = [min(lon), max(lon)];
        
        % 汇总信息
        summary_info(i, :) = {curr_file_name, lat_range(1), lat_range(2), lon_range(1), lon_range(2)};
        
        % 保存文件
        [~, name_only, ~] = fileparts(curr_file_name);
        save(fullfile(output_dir, [name_only, '.mat']), 'lat', 'lon', 'snow_depth', 'thickness');
        
        fprintf('已转换: %s | 经度范围: [%.2f, %.2f]\n', curr_file_name, lon_range(1), lon_range(2));
        
    catch ME
        fprintf('跳过文件 %s: %s\n', curr_file_name, ME.message);
    end
end

%% 4. 导出汇总报告
summary_table = cell2table(summary_info, 'VariableNames', ...
    {'FileName', 'Min_Lat', 'Max_Lat', 'Min_Lon', 'Max_Lon'});
writetable(summary_table, fullfile(input_dir, 'Data_Summary_Standard_Lon.csv'));

disp('--- 转换完成！---');