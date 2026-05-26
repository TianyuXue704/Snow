clc; clear; close all; 
z = 1; 
input_root = 'C:\Users\Jinghong Xu\Desktop\OIB\OIB\2019\0419\ATL03_mat_data\';
result_root = 'C:\Users\Jinghong Xu\Desktop\OIB\OIB\2019\0419\result2\';

if ~exist(result_root, 'dir') 
    mkdir(result_root); 
end  
a = dir(input_root); 
num = length(a); 
Dayflag = 'Day'; 
for i = 5                                                                                                                                                                                                                                    
    % 提取文件名中的月份（第29和30个字符）


        fprintf('Processing file %d: %s\n', i, a(i).name); 
        file_choose = fullfile(a(i).folder, a(i).name); 
        file_profile = f_ocean_choose(file_choose, Dayflag, result_root); 


   close all; 
 end 



clc; clear; close all;


% --- 路径设置 ---
path = 'C:\Users\Jinghong Xu\Desktop\OIB\OIB\2019\0419\result2\*.mat';
save_dir = 'C:\Users\Jinghong Xu\Desktop\OIB\OIB\2019\0419\result_snow\'; % 结果保存路径

% 如果保存文件夹不存在则创建
if ~exist(save_dir, 'dir')
    mkdir(save_dir);
end

files = dir(path);
num = length(files);

fprintf('开始处理并保存数据，共 %d 个文件...\n', num);

for i = 1:num
    fprintf('正在处理 (%d/%d): %s\n', i, num, files(i).name);
    
    file_path = fullfile(files(i).folder, files(i).name);
    
    % --- 1. 数据处理阶段 ---
    % 执行你的原始函数
    file_deconv = f_deconvolution(file_path);
    [snow, data_S, snowA] = snow_reterival(file_deconv);
    
    % --- 2. 直接保存结果 ---
    % 构造保存文件名（Result_0419_原文件名.mat）
    save_file_name = fullfile(save_dir, ['Result', files(i).name]);
    
    % 只保存反演计算出的关键变量
    save(save_file_name, 'snow', 'data_S', 'snowA'); 
end

fprintf('所有结果已成功保存至: %s\n', save_dir);