clc; clear; close all;

% 定义路径常量
path = 'C:\Users\TianyuXue\Desktop\SnowCode\Data';
SavePath = 'C:\Users\TianyuXue\Desktop\SnowCode\Data';

% 创建保存目录（如需）
if ~exist(SavePath, 'dir')
    mkdir(SavePath);
    fprintf('SavePath did not exist. Created directory: %s\n', SavePath);
end

% 获取.h5文件列表
fileList = dir(fullfile(path, '*.h5'));
num_h5 = length(fileList);


for i = 1  % 防止索引超出范围
    currentFile = fileList(i);
    h5Path = fullfile(currentFile.folder, currentFile.name);
    
    % 显示当前进度
    fprintf('\nProcessing file %d/%d: %s\n', i, num_h5, currentFile.name);
    
    try
        % 生成目标.mat文件名
        baseName = erase(currentFile.name, '.h5');
        matFile1 = fullfile(SavePath, ['Ngt_h5data_strg3_', baseName, '.mat']);
        matFile2 = fullfile(SavePath, ['Day_h5data_strg3_', baseName, '.mat']);
        
        % 跳过已存在的文件
        if exist(matFile1, 'file') || exist(matFile2, 'file')
            fprintf('Skipping existing files for %s\n', baseName);
            continue;
        end
        
        % 核心处理函数（添加try-catch内部二次保护）
        f_read_ATL03_strong(h5Path, SavePath);  % 假设此函数无返回值
        
        % 清理资源
        close all;
        fprintf('Successfully processed: %s\n', baseName);
        
    catch ME  % 捕获所有异常
        % 记录错误信息
        errorTime = datestr(now, 'yyyy-mm-dd HH:MM:SS');
        errorMsg = sprintf('[%s] ERROR in file %s:\n%s\n',...
            errorTime, currentFile.name, ME.message);
        
        % 打印到命令行和日志文件
        fprintf(2, '%s', errorMsg);  % 2表示标准错误
        fprintf(fid_log, '%s', errorMsg);
        
        % 清理失败文件的残留（可选）
        if exist(matFile1, 'file'), delete(matFile1); end
        if exist(matFile2, 'file'), delete(matFile2); end
        
        % 确保图形窗口关闭
        close all force;  
    end
    
    % 防止内存泄漏
    clearvars -except path SavePath fileList num_h5 i fid_log
end

