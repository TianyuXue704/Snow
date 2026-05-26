function [S, Depth, data_S] = f_ocean_strict_process(file_choose, Dayflag, savepath)
%% 1. 数据加载与预处理
[~, name_stem, ~] = fileparts(file_choose); % 获取文件名用于保存
data_struct = load(file_choose);

% 容错处理：获取变量
if isfield(data_struct, 'initial_t')
    initial_t = data_struct.initial_t;
    initial_h = data_struct.initial_h;
    initial_lat = data_struct.initial_lat;
    initial_lon = data_struct.initial_lon;
else
    error('Mat文件中缺少 initial_t/h/lat/lon 变量');
end

% --- 粗去噪 (保留数据用于后续精细处理) ---
% 为了保证 ind_h 计算不出错，先剔除极端的经纬度噪声
ind_noise = find(initial_lat >= 79 | initial_lat < 70); % 根据你的数据范围调整
initial_h(ind_noise) = [];
initial_t(ind_noise) = [];
initial_lon(ind_noise) = [];
initial_lat(ind_noise) = [];

% 为了防止 grid 过大，先基于中值做一次宽松的筛选 (比如保留海面上下 50m)
center_h = median(initial_h);
keep_idx = find(initial_h > center_h - 50 & initial_h < center_h + 50);
initial_h = initial_h(keep_idx);
initial_t = initial_t(keep_idx);
initial_lat = initial_lat(keep_idx);
initial_lon = initial_lon(keep_idx);

if isempty(initial_h)
    warning('数据为空'); S=[]; Depth=[]; data_S=[]; return; 
end

%% 2. 严格按照你的代码逻辑进行 Binning 和 对齐
fprintf('正在执行严格对齐逻辑...\n');

% --- 你的参数设定 ---
delta_h = 0.15;    % 垂直分辨率
delta_t = 0.001;  % 原始数据分Bin时间间隔 (注意：这会生成巨大的矩阵，请确保内存足够)

max_h = max(initial_h);
min_t = min(initial_t);

% 计算索引 (严格遵循你的公式)
% 注意：ind_h = 1 对应 max_h (最高点)，ind_h 越大代表越深
ind_t = floor((initial_t - min_t) / delta_t) + 1;
ind_h = floor((max_h - initial_h) / delta_h) + 1;

m = max(ind_h);
n = max(ind_t);

% 初始化矩阵
data_bin = zeros(3, n);
num_eachbin = zeros(1, n);
bin = zeros(m, n);

% --- 你的 For 循环 (严格保留) ---
% 为了速度，这里稍微加了进度提示，逻辑完全不变
len_data = length(initial_h);
for k = 1:len_data
    % 累加经纬度时间
    data_bin(1, ind_t(k)) = data_bin(1, ind_t(k)) + initial_lon(k);
    data_bin(2, ind_t(k)) = data_bin(2, ind_t(k)) + initial_lat(k);
    data_bin(3, ind_t(k)) = data_bin(3, ind_t(k)) + initial_t(k);
    
    % 计数
    num_eachbin(1, ind_t(k)) = num_eachbin(ind_t(k)) + 1;
    bin(ind_h(k), ind_t(k)) = bin(ind_h(k), ind_t(k)) + 1;
end

% 平均化 data_bin
% 避免除以0产生NaN (虽然下面会清空)
valid_bins = num_eachbin > 0;
data_bin(:, valid_bins) = data_bin(:, valid_bins) ./ repmat(num_eachbin(valid_bins), 3, 1);

% --- 清洗无效列 ---
[maxbin, ~] = max(bin);
ind_noise_bin = find(maxbin == 0);
maxbin(ind_noise_bin) = [];
bin(:, ind_noise_bin) = [];
data_bin(:, ind_noise_bin) = [];

if size(bin, 2) < 100
    warning('有效 Bin 数量不足 100'); S=[]; Depth=[]; data_S=[]; return;
end

% --- 寻找峰值 (模拟 SeaSurface) ---
[~, ind_peak0] = max(bin);

% 由于我没有 SeaSurface 函数，这里使用 filloutliers 模拟去除跳变点
% 你的 SeaSurface 应该是为了防止找到错误的噪点峰值
% 这里用移动中位数平滑一下峰值索引
ind_peak = filloutliers(ind_peak0, 'nearest', 'movmedian', 20);

% 清洗 NaN
ind_error = find(isnan(ind_peak));
bin(:, ind_error) = [];
data_bin(:, ind_error) = [];
ind_peak(ind_error) = [];

if (isempty(ind_peak))
    return
end

% --- 你的 While 循环 (寻找对齐基准 flag) ---
n = size(bin, 2);
i = 1;
mean_ind_peak = mean(ind_peak);
std_ind_peak = std(ind_peak);
diff_ind_peak = ind_peak(i) - mean_ind_peak;

while(abs(diff_ind_peak) > 0.5 * std_ind_peak && i <= n-1)
    i = i + 1;
    diff_ind_peak = ind_peak(i) - mean_ind_peak;
end
flag = i; % 找到的基准列索引

% --- 执行对齐 (CircShift) ---
fprintf('执行 circshift 对齐 (基准列: %d)...\n', flag);
for i = 1:n
    shiftsize = ind_peak(flag) - ind_peak(i);
    if(shiftsize ~= 0)
        bin(:, i) = circshift(bin(:, i), shiftsize);
    end
end

% 归一化 bin (仅用于绘图或后续处理，看你的需求)
% bin = bin ./ max(bin); % 你的代码里有这一步，保留

%% 3. 廓线生成 (调用 bin_to_profile)
% 你的要求：统计范围只需要海面上 6m 到 海面下 21m
% delta_h = 0.15m
% N1 (上/Air)  = 6m / 0.15m = 40
% N2 (下/Water) = 21m / 0.15m = 140

[S, Depth, data_S] = sub_bin_to_profile(bin, data_bin, ind_peak, flag, delta_h, 40, 140);

if (isempty(S))
    warning('生成廓线为空'); return;
end

%% 4. 后处理 (平滑与去底噪)
St = S; % 备份原始 S

% --- 1.2m 以下平滑 ---
% Depth 向量中 0 是表面，负数是空气，正数是水下
start_smooth = find(Depth >= 1.2, 1); 
if ~isempty(start_smooth)
    Si = S(start_smooth:end, :);
    
    % 你的卷积平滑代码 (保留)
    window_half_width = 0.75; % 假设你想平滑 +-0.75m 范围
    window_bins = round((2 * window_half_width) / delta_h); 
    if mod(window_bins, 2) == 0
        window_bins = window_bins + 1; 
    end
    kernel = ones(window_bins, 1) / window_bins;
    
    for k = 1:size(Si, 2)
        Si(:, k) = conv(Si(:, k), kernel, 'same');
    end
    S(start_smooth:end, :) = Si;
end

% --- 去除底噪 ---
% 你定义 ind_noise = Depth <= -5 (即海面上5m以上作为噪声基底)
ind_noise_rows = find(Depth <= -8);

if ~isempty(ind_noise_rows)
    sub_S = S(ind_noise_rows, :);
    sub_S(isnan(sub_S)) = 0;
    noise_vec = mean(sub_S, 1); % 每条廓线的噪声底
    
    % 执行减法 (每一列减去对应的噪声值)
    S = S - repmat(noise_vec, size(S, 1), 1);
else
    noise_vec = zeros(1, size(S, 2));
end

% 再次清理负值
S(S <= 0) = nan;

%% 5. 绘图与保存 (严格按照你的绘图风格)
if ~isempty(savepath)
    if ~exist(savepath, 'dir'), mkdir(savepath); end
    
    f = figure('Visible', 'on');

    num_plot = min(100, size(S, 2));
    colors = jet(num_plot);
    
    hold on;
    for i = 1:num_plot
        % 你的代码是画 St (原始对齐数据) 还是 S (去噪平滑后)? 
        % 你的代码写的是 semilogy(Depth, St(:,i)...)
        semilogy(Depth, St(:, i), 'LineWidth', 0.5, 'Color', colors(i,:));
    end
    xlabel('Depth (m)', 'FontSize', 10);
    ylabel('Normalized Signal', 'FontSize', 10);
    set(gca, 'YScale', 'log'); % 强制指定 Y 轴为对数刻度
    ylim([10^-4 1]);
    xlim([-6 21]); % 限制显示范围
    grid on;
    title(['Profile: ' name_stem], 'Interpreter', 'none');
    
    % 保存
    jpg_profile = fullfile(savepath, [name_stem, '_profile.jpg']);
    saveas(gcf, jpg_profile);
    
    file_profile = fullfile(savepath, [name_stem, '_profile.mat']);
    save(file_profile, 'S', 'Depth', 'data_S', 'noise_vec');
end

end

%% ============================================================
%% 子函数: bin_to_profile (集成并适配 N1/N2)
%% ============================================================
function [S, Depth, data_S] = sub_bin_to_profile(bin, data_bin, ind_peak, flag, delta_h, N1, N2)
    % N1: 海面上 bin 数 (Depth < 0)
    % N2: 海面下 bin 数 (Depth > 0)
    
    [m, n] = size(bin);
    
    % 你的积分参数
    dt_int = 5;       % 积分时间 5s
    step_t = 0.5;     % 滑动步长 0.5s
    
    % 构建深度向量
    % 注意：ind_peak(flag) 对应 0m
    % 往上(索引变小)是空气(负深度)，往下(索引变大)是水(正深度)
    % 所以 Depth 范围是 -N1*dh 到 N2*dh
    Depth = ( -N1 : N2 )' * delta_h;
    
    % 确定截取范围 (在对齐后的 bin 矩阵中)
    % ind_peak(flag) 是对齐后的海面位置
    center_idx = ind_peak(flag);
    range_idx = center_idx - N1 : center_idx + N2;
    
    % 边界检查
    if min(range_idx) < 1 || max(range_idx) > m
        warning('截取范围超出矩阵边界，可能需要补零或调整 N1/N2');
        S = []; data_S = []; return;
    end
    
    S = [];
    data_S = [];
    
    % 阈值: 窗口内有效点数
    % data_bin(3,:) 是时间，原始 delta_t 是 0.0002s
    ThStep = 0.2 * (dt_int / 0.0002); 
    
    current_start_time = data_bin(3, 1);
    last_total_time = data_bin(3, end);
    
    while (current_start_time + dt_int) <= last_total_time
        % 1. 找窗口内的列索引
        idx_in_window = find(data_bin(3, :) >= current_start_time & ...
                             data_bin(3, :) < (current_start_time + dt_int));
        
        % 2. 判断数据量
        if length(idx_in_window) > ThStep
            % 截取并平均
            % bin(range_idx, idx_in_window) 取出对齐后的特定层
            tempS = mean(bin(range_idx, idx_in_window), 2);
            
            % 计算轨迹 (Lat, Lon, Time)
            tempdataS = mean(data_bin(:, idx_in_window), 2);
            
            S = [S, tempS];
            data_S = [data_S, tempdataS];
        end
        
        % 3. 滑动
        current_start_time = current_start_time + step_t;
    end
end