function [file_profile] = f_ocean_choose(file_choose, Dayflag, savepath)
%% Memory-safe ocean profile extraction
if nargin < 3
    savepath = '';
end
load(file_choose);
file_profile = [];
Depth = [];
S = [];
data_S = [];
a = dir(file_choose);
lat_min =  80;
lat_max = 90;
lon_max = -100.63;
lon_min = -179.36;
t_start = initial_t(1);
rel_t= initial_t - t_start;
ind_noise = (initial_lat < lat_min) | (initial_lat >= lat_max) | ...
            (initial_lon < lon_min) | (initial_lon >= lon_max); 
initial_h(ind_noise) = [];
initial_t(ind_noise) = [];
initial_lon(ind_noise) = [];
back_ground_save(ind_noise) = [];
initial_lat(ind_noise) = [];
keep_mask = true(size(initial_h));

% Define the time windows (rounded to the nearest second)
t_floor = floor(initial_t);
unique_seconds = unique(t_floor);

for i = 1:length(unique_seconds)
    % Find indices for the current 1-second window
    idx = find(t_floor == unique_seconds(i));
    
    % Get the data for this segment
    segment_h = initial_h(idx);
    
    % Condition 1: Photon count (number of points) < 30,000
    cond1 = length(segment_h) < 30000;
    
    % Condition 2: Vertical fluctuation > 100m
    % (Checking the range: max - min)
    cond2 = (max(segment_h) - min(segment_h)) > 100;
    
    % If either condition is met, mark the whole segment for deletion
    if cond1 
        keep_mask(idx) = false;
    end
end

% % --- 1. 统计每一段（1秒）的数据 ---
%     t_floor = floor(initial_t);
%     unique_seconds = unique(t_floor);
% 
%     photon_counts = zeros(length(unique_seconds), 1);
%     mean_lats = zeros(length(unique_seconds), 1);
%     h_std = zeros(length(unique_seconds), 1); % 用于记录波动情况
% 
%     for j = 1:length(unique_seconds)
%         idx = (t_floor == unique_seconds(j));
%         segment_h = initial_h(idx);
% 
%         photon_counts(j) = sum(idx);
%         if photon_counts(j) > 0
%             mean_lats(j) = mean(initial_lat(idx));
%             h_std(j) = max(segment_h) - min(segment_h); % 计算波动范围
%         end
%     end
% 
%     % --- 2. 绘图监控 ---
%     fig = figure('Color', 'w');
% 
%     % 子图1：光子数分布
%     subplot(2,1,1);
%     bar(mean_lats, photon_counts, 'FaceColor', [0.2 0.6 0.8]);
%     hold on;
%     % 画出 30,000 的红色阈值线
%     yline(30000, 'r--', 'Threshold: 30k', 'LineWidth', 1.5);
%     ylabel('Photon Count');
% 
%     grid on;
% Apply the mask to all variables
initial_h = initial_h(keep_mask);
initial_t = initial_t(keep_mask);
initial_lon = initial_lon(keep_mask);
initial_lat = initial_lat(keep_mask);




% 
% ind_noise = (initial_lat >82.6 ) & (initial_lat <84.4);
% 
% initial_h(ind_noise) = [];
% initial_t(ind_noise) = [];
% initial_lon(ind_noise) = [];
% initial_lat(ind_noise) = [];


if isempty(initial_h)
    return
end

% 高度粗筛选
center = median(initial_h);
threshold_h = [50 50];
filter = find(initial_h >= center - threshold_h(1) & initial_h <= center + threshold_h(2));

% --- 强制维度对齐修复 ---
initial_h = initial_h(filter);
initial_t = initial_t(filter);
initial_lon = initial_lon(filter);
initial_lat = initial_lat(filter);

% 确保所有列向量长度一致
n_photons = length(initial_h);
initial_t = initial_t(1:n_photons);
initial_lon = initial_lon(1:n_photons);
initial_lat = initial_lat(1:n_photons);

rel_t = initial_t - initial_t(1);

if isempty(initial_h), return; end

%% 2. 表面提取与迭代差分 (绝对不平滑)
delta_h = 0.15;
delta_t = 0.001;
issue.max_t_bins_per_chunk = 5000;
[t_edges_surface, ~, raw_surface, valid_mask] = ...
    local_surface_histogram_chunked(rel_t, initial_h, delta_t, delta_h, issue.max_t_bins_per_chunk);

if ~any(valid_mask), return; end

ind_peak0 = raw_surface;
ind_peak0(~valid_mask) = NaN; 
ind_peak0(ind_peak0 < (center - 20)) = NaN; 

% 迭代差分
Threshold_diff = 0.5; 
x0 = 1:length(ind_peak0);
x_old = x0(~isnan(ind_peak0)); x_old = x_old(:)'; 
ind_peak_old = ind_peak0(~isnan(ind_peak0)); ind_peak_old = ind_peak_old(:)';

x_new = []; diffTime = 0;
while(length(x_new) ~= length(x_old) && diffTime < 2)
    if(diffTime > 0)
        x_old = x_new; 
        ind_peak_old = ind_peak_new; 
    end
    if length(ind_peak_old) < 2, break; end
    df = abs(diff(ind_peak_old));
    ind_remain = find(df <= Threshold_diff);
    ind_remain = ind_remain(:)';
    x_new = x_old([1, ind_remain + 1]);
    ind_peak_new = ind_peak_old([1, ind_remain + 1]);
    diffTime = diffTime + 1;
end

ind_peak_step2 = nan(size(ind_peak0));
ind_peak_step2(x_new) = ind_peak_new;

%% 3. 连续性筛选 (dt * 7km)
N_km = 0; % 阈值 2km
nan_mask = isnan(ind_peak_step2);
nan_mask = nan_mask(:)'; 

starts = find(diff([0, nan_mask]) == 1); 
ends = find(diff([nan_mask, 0]) == -1);   
t_bin_centers = (t_edges_surface(1:end-1) + t_edges_surface(2:end)) / 2;

ind_peak_final = ind_peak_step2; 

% 遍历缺口，如果 dt*7 > N_km，保持 NaN 状态，不插值
for i = 1:length(starts)
    s_idx = starts(i); e_idx = ends(i);
    if s_idx > 1 && e_idx < length(ind_peak_final)
        dt = t_bin_centers(e_idx + 1) - t_bin_centers(s_idx - 1);
        % gap_dist = dt * 7; 
        % 这里如果不连续且不插值，对应的光子会被 final_mask 直接物理删除
    end
end

%% 4. 数据对齐与变量输出 (锁定变量名，无平滑)
photon_bin_idx = discretize(rel_t, t_edges_surface);
valid_map = ~isnan(photon_bin_idx);
photon_surface = NaN(size(initial_h));

% 直接映射，绝对不平滑
photon_surface(valid_map) = ind_peak_final(photon_bin_idx(valid_map));

h_aligned = initial_h - photon_surface;
final_mask = ~isnan(h_aligned); 

% 严格保留变量名
rel_t_final = rel_t(final_mask);
h_final = initial_h(final_mask);
lat_final = initial_lat(final_mask);
lon_final = initial_lon(final_mask);
surf_final = photon_surface(final_mask);
h_aligned_final = h_aligned(final_mask);

if isempty(rel_t_final), return; end

%% 5. 绘图
figure('Color','w');
scatter(lat_final, h_final, 1, [0.8 0.8 0.8], 'filled'); hold on;
plot(lat_final, surf_final, 'r.', 'MarkerSize', 3); 
title('Final: Iterative Diff & Distance Gap Removal');
ylim([min(surf_final)-10, max(surf_final)+10]);

%% Build compact aligned bin matrix only in the profile depth range
delta_t_prof = 0.001;
dt = 2;
step_t = 0.5;
N1 = 100;
N2 = 200;
Depth = (-N1:N2)' * delta_h;

depth_val = -h_aligned_final;
z_edges = [Depth' - delta_h/2, Depth(end) + delta_h/2];
t_edges_prof = min(rel_t_final) : delta_t_prof : (max(rel_t_final) + delta_t_prof);
t_centers = t_edges_prof(1:end-1) + delta_t_prof/2;

if isempty(t_centers)
    return
end

% ... 前面 histcounts2 等逻辑保持不变 ...

N_prof = histcounts2(rel_t_final, depth_val, t_edges_prof, z_edges);
bin_matrix = N_prof'; % 此时 bin_matrix 每一列对应一个 t_center
time_idx = discretize(rel_t_final, t_edges_prof, 'IncludedEdge', 'right'); % 修正：应与 t_edges_prof 对应
num_t = numel(t_centers);


bin_lat = accumarray(time_idx, lat_final, [num_t, 1], @mean, NaN);
bin_lon = accumarray(time_idx, lon_final, [num_t, 1], @mean, NaN);
bin_lat = fillmissing(bin_lat, 'linear', 'EndValues', 'nearest');
bin_lon = fillmissing(bin_lon, 'linear', 'EndValues', 'nearest');

base_data_bin = [bin_lat'; bin_lon'; t_centers];

%% --- 核心修改：物理删除全 0 列 ---
% 统计每列光子总数，只保留总数 > 0 的列
has_photons = sum(bin_matrix, 1) > 0;
bin_matrix = bin_matrix(:, has_photons);
base_data_bin = base_data_bin(:, has_photons);
bin_lat = bin_lat(has_photons);
bin_lon = bin_lon(has_photons);


if isempty(base_data_bin)
    return; % 如果整段都没光子，直接跳出
end
%% ------------------------------

ThStep = 0.2* (dt / delta_t_prof);
current_start_time = base_data_bin(3, 1);
last_total_time = base_data_bin(3, end);

while (current_start_time + dt) <= last_total_time
    % 在已经剔除 0 列的 base_data_bin 中寻找符合当前时间窗的索引
    idx_in_window = find(base_data_bin(3, :) >= current_start_time & ...
                         base_data_bin(3, :) < (current_start_time + dt));
     
    
    if length(idx_in_window) > ThStep
        % 这里的平均值将只基于有数据的列进行计算
        tempS = mean(bin_matrix(:, idx_in_window), 2);
        tempdataS = mean(base_data_bin(:, idx_in_window), 2);
        S = [S, tempS]; 
        data_S = [data_S, tempdataS]; 
    end
    current_start_time = current_start_time + step_t;
end

if isempty(S)
    return
end

St = S;

start = find(Depth >= 3, 1);
if ~isempty(start)
    Si = S(start:end, :);
    S(start:end, :) = Si;
end
start_smooth = find(Depth >= 8, 1);
if ~isempty(start_smooth)
    Si = S(start_smooth:end, :);
    window_half_width = 1.5;
    window_bins = round((2 * window_half_width) / delta_h);
    Si = smoothdata(Si, 1, 'sgolay', window_bins, 'Degree', 2);
    S(start_smooth:end, :) = Si;
end

ind_noise = find(Depth >= -15 & Depth <= -8);
sub_S = S(ind_noise, :);
sub_S(isnan(sub_S)) = 0;
noise = mean(sub_S, 1);
S = S-noise;
S(S <= 0) = nan;
f = figure('Visible','on');
set(f,'units','normalized','position',[0.4 0.2 0.4 0.6]);
set(gca,'position',[0.15 0.15 0.7 0.7]);
hold on;
for i = 1:size(St, 2)
    semilogy(Depth, St(:, i)./max( St(:, i)), 'LineWidth', 0.5, 'color', 'r');
    semilogy(Depth, S(:, i)./max( S(:, i)), 'LineWidth', 0.5, 'color', 'b');
end
 set(gca, 'YScale', 'log');
xlabel('Depth','FontSize',10);
ylabel('Normalized Signal','FontSize',10);
ylim([10^-4 10^1]);

jpg_profile = [savepath, a.name(1:4), 'profile_', a.name(12:end-3), 'jpg'];
saveas(gcf, jpg_profile);
file_profile = [savepath, a.name(1:4), 'profile_', a.name(12:end-3), 'mat'];
save(file_profile, 'S', 'Depth', 'data_S', 'noise');

end

function [t_edges, t_centers, raw_surface, valid_mask] = local_surface_histogram_chunked(rel_t, initial_h, dt_fine, dh_fine, max_t_bins_per_chunk)
t_min = min(rel_t);
t_max = max(rel_t);
t_edges = t_min : dt_fine : (t_max + dt_fine);
num_t_bins = numel(t_edges) - 1;
t_centers = t_edges(1:end-1) + dt_fine / 2;

raw_surface = NaN(num_t_bins, 1);
valid_mask = false(num_t_bins, 1);

chunk_start_bin = 1;
while chunk_start_bin <= num_t_bins
    chunk_end_bin = min(chunk_start_bin + max_t_bins_per_chunk - 1, num_t_bins);
    t_left = t_edges(chunk_start_bin);
    t_right = t_edges(chunk_end_bin + 1);

    if chunk_end_bin < num_t_bins
        idx_chunk = rel_t >= t_left & rel_t < t_right;
    else
        idx_chunk = rel_t >= t_left & rel_t <= t_right;
    end

    if any(idx_chunk)
        h_chunk = initial_h(idx_chunk);
        t_chunk = rel_t(idx_chunk);

        h_min = min(h_chunk);
        h_max = max(h_chunk);
        if h_min == h_max
            raw_surface(chunk_start_bin:chunk_end_bin) = h_min;
            valid_mask(chunk_start_bin:chunk_end_bin) = true;
        else
            h_edges = h_min : dh_fine : (h_max + dh_fine);
            N = histcounts2(t_chunk, h_chunk, t_edges(chunk_start_bin:chunk_end_bin+1), h_edges);
            [max_intensity, max_idx] = max(N, [], 2);
            h_centers = h_edges(1:end-1) + dh_fine / 2;

            local_surface = NaN(chunk_end_bin - chunk_start_bin + 1, 1);
            local_valid = max_intensity > 0;
            local_surface(local_valid) = h_centers(max_idx(local_valid))';

            raw_surface(chunk_start_bin:chunk_end_bin) = local_surface;
            valid_mask(chunk_start_bin:chunk_end_bin) = local_valid;
        end
    end

    chunk_start_bin = chunk_end_bin + 1;
end
end
