function [S, Depth, data_S] = f_profile_processing(file_choose, Dayflag, savepath, lat_range, window_m)
%% 1. 数据加载与预处理
% Dayflag 目前未使用，保留接口兼容
% lat_range: 可选，[min_lat, max_lat]，只保留该纬度范围内的廓线
% window_m: 可选，沿轨统计窗口长度 [m]，默认100m
if nargin < 3
    savepath = '';
end
if nargin < 4 || isempty(lat_range)
    lat_range = [-inf, inf];
end
if nargin < 5 || isempty(window_m)
    window_m = 100;
end

data_struct = load(file_choose);
if isfield(data_struct, 'initial_t')
    initial_t = data_struct.initial_t;
    initial_h = data_struct.initial_h;
    initial_lat = data_struct.initial_lat;
    initial_lon = data_struct.initial_lon;
else
    error('缺少必要的变量 (initial_t/h/lat/lon)');
end

% 初始化
S = [];
Depth = [];
data_S = [];

% 保留原始数据，不做经纬度范围过滤
if isempty(initial_h)
    warning('数据为空');
    return;
end

% --- 相对时间与高度筛选（分段局部中位数）---
% 轨道可能跨越冰盖（~2000m）和洋面（~0m），全局中位数会偏向
% 某一端导致另一端数据被误删。改为按时段分别计算中位数。
t_start = initial_t(1);
rel_t = initial_t - t_start;
threshold_h = [200 600];

% 分段时间窗口（秒），默认10秒一段
seg_dt = 10;
t_total = max(rel_t);
if t_total <= seg_dt
    seg_edges = [0, t_total + 1];
else
    seg_edges = 0:seg_dt:(t_total + seg_dt);
end
seg_id = discretize(rel_t, seg_edges);

filter = false(size(initial_h));
n_segs = max(seg_id);
for s = 1:n_segs
    in_seg = (seg_id == s);
    n_seg = sum(in_seg);
    if n_seg < 50
        continue;  % 点数过少，跳过该段
    end
    local_center = median(initial_h(in_seg));
    if ~isfinite(local_center)
        continue;
    end
    seg_h = initial_h(in_seg);
    filter(in_seg) = seg_h >= local_center - threshold_h(1) ...
                   & seg_h <= local_center + threshold_h(2);
end

initial_h = initial_h(filter);
rel_t = rel_t(filter);
initial_lat = initial_lat(filter);
initial_lon = initial_lon(filter);

if isempty(initial_h)
    warning('高度筛选后数据为空');
    return;
end

%% 2. 阶段一：拟合对齐 (Fitting Alignment - rloess)

% 2.1 寻找粗表面（分块 histogram，避免一次性爆内存）
dt_fine = 0.001;
dh_fine = 0.01;
max_t_bins_per_chunk = 5000;  % 单块最多 2 万个时间 bin，可按内存继续调小

[t_edges, t_centers_fine, raw_surface, valid_mask] = ...
    local_surface_histogram_chunked(rel_t, initial_h, dt_fine, dh_fine, max_t_bins_per_chunk);

if ~any(valid_mask)
    warning('未能从 histogram 中找到有效海面');
    return;
end

% 对空 bin 先插值，再做去异常和平滑
surface_series = raw_surface;
surface_series(~valid_mask) = NaN;
surface_series = fillmissing(surface_series, 'linear', 'EndValues', 'nearest');
[h_despiked, ~] = filloutliers(surface_series, 'linear', 'movmedian', 50);
h_filtered = smoothdata(h_despiked, 'rloess', 150);

% ===================== 海面拟合对比绘图 =====================
hSurfaceFig = figure('Visible', 'off', 'Name', 'Sea Surface Detection', 'Color', 'w');
plot(t_centers_fine(valid_mask), raw_surface(valid_mask), '.', 'Color', [0.7 0.7 0.7], 'DisplayName', 'Initial Surface (Histogram Max)');
hold on;
plot(t_centers_fine, h_filtered, 'r-', 'LineWidth', 1.5, 'DisplayName', 'Filtered Surface (rloess)');
xlabel('Relative Time (s)'); ylabel('Height (m)');
title('Comparison of Initial and Filtered Sea Surface');
legend('Location', 'best'); grid on;
if ~isempty(savepath)
    [~, name_stem, ~] = fileparts(file_choose);
    if ~exist(savepath, 'dir')
        mkdir(savepath);
    end
    saveas(hSurfaceFig, fullfile(savepath, [name_stem, '_surface_fit.jpg']));
end
close(hSurfaceFig);
% ===========================================================

% 2.3 计算对齐后高度
num_bins = length(t_edges) - 1;
surface_lookup = NaN(num_bins, 1);
surface_lookup(:) = h_filtered(:);
photon_bin_idx = discretize(rel_t, t_edges);
valid_p_idx = ~isnan(photon_bin_idx);
photon_surface = NaN(size(initial_h));
photon_surface(valid_p_idx) = surface_lookup(photon_bin_idx(valid_p_idx));

% *** 核心：海面变为了 0 ***
h_aligned = initial_h - photon_surface;

% 剔除无效点
final_mask = ~isnan(h_aligned);
rel_t_final = rel_t(final_mask);
h_aligned_final = h_aligned(final_mask);
lat_final = initial_lat(final_mask);
lon_final = initial_lon(final_mask);

if isempty(rel_t_final)
    warning('海面对齐后无有效光子');
    return;
end

%% 3. 阶段二：构建严格网格 (Strict Binning)

delta_h = 0.15;
delta_t_prof = 0.001;
N1 = 100;
N2 = 140;
Depth = (-N1:N2) * delta_h;

depth_val = -h_aligned_final;
z_edges = [Depth - delta_h/2, Depth(end) + delta_h/2];
t_edges_prof = min(rel_t_final) : delta_t_prof : max(rel_t_final) + delta_t_prof;
t_centers = t_edges_prof(1:end-1) + delta_t_prof/2;

[N_prof, ~, ~] = histcounts2(rel_t_final, depth_val, t_edges_prof, z_edges);
bin_matrix = N_prof';

[time_idx] = discretize(rel_t_final, t_edges_prof);
num_t = length(t_centers);
valid_t = ~isnan(time_idx);
bin_lat = accumarray(time_idx(valid_t), lat_final(valid_t), [num_t, 1], @mean, NaN);
bin_lon = accumarray(time_idx(valid_t), lon_final(valid_t), [num_t, 1], @mean, NaN);
bin_lat = fillmissing(bin_lat, 'linear');
bin_lon = fillmissing(bin_lon, 'linear');
base_data_bin = [bin_lat'; bin_lon'; t_centers];

%% 4. 阶段三：廓线生成 (滑动窗口，基于沿轨距离)
% 计算卫星地面速度 [m/s]，用于将距离窗口转换为时间窗口
n_velo = min(1000, length(bin_lat));
if n_velo >= 2
    d = zeros(1, n_velo-1);
    for k = 1:n_velo-1
        d(k) = haversine(bin_lat(k), bin_lon(k), bin_lat(k+1), bin_lon(k+1));
    end
    ground_speed = sum(d) / (t_centers(n_velo) - t_centers(1));
else
    ground_speed = 6900;
end
dt_int = window_m / ground_speed;
step_t = dt_int / 2;
fprintf('[f_profile_processing] 沿轨窗口=%.0fm, 地速=%.0fm/s, dt_int=%.3fs, step_t=%.3fs\n', ...
    window_m, ground_speed, dt_int, step_t);

ThStep = max(3, 0.2 * (dt_int / delta_t_prof));
current_time = base_data_bin(3, 1);
end_time = base_data_bin(3, end);
while (current_time + dt_int) <= end_time
    idx_win = find(base_data_bin(3, :) >= current_time & ...
                   base_data_bin(3, :) < (current_time + dt_int));

    if length(idx_win) > ThStep
        win_block = bin_matrix(:, idx_win);
        tempS = mean(win_block, 2);

        mean_lat = mean(base_data_bin(1, idx_win));
        mean_lon = mean(base_data_bin(2, idx_win));
        total_ph = sum(win_block(:));

        S = [S, tempS]; %#ok<AGROW>
        data_S = [data_S, [mean_lat; mean_lon; total_ph]]; %#ok<AGROW>
    end
    current_time = current_time + step_t;
end
if isempty(S)
    warning('S为空');
    return;
end

%% 5. 阶段四：后处理 (平滑与噪声计算)
St = S;
start_smooth = find(Depth >= 5, 1);
if ~isempty(start_smooth)
    Si = S(start_smooth:end, :);
    window_half_width = 1.5;
    window_bins = round((2 * window_half_width) / delta_h);
    Si = smoothdata(Si, 1, 'sgolay', window_bins, 'Degree', 2);
    S(start_smooth:end, :) = Si;
end

ind_noise_rows = find(Depth >= -15 & Depth <= -5);
if ~isempty(ind_noise_rows)
    sub_S = S(ind_noise_rows, :);
    sub_S(isnan(sub_S)) = 0;
    noise = mean(sub_S, 1);
else
    noise = zeros(1, size(S, 2));
    warning('未找到 -15m 到 -5m 范围内的数据用于计算噪声');
end
S = S - noise;
S(S <= 0) = NaN;

%% 6. 生成Np_vec和phi_vec
% Np_vec：每廓线雪面光子数
% phi_vec：每廓线沿轨坡度（默认0度）
Np_vec = zeros(1, size(S, 2));
phi_vec = zeros(1, size(S, 2));

% ICESat-2发射频率为10kHz，即每0.001秒10发
% 需要将光子数从每0.001秒转换为每shot
shots_per_bin = 10;

for i = 1:size(S, 2)
    % 找到雪面附近的信号峰值（深度0附近）
    ind_surface = find(Depth >= -0.5 & Depth <= 0.5);
    if ~isempty(ind_surface)
        Np_vec(i) = max(S(ind_surface, i)) / shots_per_bin;
    else
        Np_vec(i) = max(S(:, i)) / shots_per_bin;
    end
    % 坡度默认为0度
    phi_vec(i) = 0;
end

%% 6.5 纬度过滤
if ~all(isfinite(lat_range)) || lat_range(1) > -inf || lat_range(2) < inf
    keep = data_S(1,:) >= lat_range(1) & data_S(1,:) <= lat_range(2);
    n_total = size(S, 2);
    S       = S(:, keep);
    St      = St(:, keep);
    data_S  = data_S(:, keep);
    noise   = noise(keep);
    Np_vec  = Np_vec(keep);
    phi_vec = phi_vec(keep);
    if isempty(S)
        warning('纬度过滤 [%.2f, %.2f] 后无廓线', lat_range(1), lat_range(2));
        return;
    end
    fprintf('纬度过滤 [%.2f, %.2f]: %d/%d 条廓线保留\n', lat_range(1), lat_range(2), sum(keep), n_total);
end

%% 7. 绘图与保存
if ~isempty(savepath)
    if ~exist(savepath, 'dir')
        mkdir(savepath);
    end
    [~, name_stem, ~] = fileparts(file_choose);

    hFig = figure('Visible', 'off');
    set(hFig, 'units', 'normalized', 'position', [0.4 0.2 0.4 0.8]);

    num_plot = min(100, size(S, 2));
    hold on;
    for i = 1:num_plot
        curr_prof = S(:, i);
        curr_max = max(curr_prof);
        if curr_max > 0
            curr_prof = curr_prof / curr_max;
        end
        curr_prof(curr_prof <= 0) = NaN;
        semilogy(Depth, curr_prof, 'LineWidth', 0.8, 'Color', 'b');
    end
    for i = 1:num_plot
        curr_prof = St(:, i);
        curr_max = max(curr_prof);
        if curr_max > 0
            curr_prof = curr_prof / curr_max;
        end
        curr_prof(curr_prof <= 0) = NaN;
        semilogy(Depth, curr_prof, 'LineWidth', 0.8, 'Color', 'r');
    end
    set(gca, 'YScale', 'log');
    grid on;
    set(gca, 'YMinorGrid', 'on', 'YMinorTick', 'on');
    xlabel('Depth (m)');
    ylabel('Normalized Signal (Log Scale)');
    title(['Profile: ' name_stem], 'Interpreter', 'none');
    xlim([-15 21]);
    ylim([1e-4 1.5]);
    saveas(hFig, fullfile(savepath, [name_stem, '_profile.jpg']));
    save(fullfile(savepath, [name_stem, '_profile.mat']), 'S', 'Depth', 'data_S', 'noise', 'Np_vec', 'phi_vec');
end
end

function [t_edges, t_centers, raw_surface, valid_mask] = local_surface_histogram_chunked(rel_t, initial_h, dt_fine, dh_fine, max_t_bins_per_chunk)
% 按时间 bin 分块做 2D histogram，避免一次性申请超大矩阵

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

function d = haversine(lat1, lon1, lat2, lon2)
R = 6371000;
dlat = deg2rad(lat2 - lat1);
dlon = deg2rad(lon2 - lon1);
a = sin(dlat/2)^2 + cos(deg2rad(lat1)) * cos(deg2rad(lat2)) * sin(dlon/2)^2;
c = 2 * atan2(sqrt(a), sqrt(1-a));
d = R * c;
end
