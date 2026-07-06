function Ta = f_Ta_from_atl09_v2(atl09_file, lat_query, lon_query, varargin)
% F_TA_FROM_ATL09_V2 从ATL09文件中提取大气光学厚度并计算Ta
%
% 输入:
%   atl09_file - ATL09文件路径（.mat或.h5）或文件路径单元格数组
%   lat_query - 查询点纬度
%   lon_query - 查询点经度
%   varargin - 可选参数
%
% 输出:
%   Ta - 大气透过率

%% 解析输入参数
p = inputParser;
addRequired(p, 'atl09_file', @(x) ischar(x) || iscell(x));
addRequired(p, 'lat_query', @isnumeric);
addRequired(p, 'lon_query', @isnumeric);
addParameter(p, 'Ta_dataset', '', @ischar);
addParameter(p, 'gt', '', @ischar);
addParameter(p, 'radius_m', 5000, @isnumeric); % 默认5km
addParameter(p, 'Ta_is_twoway', false, @islogical);
addParameter(p, 'mat_lat_var', 'initial_lat', @ischar);
addParameter(p, 'mat_lon_var', 'initial_lon', @ischar);
addParameter(p, 'mat_Ta_var', 'column_od_asr', @ischar);
addParameter(p, 'verbose', true, @islogical);
parse(p, atl09_file, lat_query, lon_query, varargin{:});

radius_m = p.Results.radius_m;
verbose = p.Results.verbose;
mat_lat_var = p.Results.mat_lat_var;
mat_lon_var = p.Results.mat_lon_var;
mat_Ta_var = p.Results.mat_Ta_var;

if verbose
    fprintf('[f_Ta_from_atl09_v2] === 诊断开始 ===\n');
    fprintf('[f_Ta_from_atl09_v2] mat_lat_var=[%s], mat_lon_var=[%s], mat_Ta_var=[%s]\n', mat_lat_var, mat_lon_var, mat_Ta_var);
end

%% 加载ATL09数据
all_lat = [];
all_lon = [];
all_od = [];

% 处理多个文件
if iscell(atl09_file)
    file_list = atl09_file;
else
    file_list = {atl09_file};
end

if verbose
    fprintf('[f_Ta_from_atl09_v2] 加载%d个ATL09文件\n', length(file_list));
end

try
    for i = 1:length(file_list)
        current_file = file_list{i};
        if verbose
            fprintf('[f_Ta_from_atl09_v2] 加载文件 %d/%d: %s\n', i, length(file_list), current_file);
        end
        
        atl09_data = load(current_file);
        
        % 检查变量是否存在
        has_lat = isfield(atl09_data, mat_lat_var);
        has_lon = isfield(atl09_data, mat_lon_var);
        has_Ta = isfield(atl09_data, mat_Ta_var);
        
        if has_lat && has_lon && has_Ta
            lat = atl09_data.(mat_lat_var);
            lon = atl09_data.(mat_lon_var);
            od = atl09_data.(mat_Ta_var);
            
            % 确保数据是列向量
            if isrow(lat); lat = lat'; end
            if isrow(lon); lon = lon'; end
            if isrow(od); od = od'; end
            
            % 合并数据
            all_lat = [all_lat; lat];
            all_lon = [all_lon; lon];
            all_od = [all_od; od];
            
            if verbose
                fprintf('[f_Ta_from_atl09_v2] 文件%d: lat(%dx1), lon(%dx1), Ta/OD(%dx1)\n', i, length(lat), length(lon), length(od));
                fprintf('[f_Ta_from_atl09_v2] 范围: lat[%.4f,%.4f], lon[%.4f,%.4f]\n', min(lat), max(lat), min(lon), max(lon));
            end
        else
            warning('ATL09文件 %s 中缺少必要的变量', current_file);
        end
    end
    
    if isempty(all_lat)
        error('所有ATL09文件中都缺少必要的变量');
    end
    
    % 使用合并后的数据
    lat = all_lat;
    lon = all_lon;
    od = all_od;
    
    if verbose
        fprintf('[f_Ta_from_atl09_v2] 合并后: lat(%dx1), lon(%dx1), Ta/OD(%dx1)\n', length(lat), length(lon), length(od));
        fprintf('[f_Ta_from_atl09_v2] 合并后范围: lat[%.4f,%.4f], lon[%.4f,%.4f]\n', min(lat), max(lat), min(lon), max(lon));
    end
    
    % 处理数据
    finite_data = od(isfinite(od));
    if isempty(finite_data)
        error('ATL09数据中无有效OD值');
    end
    
    if verbose
        fprintf('[f_Ta_from_atl09_v2] Ta/OD原始值: 范围[%.6f,%.6f], 中位数=%.6f, 有效数=%d/%d\n', ...
            min(finite_data), max(finite_data), median(finite_data), length(finite_data), length(od));
    end
    
    %% 强制将column_od_asr视为OD，直接进行 Ta=exp(-OD) 转换
    % 注释掉自动判断逻辑，直接做 Ta = exp(-od)
    if verbose
        fprintf('[f_Ta_from_atl09_v2] ★★★ 强制将column_od_asr视为OD，进行 Ta=exp(-OD) 转换 ★★★\n');
    end
    Ta_atl09 = exp(-od);
    
    if verbose
        finite_Ta = Ta_atl09(isfinite(Ta_atl09));
        fprintf('[f_Ta_from_atl09_v2] 转换后Ta: 范围[%.6f,%.6f], 中位数=%.6f\n', ...
            min(finite_Ta), max(finite_Ta), median(finite_Ta));
    end
    
    %% 地理匹配
    if verbose
        fprintf('[f_Ta_from_atl09_v2] 开始ECEF匹配: %d个segment → %d个query, radius=%.1fm\n', ...
            length(lat), length(lat_query), radius_m);
    end
    
    % 调用匹配函数
    Ta = matchTaNearest(lat, lon, Ta_atl09, lat_query, lon_query, radius_m, verbose);
    
    if verbose
        valid_Ta = Ta(isfinite(Ta));
        fprintf('[f_Ta_from_atl09_v2] === 诊断结束 ===\n');
        if isempty(valid_Ta)
            fprintf('[f_Ta_from_atl09_v2] 匹配结果: 0/%d 有效, Ta范围[NaN,NaN], med=NaN\n', length(Ta));
            fprintf('[f_Ta_from_atl09_v2] ⚠ 全部返回NaN! 将在上层回退为0.9 ⚠\n');
        else
            fprintf('[f_Ta_from_atl09_v2] 匹配结果: %d/%d 有效, Ta范围[%.4f,%.4f], med=%.4f\n', ...
                length(valid_Ta), length(Ta), min(valid_Ta), max(valid_Ta), median(valid_Ta));
        end
    end
    
catch ME
    if verbose
        fprintf('[f_Ta_from_atl09_v2] 错误: %s\n', ME.message);
        fprintf('[f_Ta_from_atl09_v2] ⚠ 返回NaN，将在上层回退为0.9 ⚠\n');
    end
    Ta = NaN(size(lat_query));
end

end

function Ta = matchTaNearest(lat_seg, lon_seg, Ta_seg, lat_query, lon_query, radius_m, verbose)
% MATCHTANEarest 匹配最近的ATL09点

Ta = NaN(size(lat_query));

if isempty(lat_seg) || isempty(lat_query)
    return;
end

% 转换为ECEF坐标
[x_seg, y_seg, z_seg] = latlon2ecef(lat_seg, lon_seg);

min_distances = zeros(size(lat_query));

for i = 1:length(lat_query)
    [x_q, y_q, z_q] = latlon2ecef(lat_query(i), lon_query(i));
    
    % 计算距离
    dx = x_seg - x_q;
    dy = y_seg - y_q;
    dz = z_seg - z_q;
    distances = sqrt(dx.^2 + dy.^2 + dz.^2);
    
    % 找到最近的点
    [min_dist, min_idx] = min(distances);
    min_distances(i) = min_dist;
    
    % 检查是否在半径范围内
    if min_dist <= radius_m
        Ta(i) = Ta_seg(min_idx);
    end
end

if verbose
    fprintf('[matchTaNearest] ⚡ 诊断: query lat范围[%.4f,%.4f], lon范围[%.4f,%.4f]\n', ...
        min(lat_query), max(lat_query), min(lon_query), max(lon_query));
    fprintf('[matchTaNearest] ⚡ 诊断: segment lat范围[%.4f,%.4f], lon范围[%.4f,%.4f]\n', ...
        min(lat_seg), max(lat_seg), min(lon_seg), max(lon_seg));
    fprintf('[matchTaNearest] ⚡ 匹配统计: %d/%d 成功, 最近距离范围[%.1f, %.1f]m\n', ...
        sum(isfinite(Ta)), length(Ta), min(min_distances), max(min_distances));
    fprintf('[matchTaNearest] ⚡ 前5个query点的最近距离(m): ');
    for i = 1:min(5, length(min_distances))
        fprintf('%.1f ', min_distances(i));
    end
    fprintf('\n');
end

end

function [x, y, z] = latlon2ecef(lat, lon, alt)
% LATLON2ECEF 将经纬度转换为ECEF坐标

if nargin < 3
    alt = zeros(size(lat));
end

% 地球参数
a = 6378137; % 赤道半径
f = 1/298.257223563; % 扁率
e2 = 2*f - f^2; % 第一偏心率平方

lat_rad = deg2rad(lat);
lon_rad = deg2rad(lon);

N = a ./ sqrt(1 - e2 * sin(lat_rad).^2);

x = (N + alt) .* cos(lat_rad) .* cos(lon_rad);
y = (N + alt) .* cos(lat_rad) .* sin(lon_rad);
z = (N .* (1 - e2) + alt) .* sin(lat_rad);

end