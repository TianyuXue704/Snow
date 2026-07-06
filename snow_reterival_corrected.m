function [snow, lat, snowA, info_out] = snow_reterival_corrected(file, varargin)

%% ===== 输入解析 =====
p = inputParser;
addRequired(p,  'file',        @ischar);
addParameter(p, 'Ta',          0.9,     @isnumeric);
% 若提供 ATL09 文件，则从 ATL09 覆盖输入的粗略 Ta（用于逐廓线/逐段）
addParameter(p, 'atl09_file',        '',     @(x) ischar(x) || iscell(x));
addParameter(p, 'atl09_Ta_dataset',  '',     @ischar); % 可选：ATL09 内 Ta 数据集路径
addParameter(p, 'atl09_gt',          '',     @ischar); % 可选：例如 gt1l/gt1r（不填则自动探测）
addParameter(p, 'atl09_radius_m',    50,    @isnumeric); % 最近邻匹配半径（米），超出则回退 Ta
addParameter(p, 'Ta_is_twoway',    false,  @islogical); % 若 ATL09 里是双程透过率，则取 sqrt 转成单程
% 若 atl09_file 是 .mat 并且自动识别失败，可手动指定变量名（完全绕过启发式）
addParameter(p, 'atl09_mat_lat_var',  '', @ischar);
addParameter(p, 'atl09_mat_lon_var',  '', @ischar);
addParameter(p, 'atl09_mat_Ta_var',   '', @ischar);
addParameter(p, 'CF',          0.52,     @isnumeric);
addParameter(p, 'beam',        'strong', @ischar);
addParameter(p, 'verbose',     true,     @islogical);
addParameter(p, 'phi_fallback',10,       @isnumeric);
addParameter(p, 'strict_phi_source', true, @islogical);
addParameter(p, 'ka_bounds',   [0, 3], @isnumeric);
addParameter(p, 'max_snow_depth', 5,  @isnumeric);
addParameter(p, 'output_all_H', true,    @islogical);
parse(p, file, varargin{:});

Ta_fallback  = p.Results.Ta;   % 粗略回退值（当 ATL09 匹配不到时使用）
Ta           = Ta_fallback;
CF           = p.Results.CF;
beam         = p.Results.beam;
verbose      = p.Results.verbose;
phi_fallback = p.Results.phi_fallback;
strict_phi_source = p.Results.strict_phi_source;
ka_bounds    = p.Results.ka_bounds;
output_all_H = p.Results.output_all_H;
max_snow_depth = p.Results.max_snow_depth;
atl09_file         = p.Results.atl09_file;
atl09_Ta_dataset   = p.Results.atl09_Ta_dataset;
atl09_gt           = p.Results.atl09_gt;
atl09_radius_m     = p.Results.atl09_radius_m;
Ta_is_twoway       = p.Results.Ta_is_twoway;
atl09_mat_lat_var = p.Results.atl09_mat_lat_var;
atl09_mat_lon_var = p.Results.atl09_mat_lon_var;
atl09_mat_Ta_var  = p.Results.atl09_mat_Ta_var;

%% ===== 加载退卷积结果 =====
% 来源：f_deconvolution.m
% 包含：Sc（退卷积廓线），Depth_o，data_So，noise
load(file);
lat    = data_So(1, :);
lon    = data_So(2, :);
[~, n] = size(Sc);

%% ===== 从ATL09获取Ta值 =====
if ~isempty(atl09_file)
    fprintf('[步骤] 从ATL09读取Ta并匹配经纬度...\n');
    try
        if iscell(atl09_file)
            fprintf('[步骤] 尝试从%d个ATL09文件读取Ta值...\n', length(atl09_file));
        else
            fprintf('[步骤] 尝试从ATL09文件 %s 读取Ta值...\n', atl09_file);
        end
        fprintf('[步骤] 使用的变量名：lat=%s, lon=%s, Ta=%s\n', atl09_mat_lat_var, atl09_mat_lon_var, atl09_mat_Ta_var);
        Ta = f_Ta_from_atl09_v2( ...
            atl09_file, lat, lon, ...
            'Ta_dataset', atl09_Ta_dataset, ...
            'gt',          atl09_gt, ...
            'radius_m',   atl09_radius_m, ...
            'Ta_is_twoway', Ta_is_twoway, ...
            'mat_lat_var', atl09_mat_lat_var, ...
            'mat_lon_var', atl09_mat_lon_var, ...
            'mat_Ta_var',  atl09_mat_Ta_var, ...
            'verbose',     true); % 启用详细输出
        % 匹配结果质量控制：无效点回退
        % （不把 Ta==1 判为无效：后续仍会把 Ta clamp 到 0.99）
        valid_Ta = isfinite(Ta) & Ta > 0;
        fprintf('[步骤] 从ATL09获取的Ta值统计：有效点数=%d/%d\n', sum(valid_Ta), length(Ta));
        if any(valid_Ta)
            % 只对无效点回退，保留有效点的实际值
            fprintf('[步骤] 有效Ta值范围：[%.4f, %.4f]\n', min(Ta(valid_Ta)), max(Ta(valid_Ta)));
            Ta(~valid_Ta) = Ta_fallback;
        else
            % 如果所有点都无效，才全部回退
            fprintf('[步骤] 未找到有效Ta值，全部回退为 %.3f\n', Ta_fallback);
            Ta = Ta_fallback * ones(size(lat));
        end
        Ta = min(max(Ta, 1e-4), 0.99);
    catch ME
        warning('[步骤] ATL09 Ta 提取失败：%s；回退为 Ta=%g', ME.message, Ta_fallback);
        Ta = Ta_fallback * ones(size(lat));
    end
else
    if verbose
        fprintf('[步骤] 使用输入Ta回退值 %.3f（未提供atl09_file）...\n', Ta_fallback);
    end
end

%% ===== 读取 Np_vec 和 phi_vec =====
% 来源：f_profile_processing.m → *_profile.mat
% Np_vec：每廓线雪面光子数 [counts/shot]
% phi_vec：每廓线沿轨坡度 [度]
a_dir = dir(file);
profile_file = fullfile(a_dir.folder, strrep(a_dir.name, 'deconv_', ''));

Np_vec = [];
phi_vec = [];      % 原始坡度（来自profile）
phi = [];          % 实际用于反演的坡度

if exist(profile_file, 'file')
    pdata = load(profile_file);
    if isfield(pdata, 'Np_vec') && isfield(pdata, 'phi_vec')
        Np_vec  = pdata.Np_vec;
        phi_vec = pdata.phi_vec;
        phi     = smoothdata(phi_vec, 'movmedian', 5);
        if verbose
            fprintf('[输入] Np_vec [%.4f, %.4f] counts/shot\n', min(Np_vec), max(Np_vec));
            fprintf('[输入] phi_vec(raw) [%.2f, %.2f] deg; phi(used) [%.2f, %.2f] deg\n', ...
                    min(phi_vec), max(phi_vec), min(phi), max(phi));
        end
    else
        warning('profile文件缺少Np_vec/phi_vec，请用新版f_profile_processing重新生成');
    end
else
    warning('未找到profile文件: %s', profile_file);
end

% 长度对齐
if ~isempty(Np_vec)
    nL = length(Np_vec);
    if nL ~= n
        nu = min(n, nL);
        Np_vec  = Np_vec(1:nu);
        phi_vec = phi_vec(1:nu);
        phi     = phi(1:nu);
        Sc      = Sc(:, 1:nu);
        lat     = lat(1:nu);
        lon     = lon(1:nu);
        n       = nu;
        warning('Np_vec长度(%d)与Sc列数(%d)不符，截断至%d', nL, n, nu);
    end
else
    % 回退：从Sc峰值估算Np，使用备用坡度（可按strict_phi_source禁止）
    if strict_phi_source
        error(['未找到有效phi_vec。当前已启用strict_phi_source=true，', ...
               '请先运行f_profile_processing生成phi_vec，或将strict_phi_source设为false。']);
    end
    ind_s = find(Depth_o >= -0.3 & Depth_o <= 0.3);
    if isempty(ind_s); [~, ind_s] = min(abs(Depth_o)); end
    Np_vec = zeros(1, n);
    for ii = 1:n
        v = Sc(ind_s, ii); v(isnan(v)) = 0;
        Np_vec(ii) = max(v);
    end
    phi_vec = phi_fallback * ones(1, n);
    phi     = phi_vec;
    if verbose
        fprintf('[回退] Sc峰值估算Np，坡度备用值%.1f度\n', phi_fallback);
    end
end

%% ====================================================================
%% 步骤 f：真实光子数修正
%%
%% 物理依据：
%%   ICESat-2 PMT 探测器存在死时间效应，导致高信号时漏计光子。
%%   论文公式(17)（经验二次回归模型）：
%%     Ns = (1.06 + 0.00714·Np)·Np + (-0.07 + 0.00152·φ)·φ + 0.58
%%   当 Ns > 3 counts 时需校正，否则 Ns = Np。
%%
%% 输入：Np_vec（实测PNPS），phi_vec（坡度）
%% 输出：Ns_vec（真实期望光子数），仅用于步骤g，不修改Sc廓线
%% ====================================================================
if verbose; fprintf('\n[步骤f] PMT死时间修正...\n'); end

% 确认Np单位：打印Np_vec的范围
fprintf('Np_vec range: %.3f - %.3f (median %.3f)\n', min(Np_vec), max(Np_vec), median(Np_vec));

[Ns_vec, f_info] = f_photon_correction( ...
    Np_vec, phi, 'beam', beam, 'verbose', verbose);

if verbose
    vm = Np_vec > 0;
    if any(vm)
        r_mean = mean(Ns_vec(vm) ./ Np_vec(vm));
    else
        r_mean = NaN;
    end
    fprintf('[步骤f] Ns范围[%.4f, %.4f], 平均校正比=%.4f\n', ...
            min(Ns_vec), max(Ns_vec), r_mean);
end

%% ====================================================================
%% 步骤 g：雪面反照率计算
%%
%% 物理依据：
%%   论文公式(8)（激光雷达方程反演）：
%%     ρ = Ns · Cs · (1/Ta²) · cos(φ)
%%   其中：
%%     Cs = Cc · CF = 0.12 × 0.52 = 0.0624   ← 系统标定常数
%%     Ta = 0.81（单程大气透过率，学长给出）
%%     φ  = phi_vec（步骤b计算的沿轨坡度）
%%
%%   同时计算 ASR = ρ·Ta² = Ns·Cs·cos(φ)（不依赖大气，更稳健）
%%
%% 输入：Ns_vec（步骤f输出），phi_vec（坡度）
%% 输出：rho_vec（真实反照率），ASR_vec（表观反照率）
%% ====================================================================
% Ta：已在前面步骤中从ATL09获取，这里直接使用
if verbose
    fprintf('\n[步骤g] 雪面反照率计算...\n');
    fprintf('[步骤g] 使用的Ta值范围：[%.4f, %.4f]\n', min(Ta), max(Ta));
end

[rho_vec, ASR_vec, ~, g_info] = f_snow_albedo( ...
    Ns_vec, phi, 'Ta', Ta, 'CF', CF, 'verbose', verbose);

%% ====================================================================
%% 步骤 h+i：严格解析雪深反演（a -> k_a -> p(L) -> H）
%%
%% 步骤h：由反照率积分式反解吸收系数 k_a
%%   a = ∫β(z)dz / ∫β(z)e^(2k_a z)dz
%%   对给定 β(z), a，使用单变量求根获得 k_a（不再做外层固定点迭代）。
%%
%% 步骤i：由 k_a 构建无吸收路径分布并直接矩法求雪深
%%   p_raw(L) = β(z)e^(k_a L), L=2z
%%   直接对观测路径分布归一化并计算矩，不再引入Gamma拟合。
%%   雪深输出：
%%     H_mean = <L>/2
%%     H_l2   = (<L²>/k_sd)^(1/3)
%%     H_l3   = (<L³>/k_sd²)^(1/5)
%% ====================================================================
ind_Depth = find(Depth_o >= 0 & Depth_o < 1);
Depth_ea  = Depth_o(ind_Depth);
if isrow(Depth_ea); Depth_ea = Depth_ea'; end

snow  = nan(1, n);
snowA = nan(1, n);
ka_vec     = nan(1, n);
H_mean_vec = nan(1, n);
H_l2_vec   = nan(1, n);
H_l3_vec   = nan(1, n);
residual_F = nan(1, n);

n_used_rho    = 0;
n_fallback    = 0;
n_solved      = 0;
n_fail        = 0;
fail_no_beta  = 0;
fail_no_brkt  = 0;
fail_root     = 0;
fail_pdist    = 0;
fail_gamma    = 0;
fail_ksd      = 0;

for i = 1:n
    % 输入廓线：退卷积后的散射廓线（不做幅值修改）
    b_z = Sc(ind_Depth, i);
    if isrow(b_z); b_z = b_z'; end
    b_z(isnan(b_z)) = 0;
    b_z(b_z < 0)    = 0;
    if sum(b_z) == 0; continue; end

    % 判断步骤g的反照率是否有效
    rho_i      = rho_vec(i);
    use_rho    = ~isnan(rho_i) && rho_i > 0 && rho_i < 1;
    z = Depth_ea;
    L = 2 * z;  % 光程 L = 2z（往返）

    % 步骤h：确定反照率 a（优先使用步骤g）
    if use_rho
        a = rho_i;
        n_used_rho = n_used_rho + 1;
    else
        % 若rho被物理约束裁剪为NaN，回退使用rho_raw并做区间裁剪
        if isfield(g_info, 'rho_raw')
            rho_raw_i = g_info.rho_raw(i);
            if isfinite(rho_raw_i) && rho_raw_i > 0
                a = rho_raw_i;
                n_fallback = n_fallback + 1;
            else
                n_fallback = n_fallback + 1;
                fail_root = fail_root + 1;
                continue;
            end
        else
            n_fallback = n_fallback + 1;
            fail_root = fail_root + 1;
            continue;
        end
    end
    a = min(max(a, 1e-4), 0.99);

    integral_beta = trapz(z, b_z);
    if integral_beta <= 0
        n_fail = n_fail + 1;
        fail_no_beta = fail_no_beta + 1;
        continue;
    end

    % 反解 k_a：F(k_a)=∫βe^(2k_a z)dz - (∫βdz)/a = 0
    target_int = integral_beta / a;
    F = @(ka) trapz(z, b_z .* exp(2 * ka * z)) - target_int;

    kmin = ka_bounds(1);
    kmax = ka_bounds(2);
    k_a = NaN;

    % 先尝试在默认区间找根号变化
    Fmin = F(kmin);
    Fmax = F(kmax);
    has_bracket = isfinite(Fmin) && isfinite(Fmax) && (Fmin * Fmax <= 0);

    % 若无括号，执行自适应扩展上限并扫描找符号变化区间
    if ~has_bracket
        kmax_try_list = [ka_bounds(2), 6, 10];
        F_grid = [];
        k_grid = [];
        for kk_try = 1:numel(kmax_try_list)
            kmax_try = kmax_try_list(kk_try);
            k_grid = linspace(kmin, kmax_try, 120);
            F_grid = nan(size(k_grid));
            for kk = 1:numel(k_grid)
                F_grid(kk) = F(k_grid(kk));
            end
            for kk = 1:numel(k_grid)-1
                if isfinite(F_grid(kk)) && isfinite(F_grid(kk+1)) && F_grid(kk) * F_grid(kk+1) <= 0
                    kmin = k_grid(kk);
                    kmax = k_grid(kk+1);
                    has_bracket = true;
                    break;
                end
            end
            if has_bracket
                break;
            end
        end

        % 仍无括号：用最小残差近似解（允许较小数值误差）
        if ~has_bracket
            [min_absF, idx_min] = min(abs(F_grid));
            if isfinite(min_absF)
                k_try = k_grid(idx_min);
                rel_res = min_absF / max(target_int, 1e-8);
                if rel_res < 1e-2
                    k_a = k_try;
                else
                    n_fail = n_fail + 1;
                    fail_no_brkt = fail_no_brkt + 1;
                    continue;
                end
            else
                n_fail = n_fail + 1;
                fail_no_brkt = fail_no_brkt + 1;
                continue;
            end
        end
    end

    if isnan(k_a)
        k_a = fzero(F, [kmin, kmax]);
    end
    % fzero仍失败时，用fminbnd最小化|F|再做残差判定
    if ~isfinite(k_a) || k_a <= 0
        obj = @(ka) abs(F(ka));
        k_try = fminbnd(obj, ka_bounds(1), 10);
        rel_res = abs(F(k_try)) / max(target_int, 1e-8);
        if isfinite(k_try) && k_try > 0 && rel_res < 1e-2
            k_a = k_try;
        end
    end
    if ~isfinite(k_a) || k_a <= 0
        n_fail = n_fail + 1;
        fail_root = fail_root + 1;
        continue;
    end

    % 步骤i：构建 p(L) 并直接计算观测矩
    p_raw = b_z .* exp(k_a * L);
    p_int = trapz(L, p_raw);
    if p_int <= 0
        n_fail = n_fail + 1;
        fail_pdist = fail_pdist + 1;
        continue;
    end
    p_obs = p_raw / p_int;

    L_mean = trapz(L, L .* p_obs);
    L2_mean = trapz(L, (L.^2) .* p_obs);
    L3_mean = trapz(L, (L.^3) .* p_obs);
    if ~(isfinite(L_mean) && isfinite(L2_mean) && isfinite(L3_mean) && L_mean > 0 && L2_mean > 0)
        n_fail = n_fail + 1;
        fail_pdist = fail_pdist + 1;
        continue;
    end

    H_mean = L_mean / 2;
    if H_mean <= 0 || H_mean > max_snow_depth
        n_fail = n_fail + 1;
        fail_pdist = fail_pdist + 1;
        continue;
    end

    % Warren 关系估算 k_sd，供 H_l2/H_l3 使用（若失败不影响H_mean主解）
    term1 = (1 - a) / 8.43;
    R = (1 / k_a) * (term1^2);
    H_l2 = NaN;
    H_l3 = NaN;
    if R > 0 && isfinite(R)
        k_d = 0.65 * sqrt(k_a / R);
        k_sd = (k_d^2) / (3 * k_a) - k_a;
        if isfinite(k_sd) && k_sd > 0 && H_mean > 0
            H_l2_try = (L2_mean / k_sd)^(1/3);
            H_l3_try = (L3_mean / (k_sd^2))^(1/5);
            if isfinite(H_l2_try) && isfinite(H_l3_try) && H_l2_try > 0 && H_l3_try > 0
                H_l2 = H_l2_try;
                H_l3 = H_l3_try;
            else
                fail_ksd = fail_ksd + 1;
            end
        else
            fail_ksd = fail_ksd + 1;
        end
    else
        fail_ksd = fail_ksd + 1;
    end

    n_solved      = n_solved + 1;
    snow(i)       = H_mean;
    snowA(i)      = a;
    ka_vec(i)     = k_a;
    H_mean_vec(i) = H_mean;
    H_l2_vec(i)   = H_l2;
    H_l3_vec(i)   = H_l3;
    residual_F(i) = F(k_a);
end

if verbose
    fprintf('\n[步骤h+i] 反演完成: 成功%d/%d, 失败%d\n', n_solved, n, n_fail);
    if n_fail > 0
        fprintf('  失败原因(no_beta/no_brkt/root/pdist/gamma/ksd): %d/%d/%d/%d/%d/%d\n', ...
            fail_no_beta, fail_no_brkt, fail_root, fail_pdist, fail_gamma, fail_ksd);
    end
    if n_solved > 0
        v = ~isnan(snow);
        fprintf('  雪深[m]范围[%.3f, %.3f], 均值%.3f\n', ...
                min(snow(v)), max(snow(v)), mean(snow(v)));
        if output_all_H
            v2 = ~isnan(H_l2_vec);
            v3 = ~isnan(H_l3_vec);
            if any(v2); fprintf('  H_l2均值=%.3f m\n', mean(H_l2_vec(v2))); end
            if any(v3); fprintf('  H_l3均值=%.3f m\n', mean(H_l3_vec(v3))); end
        end
    end
end

% 打包输出信息
info_out.f_info     = f_info;
info_out.g_info     = g_info;
info_out.rho_vec    = rho_vec;
info_out.ASR_vec    = ASR_vec;
info_out.Ns_vec     = Ns_vec;
info_out.phi        = phi;       % used in retrieval
info_out.phi_vec    = phi_vec;   % raw from profile
info_out.n_used_rho = n_used_rho;
info_out.n_fallback = n_fallback;
info_out.ka_vec     = ka_vec;
info_out.H_mean_vec = H_mean_vec;
info_out.residual_F = residual_F;
info_out.n_fail     = n_fail;
info_out.fail_no_beta = fail_no_beta;
info_out.fail_no_brkt = fail_no_brkt;
info_out.fail_root    = fail_root;
info_out.fail_pdist   = fail_pdist;
info_out.fail_gamma   = fail_gamma;
info_out.fail_ksd     = fail_ksd;
if output_all_H
    info_out.H_l2_vec = H_l2_vec;
    info_out.H_l3_vec = H_l3_vec;
end

end
