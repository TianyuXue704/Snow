function [rho, ASR, Ta2, info] = f_snow_albedo(Ns_vec, phi_vec, varargin)

%% ========== 参数解析 ==========
p = inputParser;
addRequired(p, 'Ns_vec',  @isnumeric);
addRequired(p, 'phi_vec', @isnumeric);
addParameter(p, 'Ta',      0.9, @isnumeric);   % 单程大气透过率
addParameter(p, 'CF',      0.52, @isnumeric);   % 在轨标定因子（论文IV-B2）
addParameter(p, 'Cc',      0.12, @isnumeric);   % 硬件参数（论文Table I）
addParameter(p, 'verbose', true, @islogical);
parse(p, Ns_vec, phi_vec, varargin{:});

Ta      = p.Results.Ta;
CF      = p.Results.CF;
Cc      = p.Results.Cc;
verbose = p.Results.verbose;

%% ========== 系统标定常数 ==========
% 论文公式(7)：Cs = Cc * CF
Cs  = Cc * CF;          % = 0.12 * 0.52 = 0.0624
Ta2 = Ta.^2;            % 双程大气透过率（允许Ta为向量）

%% ========== 反照率计算 ==========
% 论文公式(8)：rho = Ns * Cs * (1/Ta^2) * cos(phi)
cos_phi = cosd(phi_vec);            % cos(phi)，phi单位为度

% 允许 Ta 为标量或向量：必须用逐元素除法
rho = Ns_vec .* Cs .* (1./Ta2) .* cos_phi;

% 表观反照率：ASR = rho * Ta^2 = Ns * Cs * cos(phi)
% 注意：ASR与Ta无关，更稳健
ASR = Ns_vec .* Cs .* cos_phi;

%% ========== 物理约束 ==========
% 反照率物理上在[0,1]，超出范围标记为NaN
rho_raw = rho;
rho(rho > 1.0) = NaN;
rho(rho < 0)   = NaN;
ASR(ASR > 1.0) = NaN;
ASR(ASR < 0)   = NaN;

%% ========== 输出信息 ==========
info.Cs       = Cs;
info.Ta       = Ta;
info.Ta2      = Ta2;
info.CF       = CF;
info.Cc       = Cc;
info.rho_raw  = rho_raw;           % 未裁剪的原始反照率（用于诊断）
info.n_invalid_rho = sum(isnan(rho));
info.n_invalid_ASR = sum(isnan(ASR));

if verbose
    n = length(Ns_vec);
    valid_rho = ~isnan(rho);
    if isscalar(Ta)
        fprintf('[f_snow_albedo] Cs=%.6f, Ta=%.4f(Ta^2=%.4f), 有效rho=%d/%d\n', ...
                Cs, Ta, Ta2, sum(valid_rho), n);
    else
        Ta_ok = Ta(isfinite(Ta) & Ta > 0);
        Ta_med = NaN; Ta_min = NaN; Ta_max = NaN;
        if ~isempty(Ta_ok)
            Ta_med = median(Ta_ok);
            Ta_min = min(Ta_ok);
            Ta_max = max(Ta_ok);
        end
        fprintf('[f_snow_albedo] Cs=%.6f, Ta(向量) med=%.4f [%.4f, %.4f], 有效rho=%d/%d\n', ...
                Cs, Ta_med, Ta_min, Ta_max, sum(valid_rho), n);
    end
    if sum(valid_rho) > 0
        fprintf('  rho范围[%.4f, %.4f], 均值=%.4f\n', ...
                min(rho(valid_rho)), max(rho(valid_rho)), mean(rho(valid_rho)));
        fprintf('  ASR范围[%.4f, %.4f], 均值=%.4f\n', ...
                min(ASR(~isnan(ASR))), max(ASR(~isnan(ASR))), mean(ASR(~isnan(ASR))));
    end
    if info.n_invalid_rho > 0
        fprintf('  警告: %d 条rho超出[0,1]已置NaN\n', ...
                info.n_invalid_rho);
    end
end

end
