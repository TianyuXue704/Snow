function [Ns, correction_info] = f_photon_correction(Np, phi_deg, varargin)

p = inputParser;
addRequired(p,  'Np',         @isnumeric);
addRequired(p,  'phi_deg',    @isnumeric);
addParameter(p, 'noise_rate', 2.0,   @isnumeric);   % MHz
addParameter(p, 'beam',       'strong', @ischar);
addParameter(p, 'threshold',  true,  @islogical);
addParameter(p, 'plot',       false, @islogical);
addParameter(p, 'verbose',    true,  @islogical);
parse(p, Np, phi_deg, varargin{:});

fn       = p.Results.noise_rate;   % 噪声率 [MHz]
beam     = p.Results.beam;
do_thr   = p.Results.threshold;
do_plot  = p.Results.plot;
verbose  = p.Results.verbose;

%% ========== 探测器数量 ==========
% 强光束: 16个PMT; 弱光束: 4个PMT(每4路合并输出)
if strcmpi(beam, 'strong')
    M = 16;
else
    M = 4;
end

%% ========== 输入验证与广播 ==========
orig_size = size(Np);
Np = Np(:);          % 展平为列向量
phi_deg = phi_deg(:);

% 坡度广播: 若phi为标量则扩展
if isscalar(phi_deg)
    phi_deg = phi_deg * ones(size(Np));
end

if length(phi_deg) ~= length(Np)
    error('f_photon_correction: Np与phi_deg长度不匹配，请检查输入');
end

% 负值保护
Np(Np < 0) = 0;
phi_deg(phi_deg < 0) = 0;

N = length(Np);

%% ========== 核心校正: 论文公式(17) ==========
% Ns = (mu_p1 + mu_p2 * Np) * Np + (mu_s1 + mu_s2 * phi) * phi + mu_c
%
% 拟合参数(Zhou et al. 2023, 基于ICESat-2/ATLAS系统参数):
mu_p1 =  1.06;
mu_p2 =  7.14e-3;
mu_s1 = -0.07;
mu_s2 =  1.52e-3;   % 注意: 论文原式坡度单位为度(degree)
mu_c  =  0.58;

phi = phi_deg;       % 保持度为单位，与论文一致

Ns_raw = (mu_p1 + mu_p2 .* Np) .* Np + ...
         (mu_s1 + mu_s2 .* phi) .* phi + ...
         mu_c;

%% ========== 阈值判断 ==========
% 论文指出: 当Ns/M > 0.2时死时间效应不可忽略
% 对16通道强光束: Ns > 0.2 * 16 = 3.2 counts
threshold_Ns = 0.2 * M;   % = 3.2 for strong beam

if do_thr
    % 仅对超过阈值的点使用校正值，否则保留原始Np
    need_correction = (Np >= threshold_Ns);
    Ns = Np;                              % 默认: 不校正
    Ns(need_correction) = Ns_raw(need_correction);
else
    % 无条件应用校正
    need_correction = true(N, 1);
    Ns = Ns_raw;
end

% 物理约束: Ns必须 >= Np (死时间只会导致低估，不会高估)
Ns = max(Ns, Np);
% Ns必须为非负
Ns = max(Ns, 0);

%% ========== 构建输出信息结构体 ==========
correction_ratio        = zeros(N, 1);
valid_mask              = Np > 0;
correction_ratio(valid_mask) = Ns(valid_mask) ./ Np(valid_mask);
correction_ratio(~valid_mask) = 1.0;

delta_Ns = Ns - Np;

correction_info.Np_input         = reshape(Np, orig_size);
correction_info.phi_input         = reshape(phi_deg, orig_size);
correction_info.Ns_raw            = reshape(Ns_raw, orig_size);
correction_info.need_correction   = reshape(need_correction, orig_size);
correction_info.correction_ratio  = reshape(correction_ratio, orig_size);
correction_info.delta_Ns          = reshape(delta_Ns, orig_size);
correction_info.threshold_used    = threshold_Ns;
correction_info.M_detectors       = M;

% 统计摘要
n_corrected = sum(need_correction);
correction_info.stats.n_total         = N;
correction_info.stats.n_corrected     = n_corrected;
correction_info.stats.pct_corrected   = 100 * n_corrected / N;
correction_info.stats.mean_ratio      = mean(correction_ratio);
correction_info.stats.max_delta_Ns    = max(delta_Ns);
correction_info.stats.mean_delta_Ns   = mean(delta_Ns(need_correction));

%% ========== 打印统计 ==========
if verbose
    fprintf('[f_photon_correction] %s beam, 阈值=%.1f, 校正点=%d/%d(%.1f%%)\n', ...
            beam, threshold_Ns, n_corrected, N, correction_info.stats.pct_corrected);
    if n_corrected > 0
        fprintf('  校正比Ns/Np均值=%.4f, 最大增量=+%.3f, 平均增量=+%.3f counts/shot\n', ...
                mean(correction_ratio(need_correction)), ...
                correction_info.stats.max_delta_Ns, ...
                correction_info.stats.mean_delta_Ns);
    end
end

%% ========== 恢复原始形状 ==========
Ns = reshape(Ns, orig_size);

%% ========== 可选绘图 ==========
if do_plot
    figure('Color','w','Name','光子数死时间修正诊断','Position',[100 100 1100 400]);
    
    % --- 子图1: Np vs Ns散点图 ---
    subplot(1,3,1);
    Np_plot = correction_info.Np_input(:);
    Ns_plot = Ns(:);
    nc_mask = correction_info.need_correction(:);
    
    scatter(Np_plot(~nc_mask), Ns_plot(~nc_mask), 15, [0.6 0.6 0.6], ...
            'filled', 'DisplayName', sprintf('未校正 (Np < %.1f)', threshold_Ns));
    hold on;
    scatter(Np_plot(nc_mask),  Ns_plot(nc_mask),  15, [0.8500 0.3250 0.0980], ...
            'filled', 'DisplayName', '已校正');
    
    ref_x = linspace(0, max(Np_plot)*1.05, 50);
    plot(ref_x, ref_x, 'k--', 'LineWidth', 1, 'DisplayName', 'Ns = Np (无校正)');
    
    xlabel('实测 N_p (counts/shot)', 'FontSize', 11);
    ylabel('修正后 N_s (counts/shot)', 'FontSize', 11);
    title('校正关系', 'FontSize', 12);
    legend('Location','northwest','FontSize',9);
    grid on; box on;
    
    % --- 子图2: 校正量 ΔNs 随 Np 变化 ---
    subplot(1,3,2);
    [Np_sorted, sort_idx] = sort(Np_plot);
    delta_sorted = delta_Ns(sort_idx);
    phi_sorted   = phi_deg(sort_idx);
    
    scatter(Np_sorted, delta_sorted, 10, phi_sorted, 'filled');
    cb = colorbar;
    cb.Label.String = '坡度 φ (°)';
    colormap(gca, parula);
    xlabel('实测 N_p (counts/shot)', 'FontSize', 11);
    ylabel('校正量 ΔN_s = N_s - N_p', 'FontSize', 11);
    title('校正量 vs PNPS (颜色=坡度)', 'FontSize', 12);
    grid on; box on;
    xline(threshold_Ns, 'r--', 'LineWidth', 1.2, ...
          'Label', sprintf('阈值 %.1f', threshold_Ns), 'LabelVerticalAlignment', 'bottom');
    
    % --- 子图3: 校正比值分布直方图 ---
    subplot(1,3,3);
    ratio_corrected = correction_ratio(nc_mask);
    if ~isempty(ratio_corrected)
        histogram(ratio_corrected, 30, 'FaceColor', [0 0.4470 0.7410], ...
                  'EdgeColor', 'none', 'Normalization', 'probability');
        xline(1.0, 'r--', 'LineWidth', 1.5, 'Label', '无修正');
        xlabel('校正比值 N_s / N_p', 'FontSize', 11);
        ylabel('概率', 'FontSize', 11);
        title(sprintf('校正比值分布\n(均值=%.3f)', mean(ratio_corrected)), 'FontSize', 12);
        grid on; box on;
    else
        text(0.5, 0.5, '无需校正的点', 'HorizontalAlignment', 'center', ...
             'Units', 'normalized', 'FontSize', 12);
    end
    
    sgtitle('ICESat-2 PMT死时间修正诊断 (Zhou et al. 2023)', ...
            'FontSize', 13, 'FontWeight', 'bold');
end

end


%% =====================================================================
%% 辅助函数: 从ATL08/ATL03产品计算净PNPS
%% =====================================================================
function Np_net = compute_net_pnps(pnps_raw, bckgrd_rate, h_canopy)
% COMPUTE_NET_PNPS  从ATL08原始PNPS减去背景噪声贡献
%
% 对应论文Step(2): Np_net = Np08 - hs * NA
%
% 输入:
%   pnps_raw    - ATL08中的原始PNPS [counts/shot]
%   bckgrd_rate - ATL03的背景计数率 NA [counts/shot/m]
%   h_canopy    - ATL08的信号积分高度 hs [m]
%
% 输出:
%   Np_net      - 噪声扣除后的净PNPS
%
% 说明:
%   噪声光子在高度范围内均匀分布，每米的期望噪声计数 = bckgrd_rate
%   在积分高度hs内的总噪声期望值 = bckgrd_rate * h_canopy

Np_net = pnps_raw - bckgrd_rate .* h_canopy;
Np_net = max(Np_net, 0);  % 不允许负值

end


%% =====================================================================
%% 辅助函数: 计算总坡度
%% =====================================================================
function phi_total = compute_total_slope(s_par_deg, s_perp_deg)
% COMPUTE_TOTAL_SLOPE  由沿轨和跨轨坡度合成总坡度
%
% 对应论文公式(5):
%   phi ≈ st = arctan(sqrt(tan²(s_par) + tan²(s_perp)))
%
% 输入:
%   s_par_deg  - 沿轨坡度 [度]，从ATL08信号光子线性拟合获得
%   s_perp_deg - 跨轨坡度 [度]，从ATL06产品获得并平均到对应段
%
% 输出:
%   phi_total  - 总坡度 [度]

phi_total = atand(sqrt(tand(s_par_deg).^2 + tand(s_perp_deg).^2));

end
