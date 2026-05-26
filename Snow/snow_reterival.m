function [snow, lat, snowA] = snow_reterival(file)
    % 加载数据
    load(file); % 假设包含 data_So, Depth_o, Sc
    
    % --- 数据预处理 ---
    lat = data_So(1:2, :);
    
    % 确定深度索引 (确保 z 为正值且在合理范围内)
    ind_Depth = find(Depth_o >= 0 & Depth_o < 20); 
    
    % 强制转换为列向量，方便后续点乘和积分
    Depth_ea = Depth_o(ind_Depth); 
    if isrow(Depth_ea); Depth_ea = Depth_ea'; end
    
    alpha_init = 0.07; % 初始吸收系数猜测值 k_a (Unit: m^-1)
    
    [~, n] = size(Sc); 
    snow = nan(1, n);  % 初始化为 NaN，区分无效值
    snowA = nan(1, n); 
    
    % 遍历每一条廓线
    for i = 1:n
        % 提取第 i 条廓线数据
        b_z = Sc(ind_Depth, i); 
        
        % 强制转换为列向量
        if isrow(b_z); b_z = b_z'; end
        
        % 数据清洗：去除 NaN，如有负值需归零
        b_z(isnan(b_z)) = 0;
        b_z(b_z < 0) = 0;
        
        % 如果整条廓线为空或信号太弱，跳过
        if sum(b_z) == 0
            continue; 
        end
        
        % --- 迭代初始化 ---
        k_a = alpha_init;
        tolerance = 0.001;
        max_iter = 100;
        H_final = 0;
        a_final = 0;
        converged = false;
        
        for iter = 1:max_iter
            % --- 步骤 3: 估算无吸收光程分布 p(L) 和反照率 a ---
            % 公式: p(L) = beta(z) * exp(k_a * L)
            % 关键修正: 光程 L = 2 * z
            L = 2 * Depth_ea;
            z = Depth_ea;
            p_L = b_z .* exp(k_a * L); 
            
            % 积分计算 (利用梯形法则)
            % 注意: 积分变量如果是 dz, 上下都乘 dz 会抵消，不影响比值
            integral_beta = trapz(L, b_z);
            integral_p    = trapz(L, p_L);
            
            if integral_p == 0; break; end 
            
            % 计算反照率 a
            a = integral_beta / integral_p;
            % 物理限制: 反照率不能超过 1
            if a > 0.99; a = 0.99; end 
            
            % --- 步骤 4: 计算雪深 H ---
            % 理论公式: H = <L> / 2
            % <L> 是光程的一阶矩 (Mean Path Length)
            % <L> = ∫(L * p(L)) dL / ∫ p(L) dL
            % 换元: ∫(2z * p(2z)) * 2dz / ∫ p(2z) * 2dz
            % 分子分母中的 2dz 抵消，只剩下 L=2z
            
            integral_Lp = trapz(z, z .* p_L); % 分子
            L_mean = 2*integral_Lp / integral_p;       % <L>
            
            H_est = L_mean / 2; % H = <L>/2
            
            if H_est < 1e-3; break; end % 防止 H 过小导致除以零
            
            % --- 步骤 5: 计算扩散散射系数 k_sd ---
            % 理论公式: k_sd = <L^2> / (H^3)
            % <L^2> = ∫(L^2 * p(L)) dL / ∫ p(L) dL
            
            integral_L2p = trapz(z, (z.^2) .* p_L);
            L2_mean = 4*integral_L2p / integral_p;     % <L^2>
            
            % 注意: 此处 H_est = <L>/2, 所以公式等价于 k_sd = <L^2> / (<L>/2)^3
            k_sd = L2_mean / ((H_est)^3);
            
            % --- 步骤 6: 计算雪粒径 R 和 扩散衰减系数 k_d ---
            % R 公式 (Warren 1982)
            term1 = (1 - a) / 8.43;
            R = (1 / k_a) * (term1^2);
            
            % k_d 公式
            if R == 0; break; end
            k_d = 0.65 * sqrt(k_a / R);
            
            % --- 步骤 7: 更新吸收系数 k_a' ---
            % 公式: k_a' = k_d^2 / (3 * (k_sd + k_a))
            denominator = 3 * (k_sd + k_a);
            if denominator == 0; break; end
            
            k_a_new = (k_d^2) / denominator;
            
            % --- 步骤 8: 收敛判断 ---
            if abs(k_a_new - k_a) < tolerance
                H_final = H_est;
                a_final = a;
                converged = true;
                break;
            end
            
            % 更新 k_a 进行下一次迭代
            k_a = k_a_new;
            
            % 安全阀：防止 k_a 发散成负数或无穷大
            
        end
        
        % 仅存储收敛的结果
        if converged
            snow(i) = H_final;
            snowA(i) = a_final;
        else
            % 如果未收敛，根据需求可以存 NaN 或最后一次计算值
            snow(i) = NaN; 
            snowA(i) = NaN;
        end
    end
end