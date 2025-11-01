function n_chi = calculate_n_kld_ign(k, epsilon, delta)
    % 根据 KLD-Sampling 算法计算目标粒子数 n_chi。
    %
    % 该函数实现了 Dieter Fox 论文《Adapting the Sample Size in Particle 
    % Filters Through KLD-Sampling》中的核心公式。
    %
    % Args:
    %   k (double): 当前观测到的非空“箱子”(bin)的数量。
    %   epsilon (double): 允许的最大近似误差 (e.g., 0.05)。
    %   delta (double): 置信度参数, (1-delta) 是概率 (e.g., 0.01 对应 99% 置信度)。
    %
    % Returns:
    %   n_chi (double): 为满足 KLD 标准所需的目标粒子总数 (Inf 或一个整数)。

    % 论文中明确指出，该公式仅在 k > 1 时有意义
    if k <= 1
        % 在 k=1 时，我们无法计算卡方值。
        % 返回一个极大值（Inf），以确保循环继续，直到 k > 1。
        n_chi = Inf;
        return;
    end

    % 1. 计算标准正态分布的 z_{1-delta} 分位数
    % (norminv 是 Statistics and Machine Learning Toolbox 的一部分)
    z_1_minus_delta = norminv(1 - delta, 0, 1);
    
    % 2. 这是 Wilson-Hilferty 变换，近似计算 chi^2 分位数
    k_minus_1 = double(k - 1);
    term_1 = 2.0 / (9.0 * k_minus_1);
    term_2 = sqrt(term_1) * z_1_minus_delta;
    
    % 括号内的 { ... }^3
    inner_term_cubed = (1.0 - term_1 + term_2)^3;
    
    % 近似的卡方值 chi_squared
    chi_squared = k_minus_1 * inner_term_cubed;
    
    % 3. 这是最终的 n_chi 计算公式 (13) [cite: 338]
    n_chi = chi_squared / (2.0 * epsilon);
    
    % 我们需要的是粒子数，所以向上取整
    n_chi = ceil(n_chi);
    
end