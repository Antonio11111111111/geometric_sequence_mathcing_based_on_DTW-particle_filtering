function denoised_data = Denoise_Visushrink(data, wavelet)
%DENOISE_VISUSHRINK_MATLAB 使用 Donoho 1995 的 VisuShrink 方法对信号降噪
%
%   输入:
%   data (vector):    输入的一维地磁序列 (被噪音污染的数据 d_i)
%   wavelet (string): 要使用的小波基函数, 例如 'db4', 'sym8'
%
%   输出:
%   denoised_data (vector): 降噪后的“干净”信号 (f_hat)

    % 确保输入是行向量或列向量
    if size(data, 1) > 1 && size(data, 2) > 1
        error('输入数据必须是一维向量');
    end
    
    % 原始信号长度
    n = length(data);
    
    % --- 步骤 1: 小波分解 ---
    % 自动计算最大可能的小波分解层数
    wlevel = wmaxlev(n, wavelet);
    % [c, l] 是小波分解的系数向量和长度向量
    [c, l] = wavedec(data, wlevel, wavelet);

    % --- 步骤 2: 计算“通用门槛” (Universal Threshold) ---
    
    % 2a. 估计噪音水平 sigma
    % l 的结构是 [len(cA), len(cD_L), ..., len(cD_1), total_len]
    % 我们需要 cD1 (最精细的细节系数)
    % cD1 的长度是 l(end-1)
    % cD1 的起始索引是 l(1)+l(2)+...+l(end-2) + 1
    cD1_len = l(end-1);
    cD1_start_index = sum(l(1:end-2)) + 1;
    cD1_end_index = cD1_start_index + cD1_len - 1;
    
    cD1 = c(cD1_start_index : cD1_end_index);
    
    % 使用中位数绝对差 (MAD) 估计 sigma
    % 对应论文中当 sigma 未知时的估计方法
    sigma_est = mad(cD1, 1) / 0.6745; % mad(X, 1) 计算关于中位数的MAD

    % 2b. 计算通用门槛 t_n
    % 对应论文中的 t_n = sigma * sqrt(2*log(n))
    % 在 MATLAB 中, log() 就是自然对数 ln()
    threshold = sigma_est * sqrt(2 * log(n));

    % --- 步骤 3: 应用软阈值 (Soft-Thresholding) ---
    
    % 我们只对“细节系数”进行阈值处理
    % 近似系数 cA (c(1:l(1))) 保持不变
    cA = c(1:l(1));
    
    % 提取所有细节系数
    % 它们从 cA 之后开始
    detail_coeffs = c(l(1)+1 : end);
    
    % 使用 wthresh 进行软阈值处理 ('s' 代表 'soft')
    % 这完美实现了论文中的 eta_t(y) 公式
    denoised_details = wthresh(detail_coeffs, 's', threshold);
    
    % 组合回完整的系数向量
    denoised_c = [cA, denoised_details];

    % --- 步骤 4: 小波重构 ---
    denoised_data = waverec(denoised_c, l, wavelet);

    % 确保输出长度与输入一致 (waverec 可能因边界效应产生 1 点的差异)
    denoised_data = denoised_data(1:n);
end