function [dist] = simple_global_dtw(x, y)
    % 这是一个标准 (Global) DTW 函数，用于计算两个序列的距离
    m = length(x);
    n = length(y);
    
    D = zeros(m, n);
    local_cost = (x' - y).^2; 

    % 初始化
    D(1, 1) = local_cost(1, 1);
    for i = 2:m
        D(i, 1) = local_cost(i, 1) + D(i-1, 1);
    end
    for j = 2:n
        D(1, j) = local_cost(1, j) + D(1, j-1);
    end

    % 填充矩阵
    for i = 2:m
        for j = 2:n
            D(i, j) = local_cost(i, j) + min([D(i-1, j), D(i, j-1), D(i-1, j-1)]);
        end
    end

    % 最终距离
    dist = D(m, n);
end