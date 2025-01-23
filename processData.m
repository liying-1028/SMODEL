function [E] = processData(poolsget, gamma,  n, q, kpre)
    % processData - 处理数据并生成邻接矩阵
    % 输入参数:
    %   poolsget - 输入数据矩阵
    %   mu - 收敛阈值
    %   gamma - 权重参数
    %   highorder - 高阶参数
    %   n - 数据维度
    %   q - 数据维度
    %   kpre - 预处理参数
    % 输出参数:
    %   E - 最终的邻接矩阵
    mu = 1e-3;
    highorder=0;
    
    s = cell(1, 6);
    S = cell(1, 6);
    w = cell(1, 6);
    
    [s{1}, S{1}] = WeightingCAMatrix(poolsget, gamma); % s是指示矩阵，S是邻接矩阵
    tau = floor(n * (2000 + extractdata(relu(dlarray(q - 15000, "SSCB")))) / (3600 * (kpre + 1)));
    
    for i = 1:100
        w{i} = GeneratingWeights(poolsget, S{i}, tau, highorder);
        [s{i+1}, S{i+1}] = WeightingCAMatrix(poolsget, gamma, s{i}, w{i});
        if norm(S{i+1} - S{i}, 'fro') / norm(S{i}, 'fro') < mu
            break;
        end
    end
    
    S(cellfun(@isempty, S)) = []; % 移除空值
    E = cell2mat(S(end));
end
