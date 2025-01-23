function [LD_new,LD_new_gene,y_d,y_l] = processLD(LD, D, K, a)
    % processLD - 处理LD矩阵并生成新的关联矩阵
    % LD: 输入的LD矩阵
    % D: 细胞的相似矩阵
    % K: KNN的参数K
    % a: 权重参数a

    % 计算LD的余弦相似度矩阵
    L = pdist2(LD', LD', 'cosine');
    L = L - diag(diag(L));

    % 计算KNN矩阵
    [rL, cL] = size(L);
    KNN_L = zeros(rL, cL);
    [sort_L, idx] = sort(L, 2, 'descend');
    for i = 1 : rL
        KNN_L(i, idx(i, 1:K)) = sort_L(i, 1:K);
    end

    
    D = D - diag(diag(D));
    [rD, cD] = size(D);
    KNN_D = zeros(rD, cD);
    [sort_D, idx] = sort(D, 2, 'descend');%descend
    for i = 1 : rD
        KNN_D(i, idx(i, 1:K)) = sort_D(i, 1:K);
    end

    % 计算新的权重矩阵
    [rows, cols] = size(LD);
    y_l = zeros(rows, cols);
    y_d = zeros(rows, cols);

    knn_network_l = KNN_L;
    for i = 1 : cols
        w = zeros(1, K);
        [sort_l, idx_l] = sort(knn_network_l(i, :), 2, 'descend');
        sum_l = sum(sort_l(1, 1:K));
        for j = 1 : K
            w(1, j) = a^(j-1) * sort_l(1, j);
            y_l(:, i) = y_l(:, i) + w(1, j) * LD(:, idx_l(1, j));%y_l(i, :) = y_l(i, :) + w(1, j) * LD(idx_l(1, j), :);
        end
        y_l(:, i) = y_l(:, i) / sum_l;%y_l(i, :) = y_l(i, :) / sum_l
    end

    knn_network_d = KNN_D;
    for i = 1 : rows
        w1 = zeros(1, K);
        [sort_d, idx_d] = sort(knn_network_d(i, :), 2, 'descend');%descend
        sum_d = sum(sort_d(1, 1:K));
        for j = 1 : K
            w1(1, j) = a^(j-1) * sort_d(1, j);
            y_d(i, :) = y_d(i, :) + w1(1, j) * LD(idx_d(1, j), :);%第i列y_d(:, i) = y_d(:, i) + w1(1, j) * LD(:, idx_d(1, j));
        end
        y_d(i, :) = y_d(i, :) / sum_d;%y_d(:, i) = y_d(:, i)
    end

    % 计算新的关联矩阵
    LD_1 = zeros(rows, cols);
    for i = 1 : rows
        for j = 1 : cols
            LD_1(i, j) = max(LD(i, j), y_d(i, j));
        end
    end

    LD_new = LD_1;
    
    LD_gene = zeros(rows, cols);
    for i = 1 : rows
        for j = 1 : cols
           LD_gene(i, j) = max(LD(i, j), y_l(i, j));
        end
    end

    LD_new_gene = LD_gene;
    
end
