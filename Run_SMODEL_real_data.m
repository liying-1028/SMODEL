close all; clear all; clc
warning off;
addpath(genpath('./'));
ds = {'A1'};
resPath = './A1/';
metric = {'ACC','nmi','Purity','Fscore','Precision','Recall','AR','Entropy','ASW'};
load('./data/A1.mat')
rng('default')
for dsi = 1:length(ds)
    % load data & make folder
    dataName = ds{dsi}; disp(dataName);
    k = length(unique(Y));

    matpath = strcat(resPath,dataName);
    txtpath = strcat(resPath,strcat(dataName,'.txt'));
    if (~exist(matpath,'file'))
        mkdir(matpath);
        addpath(genpath(matpath));
    end
    dlmwrite(txtpath, strcat('Dataset:',cellstr(dataName), '  Date:',datestr(now)),'-append','delimiter','','newline','pc');
    %% para setting
    anchor = 100;
    numAnchor=100;
    iters = 1;

    mu = 1e-3;
    gamma = 0.5;
    [rows_1, cols_1] = size(X{2});
    avgS = processData(poolsget, gamma,cols_1,rows_1,k);
   
    pairwise_distances = squareform(pdist(Spot, 'euclidean'));%euclidean
      for i = 1:size(pairwise_distances, 1)
       pairwise_distances(i,i) = 1;
      end

    if pairwise_distances ~= 0
        inverse_distance = 1 ./ pairwise_distances;
    else
        inverse_distance = Inf;
    end

    [graph_K,S_K] = sparse_similarity_matrices_descend(inverse_distance,12);%12

    row_sums = sum(S_K, 2);
    K = bsxfun(@rdivide,S_K, row_sums);
    K = SPPMI(K, 1);
    LD_1=X{2}';
    
    [LD_new,LD_new_gene] = processLD(LD_1, inverse_distance,12, 0.5);%knn
    [rows, cols] = size(LD_1);
    LD_end = zeros(rows, cols);
    
     for i = 1 : rows
        for j = 1 : cols
            if LD_1(i, j) == 0
                LD_end(i, j) = max(LD_new(i, j), LD_new_gene(i, j));
            else
                LD_end(i, j) = LD_1(i, j);
            end
        end
    end

    
    X{3}=LD_end';


    r1 =10;
    beta1 =0.2;
    
    %%
    for ichor = 1:length(anchor)
        for iter=1:iters

            for lam = 1 : length(r1)
                lambda = r1(lam);
                
                fprintf('Current lambda: %f\n', lambda);
                for beta_index = 1:length(beta1)
                    beta = beta1(beta_index);
                    fprintf('Current beta: %f\n', beta);
                
                    tic;
                    [A,W,alpha,obj,P] = SMODEL(X,Y,anchor(ichor),lambda,K,beta,avgS);
                    plot(obj);
                    xlabel('Iteration');
                    ylabel('Loss');
                    title('Loss Curve');

                    timer  = toc;
                    fprintf('\niter: %d, time: %.2f', iter, timer);

                    P1=P';
                    Ypre1=kmeans(P1,10);
                    res_1= Clustering8Measure(Y,Ypre1);
                    fprintf('\nAccuracy: %.4f, NMI_kmeans: %.4f, Purity: %.4f, F-score: %.4f, Precision: %.4f, Recall: %.4f, AR: %.4f, Entropy: %.4f\n', res_1(1), res_1(2), res_1(3), res_1(4), res_1(5), res_1(6), res_1(7), res_1(8));
                    
                    loss = min(obj);
                    
                    
                    
                    
                end
            end
                
           
        end
        
    end
   
   
end
