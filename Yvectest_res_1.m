function [A,W,alpha,obj,P] = Yvectest_res_1(X,gt,d, lambda,K,beta,avgS)%改
%Yvectest(X,Y,anchor(ichor),lambda,K,beta,avgS);
% X      : n*di

%% initialize
maxIter = 100 ; % the number of iterations
rng('default')
k = 20;%length(unique(gt));%20;
m = length(X);
n = size(gt,1);
P = rand(k, n);%自己加
sign=1;


A = cell(m,1); 
W = cell(m,1); 
for i = 1:m
   %di = size(X{i},2); 
   %A{i} = rand(di,d); % di * d
   A{i} =anchor_gen(X{i}',d,sign);
   W{i} = rand(d,k); % d * k
   X{i} = mapstd(X{i},0,1); % turn into dv*n
end

alpha = ones(1,m)/m;
opt.disp = 0;

Dk = sum(K, 2); %改
Lk= diag(Dk) - K; %改

Ds = diag(sum(avgS, 2));
Ls = Ds - avgS;


flag = 1; %judge if convergence
iter = 0;
obj = [];
%%
while flag
    iter = iter + 1;

    %% optimize A_i  dv*d
    parfor iv=1:m
        [Ua,~,Va] = svds(X{iv}*P'*W{iv}',d);
        A{iv} = Ua*Va';
    end
    

    %% optimize Y  k*n
%     loss1 = zeros(n, k);
%     for ij=1:m
%         loss1 = loss1 + alpha(ij) * EuDist2(X{ij}', W{ij}'*A{ij}', 0);
%     end
    %[~, Y] = min(loss, [], 2); %minvalue of row,minvalue_position
    
    %A1=alpha(1)*W{1}'* A{1}'* A{1} * W{1} +alpha(2)* W{2}'* A{2}'* A{2} * W{2};
    %B1=lambda*(beta * Ls + (1-beta)*Lk);
    %C1=(alpha(1)*W{1}'*A{1}'*X{1} + alpha(2)*W{2}'*A{2}'*X{2});
    

    A1 = zeros(k,k);
    C1 = zeros(k,n);
    parfor ip = 1:m
        A1 = A1 + alpha(ip) * W{ip}' * A{ip}' * A{ip} * W{ip};
        C1 = C1 + alpha(ip) * W{ip}' * A{ip}' * X{ip};
    end

    B1 = lambda * (beta * Ls + (1 - beta) * Lk);


    P = sylvester(A1,B1,C1);
    %[~ , Y1] = max(P);
    %Y=Y1';
    
    
        %% optimize W_i d*k
    %=====W_i>=0
        
        tmp = cell(m, 1); % 预分配tmp作为cell数组
        %W = cell(m, 1); % 预分配W
        pp = pinv(P * P'); % 使用更稳定的伪逆

        for iw = 1:m
            tmp{iw} = A{iw}' * X{iw} * P'; % 保持原逻辑，但确保A和X的索引正确
            W{iw} = tmp{iw} * pp; % 注意这里直接使用了tmp_1{i}
        end
    
    
    

    %% optimize alpha
    aloss = zeros(1,m);
    
    parfor iv = 1:m
        aloss(iv) = norm(X{iv}-A{iv}*W{iv}*P, 'fro');
        
        %aloss(iv) = sqrt(sum(sum((X{iv}-A{iv}*W{iv}(:,Y)).^2)));
        %aloss(iv) = sqrt(sum(sum((X{iv}-A{iv}*W{iv}*P).^2)));
    end
    alpha_o=1./(2*aloss);
    alpha = alpha_o/ sum(alpha_o);

    %%
    %loss
    J_LE1 = sum(diag(P * Ls * P')); 
    J_LE2 = sum(diag(P * Lk * P')); 
    term_3=beta * J_LE1 + (1-beta)*J_LE2;
    aloss_sq = aloss .^ 2;
    %term = zeros(m, 1);
    %for iv = 1:m
        %term(iv) = sum(sum((X{iv}-A{iv}*W{iv}(:,Y)).^2));
        %term(iv) = sum(sum((X{iv}-A{iv}*W{iv}*P).^2));
    %end
    %obj(iter) = alpha*term;
    
    obj(iter) = sum(alpha .* aloss_sq)+lambda*term_3;
    
    
    
    if (iter>2) && (abs((obj(iter)-obj(iter-1)))<1e-4 ||abs((obj(iter)-obj(iter-1))/(obj(iter-1)))<1e-4 || iter>maxIter ||obj(iter) < 1e-10)%
        flag = 0;
    end
end
         
         
    
