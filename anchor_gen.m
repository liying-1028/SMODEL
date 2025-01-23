% Construct similarity matrix with Approximate Nearest Neighbor Search (ANNS). 
% It is a parameter free, distance consistent similarity.
function Anchor = anchor_gen(X,numAnchor,sign)
% X: num*dim data matrix
% numAnchor: number of anchors
% sign: 0 for random selection; 1 for kmeans; 2 for Hierarchical kmeans. 

% Anchor: each column is a data point
% if nargin < 3
%     sign = 0;
% end

[num,dim] = size(X);

%% ------------------- 1. anchor generation ------------------

% if sign == 0 % random
%     sample_anchor = randperm(num,numAnchor);
%     Anchor = X(sample_anchor,:);
% elseif sign == 1 % kmeans
%     [StartIni,~] = InitializeG(num,numAnchor);
%     [~, ~, Anchor] = kmeans_ldj(X,StartIni);
% elseif sign == 2 % Hierachical kmeans
%     num_layer = round(log2(numAnchor));
%     Anchor = Hierachical_KM(X,num_layer);
%     X = X-repmat(mean(X),[num,1]);

%{
% BKHK  ¶þ²æÊ÷
[la,Anchor] = hKM(X',[1:num],numAnchor,1);
Anchor = Anchor';
% random
sample_anchor = randperm(num,2^numAnchor);
Anchor = X(sample_anchor,:);
% Kmeans
[la,Anchor] = kmeans(X,2^numAnchor);
%}

if sign == 0 % random
    %fprintf('Random Method \n');
    sample_anchor = randperm(num,numAnchor);
    Anchor = X(sample_anchor,:);
elseif sign == 1 % kmeans
    %fprintf('K-means Method \n');
%     [StartIni,~] = InitializeG(num,2^numAnchor);
%     [~, ~, Anchor] = kmeans_ldj(X,numAnchor);
%     [~,Anchor] = kmeans(X,);
    [~, Anchor]=litekmeans(X, numAnchor);
elseif sign == 2 %  BKHK 
    %fprintf('BKHK Method \n');
    [~,Anchor] = hKM(X',[1:num],numAnchor,1);
    Anchor = Anchor';
end


