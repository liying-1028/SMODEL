function [graph_K,S_K] = sparse_similarity_matrices_descend(pairwise_distances, P1)

    KK = pairwise_distances;
    K=KK-diag(diag(KK));
    [r1,c1]=size(K);
    PNN_K = zeros(r1, c1);
    graph_K = zeros(r1, c1);
    [sort_K,idx]=sort(K,2,'descend');%[sort_K,idx]=sort(K,2,'descend');
    
     for i = 1 : r1
        PNN_K(i,idx(i,1:P1))=sort_K(i,1:P1);
    end    
     for i = 1 : r1
        idx_i=find(PNN_K(i,:));
        for j = 1 : r1           
            idx_j=find(PNN_K(j,:));
            if ismember(j,idx_i) & ismember(i,idx_j) %&& isequal(Clu_mat(i,j),1)               
                graph_K(i,j)=1;
            elseif ~ismember(j,idx_i) & ~ismember(i,idx_j) %&& ~isequal(Clu_mat(i,j),1)  
                graph_K(i,j)=0;
            else
                graph_K(i,j)=0.5;               
            end       
        end
     end

    % 稀疏化相似矩阵
    S_K = KK .* graph_K;
end
