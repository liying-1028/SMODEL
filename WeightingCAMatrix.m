function [WA,S]=WeightingCAMatrix(poolsget,gamma,A,W)
[n,m]=size(poolsget);
S=zeros(n);
WA=cell(1,m);
%%The weight matrix elements are all 1 by default
if exist('W','var')==0
    W=ones(n,m);
end
%%Generate an initial co-occurrence matrix for each base clustering
if exist('A','var')==0
    for i=1:m
        v=poolsget(:,i);
        s{i}=zeros(n,max(v));
             for j=1:n
                 for k=1:max(v)
                     if v(j)==k
                        s{i}(j,k)=1;
                     end
                 end
             end
        A{i}=s{i}*s{i}';
    end
end
%%The weighted co-association matrix is generated
for i=1:m
    WA{i}=((W(:,i)*W(:,i)').^gamma).*A{i};
    S=S+WA{i};
end
S=S/m;
%%The cell and itself must belong to the same cluster
S(logical(eye(size(S))))=1;